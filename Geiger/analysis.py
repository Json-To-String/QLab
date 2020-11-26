#Analysis for the Geiger Counting Lab
#Libraries Required
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pint import UnitRegistry
ureg = UnitRegistry()

import sys
sys.path.append("../PYTHON_LIBRARIES/")
import qLabTools as ql

matplotlib.rcParams['font.sans-serif'] = "Liberation Serif"
matplotlib.rcParams['font.family'] = "Liberation Serif"

#### DATA ####

pulseWidth = 57e-6
deadTime = 650e-6
recoveryTime = 2.5e-3

# Per 100 Seconds
backgroundCount = 22 

# Per 10 Seconds
highCount = np.array([2090, 2176, 2162, 2145, 2149, 2145, 2116, 2160, 2134, 2173, 2101, 2108, 2132, 2183, 2212, 2157, 2185, 2143, 2102, 2165, 2113, 2106, 2099, 2166, 2106, 2162, 2123, 2122, 2138, 2062, 2125, 2129, 2072, 2132, 2195, 2178, 2135, 2096, 2185, 2161, 2142, 2170, 2142, 2167, 2164, 2159, 2190, 2169, 2133, 2167, 2164, 2158, 2202, 2129, 2077, 2107, 2145, 2178, 2175, 2124, 2168, 2095, 2202, 2190, 2132, 2183, 2141, 2107, 2147, 2199, 2045, 2140, 2205, 2123, 2175, 2169, 2123, 2164, 2112, 2137, 2125, 2123, 2137, 2152, 2142, 2140, 2208, 2150, 2145, 2119, 2127, 2131, 2144, 2136, 2203, 2068, 2179, 2171, 2173, 2168])
highRate = highCount*.1
cH = 1/(1-(highRate*deadTime))

# Per .1 Seconds
lowCount = np.array([1, 4, 3, 2, 2, 6, 1, 6, 6, 4, 3, 2, 2, 5, 2, 5, 1, 1, 2, 4, 1, 5, 1, 6, 2, 1, 2, 3, 4, 0, 3, 2, 1, 0, 1, 2, 3, 5, 2, 4, 5, 3, 4, 2, 4, 0, 2, 3, 1, 3, 6, 6, 4, 5, 4, 3, 0, 0, 3, 3, 3, 4, 4, 3, 5, 6, 1, 6, 4, 3, 4, 3, 0, 5, 2, 1, 3, 9, 3, 3, 2, 3, 4, 1, 1, 3, 1, 5, 1, 2, 4, 1, 3, 4, 3, 5, 3, 2, 4, 1])
lowRate = lowCount*10
cL = 1/(1-(lowRate*deadTime))

highTrue = highCount*cH 
lowTrue = lowCount*cL 

#### DATA ####

def gaussian(x, xBar):
    normConst = (1/np.sqrt(2*np.pi*xBar))
    guassExp = np.exp((-(x - xBar)**2)/(2*xBar)) 
    return(normConst * guassExp)

def poisson(x, xBar):
    expP = np.array(np.exp(-xBar))
    powerP = np.array(np.power(xBar, x))
    gammaP = np.array([math.gamma(count + 1) for count in x])
    return((expP*powerP)/gammaP)

def basicStats(true):
    nBar = np.mean(true) 
    sigmaExp = np.std(true) 
    sigmaThe = np.sqrt(nBar) 
    rel = sigmaExp/sigmaThe
    return(nBar, sigmaExp, sigmaThe, rel)

(nBarHigh, sigmaHighExp, sigmaHighThe, relHigh) = basicStats(highTrue)
(nBarLow, sigmaLowExp, sigmaLowThe, relLow) = basicStats(lowTrue)

#### PLOTTING ####

#### HIGH COUNTS ####
figure1, axis1 = plt.subplots(figsize = (12,9))

highBinSpacing = 15
highBins = np.arange(min(highTrue), max(highTrue) + highBinSpacing, highBinSpacing)
highMiddleBins = (highBins + (highBinSpacing/2))[:-1]
highCounts, highEntries, highPatches = axis1.hist(highTrue, bins = highBins)
axis1.errorbar(highMiddleBins, highCounts, yerr=np.sqrt(highCounts),capsize = 3, fmt = 'k.')

xTheoryHigh = np.linspace(min(highTrue), max(highTrue))
gaussianScaling = highBinSpacing * len(highTrue)
gaussianTheory = gaussianScaling * gaussian(xTheoryHigh, nBarHigh) 
axis1.plot(xTheoryHigh, gaussianTheory,'r--')

axis1.set_xlabel('Counts', fontsize=30)
axis1.set_ylabel('Frequency of Counts', fontsize=30)
axis1.tick_params(labelsize = 20)
figure1.savefig('GaussHist.pdf')

#### LOW COUNTS ####
figure2, axis2 = plt.subplots(figsize = (12,9))

lowBinSpacing = 1
lowBins = np.arange(min(lowTrue), max(lowTrue) + lowBinSpacing, lowBinSpacing)
lowMiddleBins = (lowBins + (lowBinSpacing/2))[:-1]
lowCounts, lowEntries, lowPatches = axis2.hist(lowTrue, bins = lowBins)
axis2.errorbar(lowMiddleBins, lowCounts, yerr=np.sqrt(lowCounts),capsize = 3, fmt = 'k.')

xTheoryLow = np.linspace(min(lowTrue), max(lowTrue))
poissonScaling = lowBinSpacing * len(lowTrue)
poissonTheory = poissonScaling * poisson(xTheoryLow, nBarLow) 
axis2.plot(xTheoryLow, poissonTheory,'r--')

axis2.set_xlabel('Counts', fontsize=30)
axis2.set_ylabel('Frequency of Counts', fontsize=30)
axis2.tick_params(labelsize = 20)
figure2.savefig('PoissHist.pdf')

#### CHI SQUARED ANALYSIS ####
chiHigh = 0
for index, count in enumerate(highMiddleBins):
    gaussianValue = gaussianScaling * gaussian(count, nBarHigh)
    if highCounts[index] >= 5:
        chiHigh = chiHigh + (((highCounts[index] - gaussianValue)**2)/highCounts[index])

chiLow = 0
for index, count in enumerate(lowMiddleBins):
    poissonValue = poissonScaling * poisson([count], nBarLow)
    if lowCounts[index] >= 5:
        chiLow = chiLow + (((lowCounts[index] - poissonValue)**2)/lowCounts[index])

print('Gaussian Chi-Squared:', chiHigh)
print('Poisson Chi-Squared:', chiLow[0])

print('High Counts')
oneSigmaHigh = len([x for x in highTrue if nBarHigh - sigmaHighExp <= x <= nBarHigh + sigmaHighExp])/len(highTrue)
twoSigmaHigh = len([x for x in highTrue if nBarHigh - 2*sigmaHighExp <= x <= nBarHigh + 2*sigmaHighExp])/len(highTrue)

print('Low Counts')
oneSigmaLow = len([x for x in lowTrue if nBarLow - sigmaLowExp <= x <= nBarLow + sigmaLowExp])/len(lowTrue)
twoSigmaLow = len([x for x in lowTrue if nBarLow - 2*sigmaLowExp <= x <= nBarLow + 2*sigmaLowExp])/len(lowTrue)
