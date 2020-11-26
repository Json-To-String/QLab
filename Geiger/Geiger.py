import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from scipy.special import factorial
from matplotlib.ticker import FormatStrFormatter

deadTime = 606e-6
pulseWidth = 32.8e-6
recoveryTime = 1.41e-3
hiRate = 10
loRate = .1
bgRate = 100

## Correction
def nTrue(nObs,rObs):
	nObs = np.array(nObs)
	C = 1/(1-rObs*deadTime)
	return C*nObs

## Read data from file
with open('GeigerHighLow.csv','r') as gc:
	header = gc.readline()
	line = gc.readline()
	high = []
	low = []
	while line:
		line = line.split(',')
		high.append(float(line[0]))
		low.append(float(line[1]))
		line = gc.readline()

high_corr = nTrue(high,hiRate)
low_corr = nTrue(low,loRate)
high_sorted = np.sort(high_corr)
low_sorted = np.sort(low)

def gauss(x,xbar):
	sig = np.sqrt(xbar)
	g = (1/(np.sqrt(2*np.pi)*sig))*np.exp((-(x-xbar)**2)/(2*sig**2))
	return(g)

gaussian = lambda x,xbar : ((1/(np.sqrt(2*np.pi)*np.sqrt(xbar)))
							* np.exp((-(x-xbar)**2)/(2*np.sqrt(xbar)**2)))

def poisson(x,xbar):
	A = np.exp(-xbar)
	B = []
	for xVal in x:
		xVal = int(xVal)
		B.append((xbar**xVal)/math.factorial(xVal))
	B = np.array(B)
	return(A*B)

# Start fig
fig1,ax1 = plt.subplots(figsize = (12,9))

# Histogram and errorbar
n_bin1 = np.arange(high_sorted[0],high_sorted[-1]+15,15)
binspace1 = n_bin1[2] - n_bin1[1]
mid1 = n_bin1 + (binspace1/2)
mid1 = mid1[:-1]
hi_counts,hi_entries,hi_patches = ax1.hist(high_corr,bins = n_bin1,align = 'mid')
ax1.errorbar(mid1,hi_counts,yerr=np.sqrt(hi_counts),capsize = 3,fmt = 'k.')

# Gaussian Fit
xtheory1 = np.linspace(high_sorted[0],high_sorted[-1])
gaussCorr = len(mid1)*len(high)*gauss(xtheory1,np.mean(high_corr))
ax1.plot(xtheory1,gaussCorr,'r--')

# Chi Squared
chiG = 0
for index, count in enumerate(mid1):
	g = (len(mid1)*len(high)*gaussian(count,np.mean(high_corr)))
	chiG = chiG + (((hi_counts[index] - g)**2)/hi_counts[index])

# Plot formatting
ax1.set_xlabel('Counts',fontsize = 14)
ax1.set_ylabel('Frequency of Counts',fontsize = 14)
ax1.set_xticks(mid1)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
ax1.tick_params(labelsize = 12)

fig1.savefig('GaussHist.png')
# plt.show()
plt.close()

# Start fig
fig2,ax2 = plt.subplots(figsize = (12,9))

# Histogram
n_bin2 = np.arange(low_sorted[0],low_sorted[-1])
binspace2 = n_bin2[1] - n_bin2[0]
mid2 = n_bin2 + binspace2/2

lo_counts,lo_entries,lo_patches = ax2.hist(low,bins = np.arange(-0.5,7.5),range = (-.05,7.5),align = 'mid')
ax2.errorbar(mid2[:-1]-.5,lo_counts,yerr = np.sqrt(lo_counts),capsize = 2,fmt ='k.')
ax2.errorbar([7,8,9],[0,0,0],yerr = 3,capsize = 2,fmt = 'k.')

# Poisson Fit
xtheory2 = np.linspace(0,10,11)
poissonCorr = 100*1*poisson(xtheory2,np.mean(low))
ax2.plot(xtheory2,poissonCorr,'ro')

# Chi Squared
chiP = 0
for index,count in enumerate(mid2[:-1]):
	p = (100*1*poisson([count],np.mean(low)))
	if lo_counts[index] >= 5:
		chiP = chiP + (((lo_counts[index] - p)**2)/lo_counts[index])
chiP = float(chiP[0])


# ax2.set_xticks(mid2)
ax2.set_xlabel('Counts',fontsize = 14)
ax2.set_ylabel('Frequency of Counts',fontsize = 14)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
ax2.tick_params(labelsize = 12)
ax2.set_ylim(0,30)
plt.savefig('PoissonHist.png')
# plt.show()
plt.close()

hiBar = np.mean(high_corr)
hiSig = np.std(high_corr)
hi_sigTheory = np.sqrt(hiBar)
hiBar_unc = hiSig/np.sqrt(len(high_corr))
hi_rel = hiSig/hi_sigTheory

loBar = np.mean(low_corr)
loSig = np.std(low_corr)
lo_sigTheory = np.sqrt(loBar)
loBar_unc = loSig/np.sqrt(len(low_corr))
lo_rel = loSig/lo_sigTheory


## Comparison
hi_oneSig = len([x for x in high_corr
if hiBar - hiSig <= x <= hiBar + hiSig])/len(high_corr)
hi_twoSig = len([x for x in high_corr
if hiBar - 2*hiSig <= x <= hiBar + 2*hiSig])/len(high_corr)

lo_oneSig = len([x for x in low_corr
if loBar - loSig <= x <= loBar + loSig])/len(low_corr)
lo_twoSig = len([x for x in low_corr
if loBar - 2*loSig <= x <= loBar + 2*loSig])/len(low_corr)


print('hi counts var: ', hi_counts)
print('\nGaussian chi-Squared:',chiG)
print(f'N Bar High: \t {hiBar}')
print(f'{hi_oneSig} of the high values are within one sigma')
print(f'{hi_twoSig} of the high values are within two sigma')
print(f'High Sigma Exp: \t{hiSig}')
print(f'High Sigma Theory: \t{hi_sigTheory}')
print(f'High reliability factor: {hi_rel}')


print('\nPoissonian Chi-Squared:',chiP)
print(f'N Bar Low: \t {loBar}')
print(f'{lo_oneSig} of the low values are within one sigma')
print(f'{lo_twoSig} of the low values are within two sigma')
print(f'Low Sigma Exp: \t\t{loSig}')
print(f'Low Sigma Theory: \t{lo_sigTheory}')
print(f'Low reliability factor: {lo_rel}\n')

df = pd.DataFrame({
	'NBar' : [hiBar,loBar],
	'NBarUnc' : [hiBar_unc,loBar_unc],
	'sigTheory' : [hi_sigTheory,lo_sigTheory],
	'sigExp' : [hiSig,loSig],
	'ChiSq': [chiG,chiP],
	'1sig': [hi_oneSig,lo_oneSig],
	'2sig': [hi_twoSig,lo_twoSig],
	'Reliability' : [hi_rel,lo_rel],
	})

print(df.to_latex())
