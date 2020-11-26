import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from pint import UnitRegistry
import sys
sys.path.append('../PYTHON_LIBRARIES/')
import qLabTools as ql
ureg = UnitRegistry()

e = 1.602176634e-19
c = 2.99792458e8

## IV Curve
with open('photoelectricIV.csv','r') as iv:
    ivHeader = iv.readline()
    V = []
    dV = []
    I = []
    dI = []
    line = iv.readline()
    while line:
        line = line.split(',')
        V.append(float(line[0]))
        dV.append(float(line[1]))
        I.append(float(line[2]))
        dI.append(float(line[3]))
        line = iv.readline()

## Wavelengths and stopping voltages
with open('photoelectric.csv','r') as pe:
    peHeader = pe.readline()
    wave = []
    halfWave = []
    dWave = []
    vStop = []
    dVStop = []
    line = pe.readline()
    while line:
        line = line.split(',')
        wave.append(float(line[0]))
        halfWave.append(float(line[1]))
        dWave.append(float(line[2]))
        vStop.append(float(line[3]))
        dVStop.append(float(line[4]))
        line = pe.readline()

## Cast lists to arrays
wave = np.array(wave)
dWave = np.array(dWave)
halfWave = np.array(halfWave)
vStop = np.array(vStop)
dVStop = np.array(dVStop)


deltaW = wave - halfWave
doubleDelta = halfWave - deltaW
squareDelta = wave - deltaW*np.sqrt(5)
# xData = [wave,halfWave,doubleDelta,squareDelta]
xData = [wave]
# print(wave)
# print(deltaW)
# print(doubleDelta)

def func(x,A,B):
    return(A*x + B)
def plots(x,y,dVx,dVy,fitFunction):

    popt,pcov = curve_fit(func,x,y,sigma = dVx)
    perr = np.sqrt(np.diag(pcov))
    slopeErr = perr[0]
    p_weight = np.poly1d(popt)
    xtheory = np.linspace(min(x),max(x))
    yfit = fitFunction(xtheory,popt[0],popt[1])

    plt.errorbar(x,y,yerr = dVy, fmt ='b.',capsize = 4)
    plt.plot(xtheory,yfit,'k--')

    print(popt[0])

plt.figure()
plt.plot(V,I,'k.')
plt.xlabel('Voltage [V]')
plt.ylabel('Current [nA]')
plt.savefig('IVcurve.png')

plt.figure()
for xVals in xData:
    freq = c/(xVals*10**-9)
    plots(freq,vStop,dWave,dVStop,func)
plt.xlabel('Frequencies [Hz]')
plt.ylabel('Stopping Voltage [eV]')
plt.savefig('vStopvsFreqs.png')

freq = c/(halfWave*1e-9)
freqUnc = (dWave/halfWave)

def wFit(x,h,w):
	return((h*x) - w)

latexTemp, latexRets, numRets = ql.fitPlotter(wFit, freq, vStop, [ureg.hertz, ureg.eV], ["Frequency", "Stopping Voltage"], uncX = freqUnc, uncY = dVStop)
print(latexTemp)

def valFunction(wavelength, stoppingVoltage):
	return((wavelength/c)*(stoppingVoltage + numRets[1][0]))

def uncFunction(wavelength, stoppingVoltage, deltaWave, deltaVolt):
	h = valFunction(wavelength, stoppingVoltage)
	deltaH = h*np.sqrt((deltaWave/wavelength)**2 + (deltaVolt/stoppingVoltage)**2)
	return(deltaH)

latexTemp1, latexRets1, numRets1 = ql.sdomPlotter(valFunction, uncFunction, ureg.eV, "Planck's Constant", [halfWave*1e-9, vStop], [dWave*1e-9, dVStop])
print(latexTemp1)

print(ql.latexTable(["Wavelength", "Stopping Voltage"], [halfWave*1e-9, vStop], [dWave*1e-9, dVStop], [ureg.meter, ureg.eV]))
print(ql.latexTable(["Planck's Constant (h)", "Work Function"], [[numRets[0][0],numRets1[0]], [numRets[1][0]]], [[numRets[0][1], numRets1[1]], [numRets[1][1]]], [ureg.eV, ureg.eV]))

print(numRets[1])
