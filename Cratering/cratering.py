import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import qLabMods
from pint import UnitRegistry

ureg = UnitRegistry()
## logD = m*logE + logD0 - mlogE0
m =  65.31 # g
m = m/1000 # convert to kg
dm = 0.01 # g
dm = dm/1000
g = 9.8 # m/s**2

## Data taken
h = np.array([90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,
            165]) # cm
dh = 1 # cm
D = np.array([71.79,73.52,75.33,75.59,79.65,82.80,83.08,77.47,81.63,89.25,92.69,
            93.42,94.70,94.24, 98.14, 97.84]) # mm
D = D/10 # convert to cm
dD = .25 # cm

# Convert to m
h = np.array(h)/100 # m
dh = dh*np.ones(len(h))/100 # m
D = np.array(D)/100 # m
dD = dD*np.ones(len(D))/100 # m

## Calculate mgh for each value of h
E = []
for hval in h:
    E.append(m*g*hval)

## Uncertainty of E is calculated by adding uncertainties in quadrature
dE = E*np.sqrt((dm/m)**2 + (dh/h)**2)

# Natural log of Values
logE = np.log(E)
logD = np.log(D)
dlogD = (1/D[0])*dD

## Weighted fit - the good one
def func(x,A,B):
    return A*x+B

popt,pcov = curve_fit(func,logE,logD,sigma = dlogD)
perr = np.sqrt(np.diag(pcov))
slopeErr = perr[0]
p_weight = np.poly1d(popt)
yfit = func(logE, *popt)

plt.figure(1)
plt.plot(logE,logD,'k.',markersize = 15)
plt.errorbar(logE,logD,yerr=dlogD,fmt='k.',capsize = 4)
plt.plot(logE,yfit,'k--')
plt.xlabel('logE [log(J)]')
plt.ylabel('logD [log(m)]')
plt.savefig('cratering1.pdf')
plt.show()

## Method 2
m = np.log(D/D[0]) / np.log(E/E[0])
m = m[1:] # get rid of NaN

mstd = np.std(m)
sdom = mstd/np.sqrt(len(D))
mavg = np.mean(m)
counts = np.linspace(1,len(m),len(m))

plt.figure(2)
plt.axhline(y=mavg+sdom, color='k', linestyle='-')
plt.axhline(y=mavg, color='k', linestyle='--')
plt.axhline(y=mavg-sdom, color='k', linestyle='-')
plt.errorbar(counts,m,yerr = mstd,fmt = '.',capsize = 4)
plt.xlabel('Sample Number')
plt.ylabel('m value')
plt.savefig('cratering2.pdf')
plt.show()

# Find the energy of the crater
craterD = 1200 # m
# Energy found with our m value
craterEnergy = np.exp((1/np.mean(m))*np.log(craterD/D[0])+np.log(E[0]))
# What it should be
craterEnergy0 = np.exp((1/0.25)*np.log(craterD/D[0])+np.log(E[0]))

print('\nh array:',h)
print('\ndD array',dD)
print('\ndh array',dh)
print('\nD array',D)
print('\nEnergy array:',E)
print('\ndE array',dE)
print('\nlogD',logD)
print('dlogD',dlogD)
print('method 1 m value:',popt[0])
print('intercept',intercept)
print('method 2 m value',mavg)
print('standard deviation of mean',sdom)
print('Energy with our m:',craterEnergy)
print('Literature Energy:',craterEnergy0)

table = qLabMods.LaTeXTable([u'$Height\ \pm\ \delta h$','$Diameter\ \pm\ \delta D$',
                    '$Energy\ \pm\ \delta E$'],
                    [(h,dh,ureg.meter),(D,dD,ureg.meter),(E,dE,ureg.joule)],
                    'Height,Diameter, and corresponding Energy',
                    'DataTable', '!htp')
print(table.allTogether())
