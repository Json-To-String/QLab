import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


# Reading in the data
with open('speedofLight.csv','r') as sl:
    header = sl.readline()
    info = header.split(',')

    line = sl.readline()
    distList = []
    tList = []
    dt = []
    dD = []
    while line:
        line = line.split(',')
        distList.append(float(line[1]))
        tList.append(float(line[2]))
        dt.append(float(line[3]))
        dD.append(float(line[4]))
        line = sl.readline()

# Getting the data into vars
distList = 2*np.array(distList)*10**-2
tList = np.array(tList)*10**-9
dt = (np.array(dt)*10**-9)
dD = (np.array(dD))*10**-2

# The fit
def func(x,A,B):
    return A*x+B

popt,pcov = curve_fit(func,tList,distList,sigma = dt)
perr = np.sqrt(np.diag(pcov))
slopeErr = perr[0]
p_weight = np.poly1d(popt)
xtheory = np.linspace(0,5e-8)
yfit = func(xtheory,popt[0],popt[1])

# Values from line fit
c = np.format_float_scientific(popt[0],precision=3)
print('c:',c)
c = float(c)
cerr = c*np.sqrt((dD/distList)**2 + (dt/tList)**2)

# SDOM Analysis
v = distList/tList
sdom = np.std(v)/np.sqrt(len(v))
vBar = np.mean(v)
sample = np.arange(1,len(v)+1)

# Console Output
print('cerr:',cerr)
print('SDOM:',np.format_float_scientific(sdom))

# The Plots
# plt.figure(figsize = (12,9))
# plt.plot(tList,distList,'k.',markersize = 12)
# plt.errorbar(tList,distList,yerr=dD,xerr = dt,fmt='k.',capsize = 3)
# plt.plot(xtheory,yfit,'b--')
# plt.xlabel('time[s]')
# plt.ylabel('distance[m]')
# plt.savefig('SpeedOfLight.pdf')
# plt.show()

# plt.figure(figsize = (12,9))
# plt.errorbar(sample,v,yerr = cerr,fmt = '.',markersize = 10,capsize = 10)
# plt.axhline(vBar,ls = '--')
# plt.axhline(vBar + sdom)
# plt.axhline(vBar - sdom)
# plt.xlabel('Sample Number')
# plt.ylabel('Speed of Light[m/s]')
# plt.savefig('SpeedOfLightSDOM.pdf')
# plt.show()

# Pandas to LaTeX
dataDict = {'Distance': distList, 'Times': tList}
df = pd.DataFrame(data = dataDict)
print(df.to_latex())
