
import qLab as ql
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import pandas as pd

from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from pathlib import Path
import re

font = {'size' : 22}
matplotlib.rc('font', **font)

b = 3e30
r_0 = 2.82e-15 # m
ga = 662e3 / ql.massElectron 


matplotlib.rc('xtick', labelsize = 20)
matplotlib.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 25)

## Convert HH:MM:SS to seconds
def get_sec(time_str):
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + int(s)

def scatterE(E0, theta):
    # Takes theta in radians
    # Takes E0 in keV
    E0 = E0 * 1000
    den = 1 + ((E0 / ql.massElectron) * (1 - np.cos(theta)))
    return(E0 / (den * 1000))

def kleinNishina(th):
    A = b * (r_0 **2)
    B = (1 + (np.cos(th) ** 2)) / 2
    C = 1 / ( ( 1 + ga * (1 - np.cos(th)) ) ** 2 )
    DNum = (ga ** 2) * ((1 - np.cos(th))**2)
    DDen = (1 + (np.cos(th)**2)) * (1 + ga * (1 - np.cos(th)))
    D = 1 + (DNum / DDen)
    return(A * B * C * D)

def thompson(th):
    A = b * (r_0 **2)
    B = (1 + (np.cos(th) ** 2)) / 2
    return(A * B)

def chiSqr(exp, the, unc):
    chi = []
    for exp1, the1, unc1 in zip(exp, the, unc): 
        if unc1 != 0:
            numerator = (exp1 - the1)**2
            denominator = unc1 ** 2 
            chi.append(numerator / denominator)

    return(np.sum(chi))

df = ql.getCsv('https://docs.google.com/spreadsheets/d/1omtMNUBSYKocbhReMvC2iGUG-THQABjCvus0ghR2J80/edit#gid=0')
col = lambda x : ql.colVar(df, x)
sCol = lambda x : ql.stringVar(df, x)

## Calibration Curve ##
channel = col(0)
channelUnc = col(2)
energy = col(1)
calOpt, calErr, calThe, calFit = ql.fit_routine(ql.quadratic_fit,
        energy, channel, channelUnc)

fig = plt.figure(figsize = (12, 9))
ax = fig.gca()
ax.tick_params(axis = 'both', direction = 'in', 
                top = True, right = True, 
                pad = 25)
ax.plot(calThe, calFit, 'r-')
ax.errorbar(energy, channel, yerr = channelUnc, fmt = 'k.', ms = 10)
ax.set_title("Energy vs Channels Calibration Curve", fontsize = 20)
ax.set_xlabel("Energy [keV]", fontsize = 20)
ax.set_ylabel("Channel Number", fontsize = 20)
#plt.show()
plt.savefig('calibration.png')

toEnergy = lambda channel : ql.quadratic_inv(channel, *calOpt)
toChannel = lambda energy : ql.quadratic_fit(energy, *calOpt)
## Calibration Curve ##

#deg100 = scatterE(662, ql.toRadians(100))
#print(deg100)
#print(toChannel(deg100))

## Compare ##
angles = col(5)
counts = col(9)
rate = counts / np.array([get_sec(y) - get_sec(x) for x, y in zip(col(11), col(12))]) 
rateUnc = np.sqrt(rate)

xTest = np.linspace(0, 180)
xRad = ql.toRadians(xTest)
yKlein = kleinNishina(xRad)
yThomp = thompson(xRad)

fig = plt.figure( figsize = (12, 9))
ax = fig.gca()
ax.tick_params(axis = 'both', direction = 'in', 
                top = True, right = True, 
                pad = 25)

ax.plot(xTest, yKlein, 'r-')
ax.plot(xTest, yThomp, 'b-')
ax.errorbar(angles, rate, yerr = rateUnc, fmt = 'k.', capsize = 3)

ax.set_title('Klein-Nishina vs. Thompson')
ax.set_xlabel('Angle [degrees]')
ax.set_ylabel('Rate [Counts / sec]')
#plt.show()
plt.savefig('comparison.png')


