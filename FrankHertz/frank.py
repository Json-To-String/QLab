import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import qLab as ql
import pandas as pd
import os
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import sys

matplotlib.rc('xtick', labelsize = 20)
matplotlib.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 20)

hgFirst = 4.67 # eV

def quartic(x, A, B, C, D, E):
    return(A*x**4 + B*x**3 + C*x**2 + D*x + E)

def routine(filename, graphind):
    fitType = ql.quadratic_fit
    df0 = pd.read_csv(filename, index_col = None)
    x, y = ql.colVar(df0, 0), ql.colVar(df0, 1)*1e12

    x = x[5:]
    y = y[5:]

    popt, perr, xtheory, yfit = ql.fit_routine(fitType, x, y)
    perr = perr / 100
    yUncUp = np.abs(yfit + fitType(xtheory, *(popt - perr)))
    yUncDown = np.abs(yfit - fitType(xtheory, *(popt - perr)))

    yUp = yfit + yUncDown
    yDown = yfit - yUncDown
    minY = np.where(yfit == min(yfit))[0][0]
    minX = x[minY]

    ySpread = np.mean(np.abs(ql.quadratic_fit(x, *popt) - y))

    print(f'Minima: {minX} +/- {ySpread:.1f}')

    plt.figure(figsize = (12, 9))
    plt.plot(x, y, 'k.')
    if graphind != 0:
        plt.plot(xtheory, yfit, 'b-')
        plt.plot(xtheory, yUp, 'r-')
        plt.plot(xtheory, yDown, 'r-')
    plt.xlabel('Accelerating Voltage [V]', fontsize = 20)
    plt.ylabel('Current [pA]', fontsize = 20)
    plt.savefig(f'scan{graphind}.png')
    plt.close()

    return minX, ySpread

path = os.path.join('Data', '*')
files = sorted(glob.glob(path))
vActual = []
minimaUnc = []
for ind, csv in enumerate(files):
    print(f'\n{csv}')
    volt, minUnc = routine(csv, ind)
    vActual.append(volt)
    minimaUnc.append(minUnc)

vActual = vActual[1:]
minimaUnc = minimaUnc[1:]
minima = np.arange(1, 8)

def frank_fit(minima, minimaUnc, vActual, fit_func, name):
    popt, perr, xthe, yfit = ql.fit_routine(fit_func, minima, vActual, yerr = minimaUnc)

    plt.figure(figsize = (12, 9))
    plt.errorbar(minima, vActual, yerr = minimaUnc, fmt = 'b.', markersize = 10,
                capsize = 3)
    plt.plot(xthe, yfit, 'r--')
    plt.xlabel('Minima Number', fontsize = 20)
    plt.ylabel('Measured Voltage [V]', fontsize = 20)
    plt.show()
    plt.savefig(f'excitation{name}.png')
    # plt.close()

    print(f'\nExcitation Energy: {popt[0]:.2f} +/- {perr[0]:.2f}')
    print(f'Contact Potential Difference: {popt[1]:.2f} +/- {perr[1]:.2f}')

def new_fit(x, A, B):
    return(A*(1 + B*(2*x-1)))

frank_fit(minima, minimaUnc, vActual, ql.linear_fit, 'old')
frank_fit(minima, minimaUnc, vActual, new_fit, 'new')

vapor = lambda T : 8.7 * 10**(9 - (3110 / T))
meanFP = lambda x : 1 / (np.sqrt(2)*np.pi * x * (.18e-9)**2)

for temp in [154+273, 185+273, 20+273]:
    print(f'\nVapor Pressure at {temp} K : {vapor(temp):.2f} Pa')
    numDens = np.array(vapor(temp))/ (1.38e-23 *temp)
    numDensSci = np.format_float_scientific(numDens)
    print(f'\nNumber Density at {temp} K : {numDensSci}')
    mfp = meanFP(numDens)
    print(f'\nMean Free Path at {temp} K : {mfp}')
