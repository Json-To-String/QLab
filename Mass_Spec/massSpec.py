import qLab as ql
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.signal import find_peaks

def routine(filename):
    tups = os.path.splitext(filename)
    nameRun = tups[0][5:].split('_')
    name = nameRun[0]
    run = nameRun[1]
    bgfile = f'BG/{name}_bg.txt'
    
    df0 = pd.read_csv(filename, delim_whitespace = True)
    df1 = pd.read_csv(bgfile, delim_whitespace = True)
    
    x = ql.colVar(df0, 0)
    x = x[17:] 
    xbg = ql.colVar(df1, 0)
    xbg = xbg[17:]

    y = ql.colVar(df0, 1)
    y = y[17:] 
    ybg = ql.colVar(df1, 1)
    ybg = ybg[17:]

    mass = np.array([float(xval[:-1]) for xval in x])
    
    press = np.array([float(yval[:-1]) - float(bgval[:-1]) for yval,
        bgval in zip(y, ybg)])*10**6

    peaks, _ = find_peaks(press, height = .25, distance = 5)

    plt.figure(figsize = (12, 9))
    plt.plot(mass, press, 'b')
    plt.plot(mass[peaks], press[peaks], 'ro', markersize = 5)
    plt.xlabel('Mass [amu]', fontsize = 14)
    plt.ylabel('Pressure [$\mu$Torr]', fontsize = 14)
    plt.title(f'{name.capitalize()} run {run}', fontsize = 14)
    plt.savefig(f'{name}run{run}.png')
    #plt.show()
    plt.close()

for filename in glob.glob('Data/*.txt'):
    routine(filename)


