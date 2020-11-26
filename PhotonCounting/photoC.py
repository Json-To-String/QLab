import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import factorial
import qLab as ql

def poisson(k, lamb):
    return (lamb**k / factorial(k)) * np.exp(-lamb)

def boseEinstein(n, nAv):
    return ( (nAv**n) / (nAv + 1)**(n + 1) )

def chiSq(theory, exper, uncert, bins):
    chi = []
    for the, exp, unc, in zip(theory, exper, uncert):
        if unc != 0:
            chi.append( ((exp - the)**2) / (bins * (unc**2)) )
    summed = np.sum(chi)
    return(summed)

def routine(csv, fitFunc):
    name = csv[5:]
    print(name)
    df0 = pd.read_csv(csv, index_col = None)
    counts = ql.colVar(df0, 0)
    countsUnc = np.sqrt(counts)
    
    avg = np.mean(counts)

    binSize = np.arange(0, max(counts), 1)

    newCounts, newBins = np.histogram(counts, bins = binSize - 0.5)
    binErr = np.sqrt(newCounts) / len(counts)
    newCounts = newCounts / len(counts)
    newBins = newBins[1:] - 0.5
    
    yTheory = fitFunc(binSize, avg)

    plt.figure(figsize = (12,9))
    plt.bar(newBins, newCounts, width = 1)
    plt.errorbar(newBins, newCounts, yerr = binErr, fmt = 'k.', 
                capsize = 3) 
    plt.plot(binSize, yTheory, 'r.')
    plt.xlabel('Photon Counts', fontsize = 13)
    plt.ylabel('Probability of Count', fontsize = 13)
    plt.title(f'Probability Distribution: {name}')
    #plt.show()
    plt.close()

    #print(len(binSize), len(yTheory), len(newBins), len(newCounts))
   
    #for val1, val2 in zip(newCounts, yTheory[1:]):
    #    print(val1, val2)
    
    chiS = chiSq(newCounts, yTheory, binErr, len(newBins))

    print('Chi Squared:', chiS)

for csv in glob.glob(os.path.join('Runs','*')):
    if 'thermal' in csv:
        routine(csv, boseEinstein)
    else:
        routine(csv, poisson)


