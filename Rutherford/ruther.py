import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import qLab as ql
from scipy.optimize import curve_fit
from itertools import zip_longest


def sineFit(theta, A, theta0, B):
    return( A / (np.sin( (theta + theta0)*np.pi/360 ))**4 + B ) 

def grouper(n, iterable, padvalue = None):
    return(zip_longest(*[iter(iterable)]*n, fillvalue = padvalue))

def routine(row1, row2, mate):
    df0 = pd.read_csv('Rutherford.csv', index_col = None)
    
    theta = ql.colVar(df0, row1)
    counts = ql.colVar(df0, row2)

    # Remove NaNs because Google Sheets is booty
    angular = []
    number = []
    for angle, count in zip(theta, counts):
        if not np.isnan(angle):
            angular.append(angle)
            number.append(count / 5)
       
    angGroup = grouper(5, angular, padvalue = 'Bad')
    numGroup = grouper(5, number, padvalue = 'Bad')

    # Take avg of each of the counts 
    angAvg = []
    numAvg = []
    for elem1, elem2 in zip(angGroup, numGroup):
        angAvg.append(np.mean(elem1))
        numAvg.append(np.mean(elem2))
    

    popt1, pcov1 = curve_fit(sineFit, angular, number)
    print(popt1)
    perr1 = np.sqrt(np.diag(pcov1))
   
    angle = 5
    xtheory1 = np.concatenate((np.linspace(min(angular), -angle), np.linspace(angle, max(angular))), axis=None)
    yfit1 = sineFit(xtheory1, *popt1)
    
    # Raw plot
    plt.figure(figsize = (12, 9))
    plt.errorbar(angular, number, xerr = len(angular)*[1] , yerr = np.sqrt(number), fmt = 'ko', capsize = 3)
    plt.plot(xtheory1, yfit1, 'r*')
    plt.title(f'Counts/s vs Angle for {mate} Foil', fontsize = 20,
            fontweight = 'bold')
    plt.xlabel('Incident Angle [Degrees]', fontsize = 14)
    plt.ylabel('Counts/s', fontsize = 14)
    plt.savefig(f'{mate}.png')
    #plt.show()
    plt.close()
    
    popt2, pcov2 = curve_fit(sineFit, angAvg, numAvg)
    perr2 = np.sqrt(np.diag(pcov2))
    
    xtheory2 = np.linspace(min(angAvg), angAvg[len(angAvg)//2])
    yfit2 = sineFit(xtheory2, *popt2)

    # Avged plot
    plt.figure(figsize = (12, 9))
    plt.plot(angAvg, numAvg, 'k.', markersize = 14)
    plt.plot(xtheory2, yfit2, 'r--')
    plt.title(f'Averaged Counts/s vs Angle for {mate} Foil', fontsize = 20,
            fontweight = 'bold')
    plt.xlabel('Incident Angle [Degrees]', fontsize = 14)
    plt.ylabel('Avg Counts/s', fontsize = 14)
    #plt.show()
    plt.savefig(f'avg{mate}.png')
    plt.close()

routine(1, 2, 'Gold')
#routine(4, 5, 'Aluminum')
