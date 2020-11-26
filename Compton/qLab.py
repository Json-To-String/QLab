import errno
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import requests
from scipy.optimize import curve_fit

## Useful Constants ##
heVs = 4.136e-15 # Plancks Constant in eV * s
muNEvT = 3.152e-8 # Nuclear Magneton in eV / T
gProton = 5.5857 # G Factor for the proton
gFluorine = 5.25 # G Factor for Fluorine
massElectron = 511e3 # eV


toRadians = lambda deg : deg * np.pi / 180

def colVar(df, colNum):
    finalCol = []
    for data in df[df.columns[colNum]].to_list():
        if isinstance(data, str):
            finalCol.append(data)
        elif np.math.isnan(data):
            pass
        else:
            finalCol.append(data)
    
    return(np.array(finalCol))

#def colVar(df, colNum):
#    dataList = np.array([x for x in df[df.columns[colNum]].to_list()if not np.math.isnan(x)])
#    return(dataList)

def stringVar(df, colNum):
    return(df[df.columns[colNum]].to_list())

def linear_fit(x, A, B):
    return(A*x + B)

def no_intercept(x, A):
    return(A*x)

def quadratic_fit(x, A, B, C):
    return(A*x**2 + B*x + C)

def weighted_fit(fit_func, x, y, err):
    xtheory = np.linspace(min(x), max(x))
    popt, pcov = curve_fit(fit_func, x, y, sigma = err)
    perr = np.sqrt(np.diag(pcov))
    yfit = fit_func(xtheory, *popt)

    return(popt, perr, xtheory, yfit)

# This routine fits a function to a set of data using least-squares
def fit_routine(fit_func, x, y, yerr = [None], start = None, stop = None,
                step = None):

    # popt are the parameters found by the least squares
    # pcov encodes the errors/uncertainties in the popt params

    # Conditional to pass in linspace params
    if start == None:
        xtheory = np.linspace(np.sort(x)[0], np.sort(x)[-1])
    else:
        xtheory = np.linspace(start, stop, step)

    if yerr[0] == None:
        popt, pcov = curve_fit(fit_func, x, y)
    else:
        popt, pcov, = curve_fit(fit_func, x, y, sigma = yerr)
    perr = np.sqrt(np.diag(pcov))


    yfit = fit_func(xtheory, *popt)

    return(popt, perr, xtheory, yfit)

def read_mca(mcaFile):
    with open(mcaFile,'r') as mca:
        counter = 0
        mcaChannel = []
        mcaCounts = []
        for line in mca:
            if 'Channel' in line:
                counter +=1
                continue
            if counter > 0:
                mcaChannel.append(float(line.split('\t')[0]))
                mcaCounts.append(float(line.split('\t')[1]))
        mcaChannel = np.array(mcaChannel)
        mcaCounts = np.array(mcaCounts)

    # print(mcaChannel, mcaCounts)
    return mcaChannel, mcaCounts

def getCsv(url):
    # Take in the url for the PUBLIC google spreadsheet and saves it in
    # the corresponding path and returns a dataframe with it inside
    # Example usage:
    # df = qLab.getCsv("url", "path")

    # Finds the spreadsheet's key
    key = re.search(r"/spreadsheets/d/([a-zA-Z0-9-_]+)", url).group(1)

    # Grabs the first sheet from the spreadsheet determined by the url
    response = requests.get(f'https://docs.google.com/spreadsheets/d/{key}/gviz/tq?tqx=out:csv&sheet=Sheet1')
    assert response.status_code == 200, 'Wrong status code'

    path = f'./Data/{key}.csv'
    if not os.path.exists(os.path.dirname(path)):
        try:
            os.makedirs(os.path.dirname(path))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

    # Writes the csv into a file as strings, may cause problems because it writes floats as strings
    with open(path, 'w+') as csv:
        csv.write(response.text)

    df = pd.read_csv(path)

    return(df)
