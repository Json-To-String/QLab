import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

heVs = 4.136e-15 # Plancks Constant in eV * s
muNEvT = 3.152e-8 # Nuclear Magneton in eV / T
gProton = 5.5857 # G Factor for the proton
gFluorine = 5.25 # G Factor for Fluorine

def colVar(df, colNum):
    return(np.array([x for x in df[df.columns[colNum]].to_list()
            if not np.math.isnan(x)]))

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
def fit_routine(x, y, fit_func):

    # popt are the parameters found by the least squares
    # pcov encodes the errors/uncertainties in the popt params
    popt, pcov = curve_fit(fit_func, x, y)
    perr = np.sqrt(np.diag(pcov))

    xtheory = np.linspace(np.sort(x)[0], np.sort(x)[-1])
    yfit = fit_func(xtheory, *popt)

    return(xtheory, yfit)
