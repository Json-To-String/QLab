import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def colVar(df, colNum):
    return(np.array(df[df.columns[colNum]].to_list()))

def linear_fit(x, A, B):
    return(A*x + B)

def quadratic_fit(x, A, B, C):
    return(A*x**2 + B*x + C)

def weighted_fit(x, y, err, fit_func):
    xtheory = np.linspace(min(x), max(x))
    popt, pcov = curve_fit(fit_func, x, y, sigma = err)
    perr = np.sqrt(np.diag(pcov))
    yfit = fit_function(xtheory, *popt)

    return(popt, perr, yfit, xtheory)
    

