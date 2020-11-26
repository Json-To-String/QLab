import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


h = 6.626e-34
mu_B = 9.274e-24
charge = 1.602e-19
mass_e = 9.11e-31

def h_coil(V, turns, R):
    B = 8.991e-3 * V * turns / R
    return B

def b_err(B, N, R, I_unc):
    fracI = 8.991e-3 * N / R
    del_B = np.sqrt((fracI**2)*(I_unc)**2)
    return del_B

def linear_fit(x,A,B):
    return(A*x+B)

def plots(freq, B, B_err, fit_function, name):
    # freq = freq * 1000
    b_zero = h_coil(.1778, 11, 16.3e-2)
    B = (B - b_zero)

    popt, pcov = curve_fit(fit_function, freq, B, sigma = B_err)
    perr = np.sqrt(np.diag(pcov))
    div_out = charge / (4*np.pi * mass_e)
    print(div_out)
    slope = popt[0] /div_out
    slopeErr = perr[0]
    xtheory = np.linspace(min(freq), max(freq))
    yfit = fit_function(xtheory, *popt)

    plt.errorbar(freq, B, yerr = B_err, fmt = 'b.', capsize = 3)
    plt.plot(xtheory, yfit, 'r--')
    plt.title(f'{name} Field vs Rf Frequency')
    plt.xlabel('Rf Freq [kHz]')
    plt.ylabel(f'{name} Field [G]')

    print(slope, '+/-', slopeErr)
    return slope, slopeErr

## Raw
df1 = pd.read_csv('OpticalPumping.csv', index_col = None)
freq = np.array(df1[df1.columns[0]].to_list())
v_87 = np.array(df1[df1.columns[1]].to_list())
v_85 = np.array(df1[df1.columns[2]].to_list())
v_unc = np.array(df1[df1.columns[3]].to_list())
f_unc = np.array(df1[df1.columns[4]].to_list())

## Calculated
b_87 = h_coil(v_87, 3, 3.22e-2)
b_87_err = b_err(b_87, 3, 3.22e-2, v_unc)
b_85 = h_coil(v_85, 3, 3.22e-2)
b_85_err = b_err(b_85, 3, 3.22e-2, v_unc)

plt.figure(figsize = (12,9))
slope_87, slope_87_err = plots(freq, b_87, b_87_err, linear_fit, '$Rb_{87}$')
slope_85, slope_85_err = plots(freq, b_85, b_85_err, linear_fit, '$Rb_{85}$')
plt.show()
