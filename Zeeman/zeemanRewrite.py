import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import qLab as ql

wave = 644e-9
speedLight = 3e8
n = 1.46
h = 6.626e-34

def deltaE(alpha1, alpha2):
    beta = lambda alpha: np.arcsin(np.sin(alpha)/n)

    A = - h * speedLight/wave
    B = (np.cos(beta(alpha1)) / np.cos(beta(alpha2))) - 1
    return(A*B)

df0 = ql.getCsv('https://docs.google.com/spreadsheets/d/1jBiHycLMfRqhLOjewJboAS7nLtQQeuBogt-y-UEfge0/edit#gid=0')
col = lambda num : ql.colVar(df0, num)

calCurr, calCurrUnc = col(0), col(2)
calField, calFieldUnc = col(1), col(3)
popt, pcov, xthe, yfit = ql.fit_routine(ql.third_order, calCurr, calField, yerr = calField)

# plt.figure()
# plt.plot(calCurr, calField, 'k.' )
# plt.plot(xthe, yfit, 'r--')
# plt.show()

currents = col(4)
fields = ql.third_order(currents, *popt)

up_alpha_left = np.radians(col(5))
up_alpha_right = np.radians(col(6))
up_angle = np.radians(col(8))

up_leftE = deltaE(up_angle, up_alpha_left)
up_rightE = deltaE(up_angle, up_alpha_right)
