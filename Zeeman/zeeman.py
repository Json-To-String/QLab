import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

h = 6.626e-34
c = 3e8
n = 1.46
lambdaCd = 644e-9

def cal_curve(x, A, B, C, D):
    return(A*x**3 + B*x**2 + C*x + D)

def beta(alpha):
    b = np.arcsin(np.sin(alpha)/n)
    return b

def delta_Lambda(alpha1, alpha2):
    cb1 = np.cos(beta(alpha1))
    cb2 = np.cos(beta(alpha2))
    return((cb2/cb1) - 1)

def delta_E(al1, al2):
    dLaByLa2 = delta_Lambda(al1, al2)/lambdaCd
    return(-h*c*dLaByLa2)

## Calibration Curves
df1 = pd.read_csv('Zeeman.csv', index_col = None)
calCurr = np.array(df1[df1.columns[0]].to_list())
calField = np.array(df1[df1.columns[1]].to_list())

popt, pcov = curve_fit(cal_curve, calCurr, calField)
perr = np.sqrt(np.diag(pcov))
slopeErr = perr[0]
yfit = cal_curve(calCurr, *popt)

#plt.figure(figsize = (12,9))
#plt.plot(calCurr, calField, 'b.')
#plt.plot(calCurr, yfit, 'r--')
#plt.show()

##########

# Convert current into B-field
currents = np.array(df1[df1.columns[4]].to_list())
currents = currents[:4]
fields = cal_curve(currents, *popt)

up_alpha_left = np.radians(np.array(df1[df1.columns[5]].to_list()[:4]))
up_alpha_right =  np.radians(np.array(df1[df1.columns[6]].to_list()[:4]))
up_angle = np.radians(np.array(df1[df1.columns[8]].to_list()[:4]))

up_leftE = delta_E(up_angle, up_alpha_left)
up_rightE = delta_E(up_angle, up_alpha_right)

totalFields = np.append(fields, fields)
totalUp = np.append(np.abs(up_leftE), up_rightE)

plt.figure()
plt.plot(totalFields, totalUp, "b.")
plt.show()

down_alpha_left =  np.array(df1[df1.columns[9]].to_list()[:4])
down_alpha_right =  np.array(df1[df1.columns[10]].to_list()[:4])
down_angle = np.array(df1[df1.columns[12]].to_list()[:4])

down_leftE = delta_E(down_angle, down_alpha_left)
down_rightE = delta_E(down_angle, down_alpha_right)

totalFields = np.append(fields, fields)
totalDown = np.append(down_leftE, np.abs(down_rightE))

plt.figure()
plt.plot(totalFields, totalDown, "b.")
plt.show()

