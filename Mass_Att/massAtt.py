import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import interpolate
from scipy.optimize import curve_fit

## Notes for Jason:
# For Co-60
# peak1 is at 1173 keV
# peak2 is at 1333 keV

df0 = pd.read_csv('massNotes.csv',index_col = None)

df1 = pd.read_csv('Al_Log.csv',index_col = None)
al_thick = df1[df1.columns[0]].to_list()
al_p1na = df1[df1.columns[1]].to_list()
al_p2na = df1[df1.columns[2]].to_list()

df2 = pd.read_csv('Pb_Log.csv',index_col = None)
pb_thick = df2[df2.columns[0]].to_list()
pb_p1na = df2[df2.columns[1]].to_list()
pb_p2na = df2[df2.columns[2]].to_list()

df3 = pd.read_csv('nist_al_coeff.txt',index_col = None,delimiter = ' ')
al_energy = df3[df3.columns[0]].to_list()
al_muRho = df3[df3.columns[2]].to_list()

df4 = pd.read_csv('nist_pb_coeff.txt',index_col = None,delimiter = ' ')
pb_energy = df4[df4.columns[0]].to_list()
pb_muRho = df4[df4.columns[2]].to_list()

def linear_fit(x,A,B):
    return(A*x+B)

def theory_val(en,mu):
    peak1 = 1.173 # MeV
    peak2 = 1.333 # MeV
    # print('energy',en)
    # print('muval',mu)

    f = interpolate.interp1d(en,mu)
    xnew = np.linspace(min(en),max(en),10000)
    ynew = f(xnew)
    muTheory1 = f(peak1)
    muTheory2 = f(peak2)

    # plt.figure()
    # plt.loglog(en,mu,'b.')
    # plt.loglog(xnew,ynew,'r-')
    # plt.axvline(peak1)
    # plt.axvline(peak2)
    # plt.show()

    return(muTheory1,muTheory2)

def log_plot(thick,N,fit_function,name):

    time = 4*60 # seconds
    log_intensity = np.log(np.array(N)/time)

    nErr = log_intensity*np.sqrt(np.array(N))/N
    popt,pcov = curve_fit(fit_function,thick,log_intensity,sigma = nErr)
    perr = np.sqrt(np.diag(pcov))
    slopeErr = perr[0]
    # p_weight = np.poly1d(popt)
    xtheory = np.linspace(min(thick),max(thick))
    yfit = fit_function(xtheory,*popt)

    plt.figure(figsize = (12,9))
    plt.errorbar(thick,log_intensity,yerr = nErr,fmt = 'b.',capsize = 3)
    plt.plot(xtheory,yfit,'r--')
    plt.xlabel('Thickness [cm]',fontsize = 14)
    plt.ylabel('ln(I) [ln(counts/s)]',fontsize = 14)
    # plt.show()
    plt.savefig(f'{name}.png')
    plt.close()
    thickfun = np.log(2)/(-1*popt[0])
    dthickdmu = -np.log(2)/np.sqrt(-1*popt[0])
    dthick = np.sqrt((dthickdmu**2)*(slopeErr**2))
    print(thickfun,'+/-',dthick)
    return(popt[0],slopeErr)

almuExp1,alUnc1 = log_plot(al_thick,al_p1na,linear_fit,'AlPeak1')
almuExp2,alUnc2 = log_plot(al_thick,al_p2na,linear_fit,'AlPeak2')
almuTh1,almuTh2 = theory_val(al_energy, al_muRho)
pbmuExp1,pbUnc1 = log_plot(pb_thick,pb_p1na,linear_fit,'PbPeak1')
pbmuExp2,pbUnc2 = log_plot(pb_thick,pb_p2na,linear_fit,'PbPeak2')
pbmuTh1,pbmuTh2 = theory_val(pb_energy, pb_muRho)


al_rawData = pd.DataFrame({
    'Al Thickness': al_thick,
    'Al Peak1 Net Area': al_p1na,
    'Al Peak2 Net Area': al_p2na,
    })
pb_rawData = pd.DataFrame({
    'Pb Thickness': pb_thick,
    'Pb Peak1 Net Area': pb_p1na,
    'Pb Peak2 Net Area': pb_p2na,
    })
measData = pd.DataFrame({
    'Measured mu peak 1': [almuExp1*(-1),pbmuExp1*(-1)],
    'Peak1 Uncs': [alUnc1,pbUnc1],
    'Measured mu peak 2': [almuExp2*(-1),pbmuExp2*(-1)],
    'Peak2 Uncs': [alUnc2,pbUnc2],
    'Accepted mu peak 1': [almuTh1*2.7,pbmuTh1*11.3],
    'Accepted mu peak 2': [almuTh2*2.7,pbmuTh2*11.3]
})

# print(al_rawData.to_latex())
# print(pb_rawData.to_latex())
# print(measData.to_latex())

# print('\n\t\t\tAluminum Peak 1:')
# print('Meas:\t',almuExp1*(-1),'\tUnc:',alUnc1,'\tAcc:\t',almuTh1*2.7)
# print('\n\t\t\tAluminum Peak 2:')
# print('Meas:\t',almuExp2*(-1),'\tUnc:',alUnc2,'\tAcc:\t',almuTh2*2.7)
# print('\n\t\t\tLead Peak 1:')
# print('Meas:\t',pbmuExp1*(-1),'\tUnc:',pbUnc1,'\tAcc:\t',pbmuTh1*11.3)
# print('\n\t\t\tLead Peak 2:')
# print('Meas:\t',pbmuExp2*(-1),'\tUnc:',pbUnc2,'\tAcc:\t',pbmuTh2*11.3)
