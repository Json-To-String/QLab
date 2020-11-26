import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit
import sys
import glob

## Convert HH:MM:SS to seconds
def get_sec(time_str):
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + int(s)

def linear_fit(x,A,B):
    return(A*x + B)

def log_plot(x,y,err,fit_function,sample):
    popt,pcov = curve_fit(fit_function,x,y,sigma = err)
    perr = np.sqrt(np.diag(pcov))
    xtheory = np.linspace(min(x),max(x))
    yfit = fit_function(xtheory,*popt)

    plt.figure(figsize = (12,9))
    plt.errorbar(x,y,yerr = err,fmt = 'b.',capsize = 3)
    plt.plot(xtheory,yfit,'r--')
    plt.xlabel('Time [s]',fontsize = 14)
    plt.ylabel('ln(C) [ln(counts)]',fontsize = 14)
    plt.savefig(f'{sample}_logcountsvstime.png')
    # plt.show()
    plt.close()

    slope = popt[0]
    slopeUnc = perr[0]

    return(slope,slopeUnc)

def nonLinear_plot(x,y,err,sample):
    plt.figure(figsize = (12,9))
    plt.errorbar(x,y,yerr = err,fmt = 'b.',capsize = 3)
    plt.xlabel('Time [s]',fontsize = 14)
    plt.ylabel('Counts',fontsize = 14)
    plt.savefig(f'{sample}_countsvstime.png')
    # plt.show()
    plt.close()

# Count Correction for Barium
def nTrue(nObs,rObs,deadTime):
    nObs = np.array(nObs)
    C = 1/(1-rObs*deadTime)
    return C*nObs

def compute_halfLife(factor,factorUnc):
    halfLife = np.log(2)/factor
    fracUnc = -1*np.log(2)/(factor**2)
    halfLifeUnc = np.sqrt( (fracUnc**2)*(factorUnc**2) )
    return halfLife,halfLifeUnc

def read_mca(filename,iters):
    with open(filename,'r') as mca:
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

    rawCounts_unc = np.sqrt(mcaCounts)
    timeSpan = 2
    rateObs = mcaChannel/timeSpan
    deadTime = 400e-6
    corr_counts = nTrue(mcaCounts,rateObs,deadTime)
    fracUnc = 1/(1-rateObs*deadTime)
    corr_counts_unc = np.sqrt((fracUnc**2)*(rawCounts_unc**2))
    time = mcaChannel*timeSpan

    lnCtrue = np.log(corr_counts)
    lnCtrue_unc = np.sqrt( ((1/corr_counts)**2) * (corr_counts_unc**2) )
    nonLinear_plot(time,corr_counts,corr_counts_unc,f'barium{iters}')
    slope,slopeUnc = log_plot(time,lnCtrue,lnCtrue_unc,linear_fit,f'barium{iters}')
    # print('Barium Slope:',slope,'+/-',slopeUnc)
    halfLifeBa,halfLifeBaUnc = compute_halfLife(slope*-1,slopeUnc)
    print(halfLifeBa,'+/-',halfLifeBaUnc)

## Indium
# Background counts data
bg_ind = np.array([95,100])
bg_counts = np.mean(bg_ind)
bg_counts_unc = np.sqrt(bg_ind[0]) + np.sqrt(bg_ind[1])

# Measured counts data
df0 = pd.read_csv('halfLife.csv',index_col = None)
ind_times_inconvenient = df0[df0.columns[0]].to_list()

## Convert times
ind_times = []
for item in ind_times_inconvenient:
    ind_times.append(get_sec(item))

ind_times = np.array(ind_times) - ind_times[0]

# Subtract out background counts
ind_meas_counts = df0[df0.columns[1]].to_list()
ind_meas_counts_unc = []
for item in ind_meas_counts:
    ind_meas_counts_unc.append(np.sqrt(float(item)))

ind_net_counts = np.array(ind_meas_counts) - np.array(bg_counts)
ind_net_counts_unc = bg_counts_unc + ind_meas_counts_unc

# Get ln of counts for second plot
ind_lnC = np.log(ind_net_counts)
ind_lnC_unc = np.sqrt( ((1/np.array(ind_net_counts))**2)
                        * (np.array(ind_net_counts_unc)**2) )

nonLinear_plot(ind_times,ind_net_counts,ind_net_counts_unc,'indium')
slopeInd,slopeIndUnc = log_plot(ind_times,ind_lnC,ind_lnC_unc,linear_fit,'indium')
# print('Indium Slope:',slopeInd,'+/-',slopeIndUnc)
halfLifeInd,halfLifeIndUnc = compute_halfLife(slopeInd*-1,slopeIndUnc)
print(halfLifeInd,'+/-',halfLifeIndUnc)

## Barium
textFiles = glob.glob(os.path.join('*.txt'))
for ind,file in enumerate(textFiles):
    print(file,ind+1)
    read_mca(file,ind+1)
