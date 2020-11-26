import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import glob
import pandas as pd
import qLab as ql
from scipy.signal import find_peaks

matplotlib.rc('xtick', labelsize = 20)
matplotlib.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 15)

def moseley(x, A, B):
    return(A*(x-B))

def spectrum(calSample, title, bgFile):
    channel, counts = ql.read_mca(calSample)
    bgChannel, bgCounts = ql.read_mca(bgFile)

    energy = toEnergyA(channel[5:])
    corrCounts = []
    corrEnergy = []
    if not 'lead' in calSample:
        correction = counts[5:] - bgCounts[5:]
    else:
        correction = counts
    for val1, val2 in zip(correction, energy):
        if val1 > 30:
            corrCounts.append(val1)
            corrEnergy.append(val2)
        else:
            corrCounts.append(0)
            corrEnergy.append(val2)

    peaks, properties = find_peaks(corrCounts, height = 50)
    heights = properties["peak_heights"]

    # print(peakLocs)
    # print(len(peakLocs))
    # print(len(corrEnergy))
    # for ind in peakLocs:
    #     print('ind:', ind)
    #     print(corrEnergy[int(ind)])

    fig = plt.figure(figsize = (12, 9))
    ax = fig.gca()
    ax.tick_params(axis = 'both', direction = 'in',
    top = True, right = True, pad = 25)

    ax.plot(energy, corrCounts, 'b')

    ax.set_title(f'{title}: Counts vs Energy', fontsize = 20)
    ax.set_xlabel('Energy [keV]', fontsize = 20)
    ax.set_ylabel('Counts', fontsize = 20)
    # plt.show()
    plt.savefig(f'{title}.png')
    plt.close()


# df0 = pd.read_csv('Xray.csv', index_col = None)
df0 = ql.getCsv(r'https://docs.google.com/spreadsheets/d/1LOnJFSP1Hn6V9N8JQve2ok9sezGbgc4ToZTjgjdCXxc/edit#gid=0')
energy1, energy2 = ql.colVar(df0, 1), ql.colVar(df0, 2)
peak1, peak1Unc = ql.colVar(df0, 3), ql.colVar(df0, 4)
peak2, peak2Unc = ql.colVar(df0, 5), ql.colVar(df0, 6)

calSamples = ql.stringVar(df0, 7)
cPeak1, cPeak1Unc = ql.colVar(df0, 8), ql.colVar(df0, 9)
cPeak2, cPeak2Unc = ql.colVar(df0, 10), ql.colVar(df0, 11)

poptA, perrA, xtheA, yfitA = ql.fit_routine(ql.linear_fit, peak1, energy1,
                                            start = 0, stop = 320, step = 321)
poptB, perrB, xtheB, yfitB = ql.fit_routine(ql.linear_fit, peak2, energy2,
                                            start = 0, stop = 400, step = 400)

toEnergyA = lambda channels : ql.linear_fit(channels, *poptA)
toEnergyB = lambda channels : ql.linear_fit(channels, *poptB)

fig = plt.figure(figsize = (12, 9))
ax = fig.gca()
ax.tick_params(axis = 'both', direction = 'in',
               top = True, right = True,
               pad = 25)
ax.plot(peak1, energy1, 'k.', markersize = 10, label = 'Calibration')
for item1, item2, item3 in zip(cPeak1, toEnergyA(cPeak1), calSamples):
    ax.plot(np.linspace(0, item1), np.ones(50)*item2,
            label = f'{item3} - {item2:.2f} keV')
ax.plot(xtheA, yfitA, 'r--')
ax.set_title('Channels to Energy Calibration for $K_{\\alpha}$ ', fontsize = 20)
ax.set_xlabel('Channel', fontsize = 20)
ax.set_ylabel('Energy [keV]', fontsize = 20)
ax.legend()
# plt.show()
plt.savefig('kAlphaCalib.png')
plt.close()


fig = plt.figure(figsize = (12, 9))
ax = fig.gca()
ax.tick_params(axis = 'both', direction = 'in',
               top = True, right = True,
               pad = 25)
ax.plot(peak2, energy2, 'k.', markersize = 10, label = 'Calibration')
for item1, item2, item3 in zip(cPeak2, toEnergyA(cPeak2), calSamples):
    ax.plot(np.linspace(0, item1), np.ones(50)*item2,
            label = f'{item3} - {item2:.2f} keV')
ax.plot(xtheB, yfitB, 'r--')
ax.set_title('Channels to Energy Calibration for $K_{\\beta}$', fontsize = 20)
ax.set_xlabel('Channel', fontsize = 20)
ax.set_ylabel('Energy [keV]', fontsize = 20)
ax.legend()
# plt.show()
plt.savefig('kBetaCalib.png')
plt.close()

samples = os.path.join('Data', 'Spectra', '*.txt')
bgPath = os.path.join('Data','background.txt')
for item in glob.glob(samples):
    name = os.path.basename(item).replace('.txt', '').capitalize()
    print(f'{name}: ')
    spectrum(item, name, bgPath)

atomicNum = [26, 29, 28, 30, 40, 42]
moseEn = np.sqrt(toEnergyA(peak1))
moseUnc = np.sqrt(toEnergyA(peak1Unc))
poptM, perrM, xtheM, yfitM = ql.weighted_fit(moseley, atomicNum, moseEn,
                                             moseUnc)
fig = plt.figure(figsize = (12, 9))
ax = fig.gca()
ax.tick_params(axis = 'both', direction = 'in',
               top = True, right = True,
               pad = 25)
ax.errorbar(atomicNum, moseEn , yerr = moseUnc, capsize = 3, fmt = 'k.')
ax.plot(xtheM, yfitM, 'r--')
ax.set_xlabel('Atomic Number', fontsize = 20)
ax.set_ylabel('$E^{1/2}$ [$keV^{1/2}$]', fontsize = 20)
# plt.show()
plt.savefig('Moseley.png')

print(f'C: {poptM[0]:.5f} +/- {perrM[0]:.5f}')
print(f'S: {poptM[1]:.2f} +/- {perrM[1]:.2f}')
