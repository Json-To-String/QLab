import numpy as np
import matplotlib.pyplot as plt
import qLab as ql
import pandas as pd

## Calculation for population ratios
popRatio = np.exp(-(5.5857) * (3.152e-8) * (3000e-4) / (1/40))
print(popRatio)

def routine(*args):
    # Converts Gauss to Tesla
    toTesla = lambda G : G * 1e-4

    df0 = pd.read_csv('Resonance.csv', index_col = None)

    calCurr = ql.colVar(df0, 0)
    calField = toTesla(np.abs(ql.colVar(df0, 1)))

    calOpt, calErr, calFit, calTh = ql.weighted_fit(ql.linear_fit, calCurr,
                                                    calField,
                                                    np.array(len(calCurr)*[5e-4]))

    # Lambdas to find frequencies to test around for a given current value
    optCurr = lambda B : (toTesla(B) - calOpt[1]) / calOpt[0]
    optBField = lambda I : (calOpt[0] * I) + calOpt[1]
    optFreqB = lambda B, g : g * ql.muNEvT * toTesla(B) / ql.heVs
    optFreqI = lambda I, g : g * ql.muNEvT * optBField(I) / ql.heVs

    # Graphs where slopes are related to the g-factor by a factor of 1/muN
    curr, freq = ql.colVar(df0, args[0]), ql.colVar(df0, args[1])*1e6
    ener = freq * ql.heVs
    currUnc, freqUnc = ql.colVar(df0, args[2]), ql.colVar(df0, args[3])
    field = optBField(curr)

    popt, perr, xtheory, yfit = ql.weighted_fit(ql.no_intercept, field, ener,
                                                freqUnc)

    plt.figure(figsize = (12, 9))
    plt.errorbar(field, ener, yerr = freqUnc*ql.heVs, fmt = 'k.', capsize = 5)
    plt.plot(xtheory, yfit, 'r--')
    plt.xlabel('Field [T]', fontsize = 20)
    plt.ylabel('Resonance Frequency [MHz]', fontsize = 20)
    plt.savefig(f'{args[4]}.png')

    gFac = popt[0] / ql.muNEvT
    gUnc = perr[0] / ql.muNEvT

    print('\nG-Factor: ', gFac, '+/-', gUnc)

routine(2, 3, 4, 5, 'water')
routine(6, 7, 8, 9, 'fluorine')
