import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qLab as ql

df1 = pd.read_csv('PhotonCount.csv', index_col = None)

discLevel = np.abs(ql.colVar(df1, 0))
discUnc = ql.colVar(df1, 1)

countSignal = ql.colVar(df1, 2)
signalUnc = np.sqrt(countSignal)

countBg = ql.colVar(df1, 3) / 10
bgUnc = np.sqrt(countBg)

signalNoise = (countSignal - countBg) / countBg
signalNoiseUnc = signalNoise * np.sqrt((signalUnc/countSignal)**2 + (bgUnc/countBg)**2)

plt.figure(figsize = (12,9))
plt.errorbar(discLevel, signalNoise, yerr = signalNoiseUnc, fmt = 'b.', capsize = 3)
plt.plot(discLevel, signalNoise, 'b.')
plt.show()

