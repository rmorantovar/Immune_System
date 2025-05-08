import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import pandas as pd
import numpy as np
from pymuller import muller

# Read in the data.
df = pd.read_csv('../../out/affinity_maturation/trajectory_0.csv', index_col=0)
df = df.sample(frac=1)
metadf = pd.read_csv('../../out/affinity_maturation/metadata_0.csv')

times = [float(i) for i in df.keys()]
# for t in df.keys():
#     df[t]/=np.sum(df[t])

# print(np.sum(df[df.keys()[-1]]))

fig, ax = plt.subplots()
ax.stackplot(times, df)
fig.savefig('../../../Figures/affinity_maturation/muller.png')

fig, ax = plt.subplots()
ax.plot(times, np.mean((metadf['fitness']*df.T).T, axis = 0))
fig.savefig('../../../Figures/affinity_maturation/trajectories.png')

fig, ax = plt.subplots()
freqs = df[df.keys()[-1]].tolist()
freqs = [i for i in freqs if i != 0]

bins = np.logspace(np.log10(np.min(freqs))/1.1, np.log10(np.max(freqs))*1.1, 30)
data_spectra = np.histogram(freqs, bins = bins, density = False)
ax.plot(data_spectra[1][:-1], 1-np.cumsum(data_spectra[0])/(np.sum(data_spectra[0])))
ax.plot(data_spectra[1], data_spectra[1]**-0.75/data_spectra[1][0]**-0.75)
ax.set_xscale('log')
ax.set_yscale('log')
fig.savefig('../../../Figures/affinity_maturation/sprectrum.png')


# populations_df = pd.read_csv('../../out/affinity_maturation/populations_df.csv', index_col=0)
# adjacency_df = pd.read_csv('../../out/affinity_maturation/adjacency_df.csv', index_col=0)
# identities = populations_df['Identity'].unique()
# color_by_dict = {i:i for i in list(identities)}
# color_by_df = pd.Series(data=color_by_dict, index=[i for i in list(identities)])
# fig, ax = plt.subplots()

# muller(populations_df, adjacency_df, ax = ax, color_by = color_by_df, smoothing_std = .01, normalize = True, colormap = 'turbo')
# ax.set_yscale('linear')
# fig.savefig('../../../Figures/affinity_maturation/muller_2.png')
