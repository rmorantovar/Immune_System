import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as scipy
import pickle
import time
from tqdm import tqdm

def neutralization(a, b, x):
	return 1/(1+((np.exp(a*(x-1.2)))/(b)))

N_epitopes = np.arange(1, 12, 1)

colors = plt.cm.plasma(np.linspace(0, 1, len(N_epitopes)+3))

N_ensemble = 4000
T = 10

fig_fast1, ax_fast1 = plt.subplots(figsize=(12, 8))
fig_fast2, ax_fast2 = plt.subplots(figsize=(12, 8), gridspec_kw={'bottom':.12, 'left':.10, 'right':.90})
fig_slow1, ax_slow1 = plt.subplots(figsize=(12, 8))

gain_function_fast = np.array([])
gain_function_slow = np.array([])

for l in np.arange(len(N_epitopes)): #Running loop over number of epitopes
	time = np.linspace(0, T, T*N_epitopes[l]+1)
	cross_protection = np.zeros_like(time)

	for n in np.arange(N_ensemble): #Running ensemble
		mutations = np.random.randint(low=1, high=N_epitopes[l]+1, size = T*N_epitopes[l])
		cross_protection_temp = 1-np.array([np.product([1-neutralization(2.5, 10, np.count_nonzero(mutations[0:t]==i)) for i in N_epitopes[0: l+1]]) for t in np.arange(len(time))])
		cross_protection = cross_protection + cross_protection_temp

	gain_function_fast = np.append(gain_function_fast, (cross_protection[time<=4][-1])/N_ensemble)
	gain_function_slow = np.append(gain_function_slow, (cross_protection[time<=1][-1])/N_ensemble)


	ax_fast1.plot(time, cross_protection/N_ensemble, label = '%d'%(N_epitopes[l]), color = colors[l+1], linewidth = 3)
	ax_fast1.vlines(4, 0, 1, linestyle = 'dashed', color='black')
	ax_slow1.plot(time*4, cross_protection/N_ensemble, label = '%d'%(N_epitopes[l]), color = colors[l+1], linewidth = 3)
	ax_slow1.vlines(4, 0, 1, linestyle = 'dashed', color='black')
	

my_plot_layout(ax_fast1, xlabel = 'Time [$1/\mu$]', ylabel = 'Cross-protection')
ax_fast1.set_xlim(right=10.5)
fig_fast1.savefig('../../Figures/4_Complexity/Cross-protection_nEpitopes_fast.pdf')	

my_plot_layout(ax_slow1, xlabel = 'Time [$1/\mu$]', ylabel = 'Cross-protection')
ax_slow1.set_xlim(right=10.5)
fig_slow1.savefig('../../Figures/4_Complexity/Cross-protection_nEpitopes_slow.pdf')

cost_function = my_linear_func(N_epitopes, 0, .05)
fitness_fast = gain_function_fast - cost_function
fitness_slow = gain_function_slow - cost_function

df_fast = pd.DataFrame(list(zip(N_epitopes, gain_function_fast)), columns = ['l', 'F_r'])
df_fast.to_csv('/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Complexity/F_fast.txt', sep = '\t', index = False)

df_slow = pd.DataFrame(list(zip(N_epitopes, gain_function_slow)), columns = ['l', 'F_r'])
df_slow.to_csv('/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Complexity/F_slow.txt', sep = '\t', index = False)


