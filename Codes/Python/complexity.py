import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import pickle
import time
from tqdm import tqdm

def neutralization(a, b, x):
	return 1/(1+((np.exp(a*(x-1.2)))/(b)))

N_epitopes = np.arange(1, 15, 1)

colors = plt.cm.plasma(np.linspace(0, 1, len(N_epitopes)+3))

N_ensemble = 1000

T = 10

fig, ax = plt.subplots(figsize=(12, 8) )
fig2, ax2 = plt.subplots(figsize=(12, 8), gridspec_kw={'bottom':.12, 'left':.10, 'right':.90})

gain_function = np.array([])

for l in np.arange(len(N_epitopes)):
	time = np.linspace(0, T, T*N_epitopes[l]+1)
	cross_protection = np.zeros_like(time)
	for n in np.arange(N_ensemble):
		mutations = np.random.randint(low=1, high=N_epitopes[l]+1, size = T*N_epitopes[l])
		cross_protection_temp = 1-np.array([np.product([1-neutralization(2.5, 10, np.count_nonzero(mutations[0:t]==i)) for i in N_epitopes[0: l+1]]) for t in np.arange(len(time))])
		cross_protection = cross_protection + cross_protection_temp

	
	
	#if(N_epitopes[l]==1):
	#	cross_protection_1 = 1-np.array([(1-neutralization(2.5, 10, np.count_nonzero(mutations[0:t]==1))) for t in np.arange(len(time))])
	#	time_1 = time

	gain_function = np.append(gain_function, (cross_protection[time<=4][-1])/N_ensemble)


	plot1 = ax.plot(time, cross_protection/N_ensemble, label = '%d'%(N_epitopes[l]), color = colors[l+1], linewidth = 3)
	ax.vlines(4, 0, 1, linestyle = 'dashed', color='black')
	#ax.plot(time_1, 1-(-cross_protection_1+1)**N_epitopes[l], linestyle = 'dashed', color = plot1[0].get_color())

	#x = np.linspace(0, T, 100)
	#ax.plot(x, 1-(1-neutralization(2.5, 10, x/np.sqrt(N_epitopes[l])))**N_epitopes[l], linestyle = 'dotted', color = plot1[0].get_color())

my_plot_layout(ax, xlabel = 'Time [$1/\mu$]', ylabel = 'Cross-protection')
#ax.legend(fontsize = 15, title = '# of epitopes $l$', title_fontsize = 24)
fig.savefig('../../Figures/4_Complexity/Cross-protection_nEpitopes.pdf')	


ax2.plot(N_epitopes, gain_function, color = 'darkgreen', lw = 4, ls = 'dotted', alpha = .8)
ax2.tick_params(axis='y', labelcolor='darkgreen')
ax2_1 = ax2.twinx()
cost_function = my_linear_func(N_epitopes, 1, .3)
ax2_1.plot(N_epitopes, -cost_function, color = 'darkred', lw = 4, ls = 'dotted', alpha = .8)
ax2_1.tick_params(axis='y', labelcolor='darkred')

my_plot_layout(ax2)
my_plot_layout(ax2_1)
fig2.savefig('../../Figures/4_Complexity/fig5-2.pdf')	