import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy



path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Single_trajectory/"

Ns = np.array([1e3, 1e4, 1e5, 1e6])
colors = plt.cm.Paired(np.linspace(0,1,len(Ns)))
fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.15, 'right':.95, 'bottom':.15})
for n, N in enumerate(Ns):
	my_response = Immune_response(L=15, N=int(N), alpha = 1, beta=.5, antigen_str = 'FMLFMAVFVMTSWYC', text_files_path=path, energy_model = 'MJ', bcells_filter = True)
	my_response.run(T = 25)

	ax.scatter(my_response.energies, my_response.hamming_distances, marker = '.', color = colors[n], alpha = .4)

	energy_array = np.linspace(np.min(my_response.energies), np.max(my_response.energies), 10)

	mean_l = np.array([])
	for i in range(1, len(energy_array)):
		a=np.where((my_response.energies >= energy_array[i-1]) & (my_response.energies < energy_array[i]))
		mean_l = np.append(mean_l, np.mean(my_response.hamming_distances[a]))

	ax.plot(energy_array[:-1], mean_l, marker = '*', linestyle = '-', linewidth = 2, color = colors[n], label = '%0.0e'%N)

	my_plot_layout(ax=ax, xlabel = r'$\epsilon$', ylabel = r'$l$')

	#fig.savefig('../../Figures/6_Broadness_response/ex_N-%0.0e.pdf'%N)
ax.legend(title=r'$N_r$', title_fontsize = 24, fontsize = 22, loc = 4)
fig.savefig('../../Figures/6_Broadness_response/ex_N.pdf')
	
#print(np.sum(my_response.E_matrix[my_response.antigen, my_response.master_sequence]))
