import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import pickle

colors = ['tab:blue','tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = ['exponential', 'linear']
text_pos = np.array([[4e0, 5e-7],[1e2, 6e-2]], dtype = object)
text = [r'$\sim N_{b}^{-1-\frac{\lambda\alpha}{\beta}}$', r'$\sim N_{b}^{-1}$']

path = "../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Python/ensemble/"

alphas = [0.2, 0.5, 1.0]
beta = 0.5
lamda = 0.62
fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'bottom':.15})

for i, alpha in enumerate(alphas):
	exponents = [(-((alpha*lamda)/beta)-1), -1]
	file_clone_sizes = open(path+'clone_size_data_alpha-%.2f_beta-%.2f_'%(alpha, beta)+growth_models[0]+'.pkl','rb')

	clone_sizes = pickle.load(file_clone_sizes)
	file_clone_sizes.close()

	clone_size_distribution = np.histogram(clone_sizes, bins = np.logspace(np.log10(np.min(clone_sizes)),np.log10(np.max(clone_sizes-1)),14), density = True)
	clone_size = ((clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]+clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]))/2
	clone_size_counts = clone_size_distribution[0][np.where(clone_size_distribution[0]!=0)]

	delta_clone_size = clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]-clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
	cumsum_clone_size_counts = np.cumsum(clone_size_counts*delta_clone_size)
	
	plaw_fit_csd = clone_size**(exponents[0])#*(np.log(clone_size**(1/nu)))**(1-lambd)
	plaw_fit_csd /= (plaw_fit_csd[-3]/(clone_size_counts[-3]))



	ax.plot(clone_size[:-1][:], plaw_fit_csd[:-1][:], linestyle = '--', marker = '', ms = 8, linewidth = 3, alpha = .9, color = colors[0])
	ax.plot(clone_size[:-1], clone_size_counts[:-1], linestyle = '-', marker = '', ms = 4, linewidth = 2, color = colors_fit[0])


	#ax.text(x=text_pos[i][0], y=text_pos[i][1], s = text[i], fontsize=44, color = colors_fit[i])


my_plot_layout(ax = ax, xscale='log', yscale= 'log', xlabel=r'Clone size $N_{b}$', ylabel=r'$p(N_{b})$', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
fig.savefig('../../Figures/1_Dynamics/CSD_alphas.pdf')

