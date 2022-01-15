import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import pickle

N_A = 6.02214076e23
colors = ['tab:blue','tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = ['exponential', 'linear']
text_pos = np.array([[4e0, 5e-7],[1e2, 6e-2]], dtype = object)
text = [r'$\sim N_{b}^{-1-\frac{\lambda\alpha}{\beta}}$', r'$\sim N_{b}^{-1}$']

path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Python/ensemble/"

alpha = 1
betas = np.array([[0.1, 0.2, 0.5, 1, 2], [0.1, 0.2, 0.5]], dtype=object)
lamda = 0.6
colors_gm = np.array([plt.cm.Blues(np.linspace(0,1,len(betas[0])+2)),plt.cm.Reds(np.linspace(0,1,len(betas[1])+2)) ], dtype=object)

for j, gm in enumerate(growth_models):
	fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
	for i, beta in enumerate(betas[j]):
		exponents = [(-((alpha*lamda)/beta)-1), -1]
		file_clone_sizes = open(path+'clone_size_data_alpha-%.2f_beta-%.2f_'%(alpha, beta)+gm+'.pkl','rb')

		clone_sizes = pickle.load(file_clone_sizes)
		file_clone_sizes.close()

		clone_size_distribution = np.histogram(clone_sizes, bins = np.logspace(np.log10(np.min(clone_sizes)),np.log10(np.max(clone_sizes-1)),14), density = True)
		clone_size = ((clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]+clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]))/2
		clone_size_counts = clone_size_distribution[0][np.where(clone_size_distribution[0]!=0)]

		delta_clone_size = clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]-clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
		cumsum_clone_size_counts = np.cumsum(clone_size_counts*delta_clone_size)
		
		plaw_fit_csd = clone_size**(exponents[j])#*(np.log(clone_size**(1/nu)))**(1-lambd)
		plaw_fit_csd /= (plaw_fit_csd[-7]/(clone_size_counts[-7]))

		ax.plot(clone_size[:-1], clone_size_counts[:-1], linestyle = '', marker = 's', ms = 5, linewidth = 3, color = colors_gm[j][i+2], label = '%.2f'%(beta))

		if (gm=='exponential'):
			plaw_fit_csd = clone_size**(exponents[j])
			plaw_fit_csd /= (plaw_fit_csd[-7]/(clone_size_counts[-7]))
			ax.plot(clone_size[:][:], plaw_fit_csd[:][:], linestyle = '--', marker = '', ms = 5, linewidth = 3, alpha = .6, color = colors_gm[j][i+2])
		if (gm=='linear'):
			plaw_fit_csd = clone_size**(exponents[j])#*(((2000*alpha/N_A)**lamda)/(beta))*np.log(clone_size**(1/beta))**(lamda-1)
			plaw_fit_csd /= (plaw_fit_csd[4]/(clone_size_counts[4]))
			ax.plot(clone_size[:][:], plaw_fit_csd[:][:], linestyle = '--', marker = '', ms = 5, linewidth = 3, alpha = .6, color = colors_gm[j][i+2])

	if (gm=='exponential'):
			ax.text(x=1e0, y=2e-19, s = r'$\sim N_{b}^{-1-\frac{\lambda\alpha}{\beta}}$', fontsize=44, color = colors_gm[j][-1])
	if (gm=='linear'):
			ax.text(x=1e0, y=3e-14, s = r'$\sim N_{b}^{-1}$', fontsize=44, color = colors_gm[j][-1])
	
	my_plot_layout(ax = ax, xscale='log', yscale= 'log', xlabel=r'Clone size $N_{b}$', ylabel=r'$p(N_{b})$', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	ax.legend(title=r'$\beta/\alpha$', fontsize = 24, title_fontsize = 24)
	fig.savefig('../../Figures/1_Dynamics/CSD_betas_'+gm+'.pdf')

	#ax.text(x=text_pos[i][0], y=text_pos[i][1], s = text[i], fontsize=44, color = colors_fit[i])




