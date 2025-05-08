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

Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/"


antigen = 'FMLFMAVFVMTSWYC'
L=len(antigen)
N=2e3
alphas = np.array([[0.5, 1.0, 2.0, 5.0], [0.2, 0.5, 1.0, 2.0]], dtype=object)
beta = 1
gamma = 0
lamda = 0.6
d=20
energy_model = 'MJ'
colors_gm = np.array([plt.cm.Blues(np.linspace(0,1,len(alphas[0])+2)),plt.cm.Reds(np.linspace(0,1,len(alphas[1])+2)) ], dtype=object)

for j, gm in enumerate(growth_models):
	fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
	for i, alpha in enumerate(alphas[j]):
		exponents = [(-((alpha*lamda)/beta)-1), -1]
		
		#file_clone_sizes = open(path+'/Dynamics/Python/ensemble/clone_size_data_d-%d_alpha-%.2f_beta-%.2f_N-%.1e_'%(d, alpha, beta, N)+gm+'_'+energy_model+'.pkl','rb')
		#clone_sizes = pickle.load(file_clone_sizes)
		clone_sizes = np.loadtxt(Text_files_path + 'Dynamics/Ensemble/bcells_ensemble_L-%d_N-%d_Antigen-'%(L, N)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_linear-%d_'%(alpha, beta, gamma, j)+energy_model+'.txt')
		clone_sizes = clone_sizes/np.max(clone_sizes)
		#file_clone_sizes.close()

		clone_size_distribution = np.histogram(clone_sizes, bins = np.logspace(np.log10(np.min(clone_sizes)),np.log10(np.max(clone_sizes)),14), density = True)
		clone_size = ((clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]+clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]))/2
		clone_size_counts = clone_size_distribution[0][np.where(clone_size_distribution[0]!=0)]

		delta_clone_size = clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]-clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
		cumsum_clone_size_counts = np.cumsum(clone_size_counts*delta_clone_size)
		
		plaw_fit_csd = clone_size**(exponents[j])#*(np.log(clone_size**(1/nu)))**(1-lambd)
		plaw_fit_csd /= (plaw_fit_csd[-7]/(clone_size_counts[-7]))

		ax.plot(clone_size[:-1], clone_size_counts[:-1], linestyle = '', marker = 's', ms = 5, linewidth = 3, color = colors_gm[j][i+2], label = '%.1f'%(alpha/beta))

		if (gm=='exponential'):
			plaw_fit_csd = clone_size**(exponents[j])
			plaw_fit_csd /= (plaw_fit_csd[-3]/(clone_size_counts[-3]))
			ax.plot(clone_size[:][:], plaw_fit_csd[:][:], linestyle = '--', marker = '', ms = 5, linewidth = 3, alpha = .6, color = colors_gm[j][i+2])
		if (gm=='linear'):
			plaw_fit_csd = clone_size**(exponents[j])#*(((2000*alpha/N_A)*lamda)/(beta))*np.log(clone_size**(1/beta))**(lamda-1)
			plaw_fit_csd /= (plaw_fit_csd[3]/(clone_size_counts[3]))
			ax.plot(clone_size[:][:], plaw_fit_csd[:][:], linestyle = '--', marker = '', ms = 5, linewidth = 3, alpha = .6, color = colors_gm[j][i+2])

	if (gm=='exponential'):
			ax.text(x=1e0, y=2e-19, s = r'$\sim N_{b}^{-1-\frac{\lambda\alpha}{\beta}}$', fontsize=44, color = colors_gm[j][-1])
	if (gm=='linear'):
			ax.text(x=1e0, y=3e-14, s = r'$\sim N_{b}^{-1}$', fontsize=44, color = colors_gm[j][-1])
	
	my_plot_layout(ax = ax, xscale='log', yscale= 'log', xlabel=r'Clone size $N_{b}$', ylabel=r'$p(N_{b})$', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	ax.legend(title=r'$\alpha/\beta$', fontsize = 24, title_fontsize = 24)
	fig.savefig('../../Figures/1_Dynamics/CSD_alphas_'+gm+'.pdf')

	#ax.text(x=text_pos[i][0], y=text_pos[i][1], s = text[i], fontsize=44, color = colors_fit[i])




