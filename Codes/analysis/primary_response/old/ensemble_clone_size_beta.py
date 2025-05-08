import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['text.usetex'] = True
import scipy as scipy
import pickle

N_A = 6.02214076e23
colors = ['tab:blue','tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = ['exponential']#, 'linear']

text_pos = np.array([[4e0, 5e-7],[1e2, 6e-2]], dtype = object)
text = [r'$\sim N_{b}^{-1-\frac{\lambda\lambda_A}{\lambda_B}}$', r'$\sim N_{b}^{-1}$']

Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/"

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'TACNSEYPNTTK'
#antigen = 'FTSENAYCGR'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'

L=len(antigen)

N_ens = 100
N_r = 1e4
T0 = 0
Tf = 8
dT = .1
days = np.arange(0, Tf, 1)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_act = .1 # hour^-1
k_act = k_act*24 #days^-1
q = 2
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e3
E_ms = -27
time = np.linspace(T0, Tf, int((Tf-T0)/dT))

print('k_on/k_act = %.1e'%(k_on/k_act))

#lambda_A = 1
lambda_Bs = np.array([np.flip([1.2])*lambda_A, np.flip([1.2])*lambda_A], dtype=object)

lamda = 0.82
lambd = 1.5

d=20
energy_model = 'MJ'
colors_gm = np.array([plt.cm.Oranges(np.linspace(0,1,len(lambda_Bs[0])+2)),plt.cm.Reds(np.linspace(0,1,len(lambda_Bs[1])+2)) ], dtype=object)

for j, gm in enumerate(growth_models):
	fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
	fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
	for i, lambda_B in enumerate(lambda_Bs[j]):
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_q-%d_linear-%d_'%(lambda_A, 0.5, k_act/24, q, j)+energy_model
		data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		data_active = data.loc[data[1]==1]
		activations_times = np.array(data_active[3])

		clone_sizes = np.exp(lambda_B*(Tf - activations_times))
		min_e_data = np.min(data[0])
		max_e_data = np.max(data[0])

		Kds0 = np.exp(data[0])
		data_Kds0 = ax0.hist(Kds0, bins = np.logspace(np.log10(np.exp(min_e_data-1)), np.log10(np.exp(max_e_data+1)), 20), density = False, histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
		counts0 = data_Kds0[0][np.where(data_Kds0[0]!=0)]
		Kds_array_data0 = (data_Kds0[1][np.where(data_Kds0[0]!=0)])
		popt, pcov = curve_fit(f = my_linear_func , xdata = np.log(Kds_array_data0[:15]), ydata= np.log(counts0)[:15] )
		print('lambda = %.2f'%(popt[1]))
		lambd = popt[1]

		exponents = [(-((lambda_A*0.8)/lambda_B)-1), -1]
		exponents2 = [(-((lambda_A*0.8)/lambda_B)-1), -1]

		clone_size_distribution = np.histogram(clone_sizes, bins = np.logspace(np.log10(np.min(clone_sizes)),np.log10(np.max(clone_sizes)),15), density = True)
		clone_size = ((clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]+clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]))/2
		clone_size_counts = clone_size_distribution[0][np.where(clone_size_distribution[0]!=0)]

		delta_clone_size = clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]-clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
		cumsum_clone_size_counts = np.cumsum(clone_size_counts*delta_clone_size)
		
		plaw_fit_csd = clone_size**(exponents[j])#*(np.log(clone_size**(1/nu)))**(1-lambd)
		plaw_fit_csd /= (plaw_fit_csd[-3]/(clone_size_counts[-3]))

		plaw_fit_csd2 = clone_size**(exponents2[j])#*(np.log(clone_size**(1/nu)))**(1-lambd)
		plaw_fit_csd2 /= (plaw_fit_csd2[-1]/(clone_size_counts[-1]))

		ax.plot(clone_size[:], clone_size_counts[:], linestyle = '', marker = 's', ms = 8, linewidth = 5, color = 'orange', label = '%.1f'%(lambda_A/lambda_B))

		if (gm=='exponential'):
			
			ax.plot(clone_size[:], plaw_fit_csd[:], linestyle = '--', marker = '', ms = 5, linewidth = 3, alpha = .8, color = 'orange')
			#ax.plot(clone_size[-12:], plaw_fit_csd2[-12:], linestyle = ':', marker = '', ms = 5, linewidth = 3, alpha = .6, color = colors_gm[j][i+2])
		if (gm=='linear'):
			ax.plot(clone_size[:][:], plaw_fit_csd[:][:], linestyle = '--', marker = '', ms = 5, linewidth = 3, alpha = .6, color = 'orange')
	#ax.plot(clone_size[:], (clone_size[:]**(-1)/clone_size[-1]**(-1))*clone_size_counts[-1], linestyle = '--', marker = '', ms = 5, linewidth = 4, alpha = .6, color = 'darkred')
	
	
	my_plot_layout(ax = ax, xscale='log', yscale= 'log', y_fontsize=30 )
	ax.set_xlim(clone_size[0]*.5, clone_size[-1]*2)
	ax.set_ylim(clone_size_counts[-1]*.1, clone_size_counts[0]*10)
	ax.legend(title=r'$\lambda_A/\lambda_B$', fontsize = 30, title_fontsize = 35)
	fig.savefig('../../Figures/1_Dynamics/Ensemble/CSD_betas_'+gm+'.pdf')

	#ax.text(x=text_pos[i][0], y=text_pos[i][1], s = text[i], fontsize=44, color = colors_fit[i])




