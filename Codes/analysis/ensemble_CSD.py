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

M1 = np.loadtxt(Text_files_path+'MJ.txt', skiprows= 1, usecols=range(1,21)).tolist()
M2 = (np.loadtxt(Text_files_path+'MJ2.txt', skiprows= 1, usecols=range(1,21))).tolist()
M3 = np.loadtxt(Text_files_path+'BLOSUM62.txt', skiprows= 1, max_rows = 23, usecols=range(1,24)).tolist()
Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
Alphabet2 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w']
Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)

Matrix = 'MJ2'
#Matrix = 'MM'

if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'TACNSEYPNTTK'
#antigen = 'FTSENAYCGR'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'

L=len(antigen)

N_ens = 100
N_r = 5e4
T0 = 0
Tf = 6
dT = .1
days = np.arange(0, Tf, 1)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_pr = .1 # hour^-1
k_pr = k_pr*24 #days^-1
qs = [1, 2, 3]
lambda_B = .5*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e3
E_ms = -27
time = np.linspace(T0, Tf, int((Tf-T0)/dT))

print('k_on/k_pr = %.1e'%(k_on/k_pr))

#----------------------------------------------------------------

antigen_list = [i for i in antigen]
antigen_seq = np.array([], dtype = int)
for i, aa in enumerate(antigen_list):
    index = Alphabet_list.index(aa)
    antigen_seq = np.append(antigen_seq, int(index))
PWM_data = M2[:,antigen_seq]

#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

Es, dE, Q0, lambdas = calculate_Q0(0.01, 50, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

beta_r = lambdas[:-1][np.cumsum(Q0*dE)<(1/(N_r*N_ens))][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r*N_ens))][-1]
Kd_r = np.exp(E_r)
#----------------------------------------------------------------

lambda_Bs = np.array([np.flip([.5])*lambda_A, np.flip([.5])*lambda_A], dtype=object)

d=20
energy_model = 'MJ'
colors_gm = np.array([plt.cm.Oranges(np.linspace(0,1,len(lambda_Bs[0])+2)),plt.cm.Reds(np.linspace(0,1,len(lambda_Bs[1])+2)) ], dtype=object)
for q in qs:

	beta_q = lambdas[lambdas>q][-1]
	E_q = Es[lambdas>q][-1]
	K_d_q = np.exp(E_q)
	Kd_act = np.max([K_d_q, Kd_r])
	#----------------------------------------------------------------
	for j, gm in enumerate(growth_models):
		fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		for i, lambda_B in enumerate(lambda_Bs[j]):
			parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_q-%d_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, q, j, N_ens)+energy_model
			data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
			data_active = data.loc[data[1]==1]
			activations_times = np.array(data_active[3])

			clone_sizes = np.exp(lambda_B*(Tf - activations_times))
			min_e_data = np.min(data[0])
			max_e_data = np.max(data[0])

			Kds0 = np.exp(data[0])
			data_Kds0 = ax0.hist(Kds0, bins = np.logspace(np.log10(np.exp(min_e_data)), np.log10(np.exp(max_e_data)), 20), density = False, histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
			counts0 = data_Kds0[0][np.where(data_Kds0[0]!=0)]
			Kds_array_data0 = (data_Kds0[1][np.where(data_Kds0[0]!=0)])
			popt, pcov = curve_fit(f = my_linear_func , xdata = np.log(Kds_array_data0[0:4]), ydata= np.log(counts0)[0:4] )
			beta_act2 = popt[1]
			beta_act = np.min([q, beta_r])
			print('beta_act = %.2f'%(beta_act))

			exponents = [(((lambda_A*beta_act)/(lambda_B*q))+1), -1]
			exponents2 = [(((lambda_A*beta_act2)/(lambda_B))+1), -1]

			clone_size_distribution = np.histogram(clone_sizes, bins = np.logspace(np.log10(np.min(clone_sizes)),np.log10(np.max(clone_sizes)),15), density = True)
			clone_size = ((clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]+clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]))/2
			clone_size_counts = clone_size_distribution[0][np.where(clone_size_distribution[0]!=0)]

			delta_clone_size = clone_size_distribution[1][1:][np.where(clone_size_distribution[0]!=0)]-clone_size_distribution[1][:-1][np.where(clone_size_distribution[0]!=0)]
			cumsum_clone_size_counts = np.cumsum(clone_size_counts*delta_clone_size)
			
			plaw_fit_csd = clone_size**(-exponents[j])#*(np.log(clone_size**(1/nu)))**(1-beta_act)
			plaw_fit_csd /= (plaw_fit_csd[-3]/(clone_size_counts[-3]))

			plaw_fit_csd2 = clone_size**(-exponents2[j])#*(np.log(clone_size**(1/nu)))**(1-beta_act)
			plaw_fit_csd2 /= (plaw_fit_csd2[-1]/(clone_size_counts[-1]))

			ax.plot(clone_size[:], clone_size_counts[:], linestyle = '', marker = '^', ms = 10, linewidth = 5, color = 'orange', label = '%.1f'%(lambda_A/(lambda_B*q)))

			if (gm=='exponential'):
				ax.plot(clone_size[:], plaw_fit_csd[:], linestyle = '--', marker = '', ms = 5, linewidth = 3, alpha = .8, color = 'orange')
			if (gm=='linear'):
				ax.plot(clone_size[:][:], plaw_fit_csd[:][:], linestyle = '--', marker = '', ms = 5, linewidth = 3, alpha = .6, color = 'orange')
		#ax.plot(clone_size[:], (clone_size[:]**(-1)/clone_size[-1]**(-1))*clone_size_counts[-1], linestyle = '--', marker = '', ms = 5, linewidth = 4, alpha = .6, color = 'darkred')
		
		
		my_plot_layout(ax = ax, xscale='log', yscale= 'log', y_fontsize=30 )
		ax.set_xlim(clone_size[0]*.5, clone_size[-1]*2)
		ax.set_ylim(clone_size_counts[-1]*.1, clone_size_counts[0]*10)
		ax.legend(title=r'$\frac{\lambda_A}{q\lambda_B}$', fontsize = 30, title_fontsize = 35)
		fig.savefig('../../Figures/1_Dynamics/Ensemble/CSD_q-%d.pdf'%q)





