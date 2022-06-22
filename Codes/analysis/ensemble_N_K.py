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
text = [r'$\sim N_{b}^{-1-\frac{\beta\lambda_A}{\lambda_B}}$', r'$\sim N_{b}^{-1}$']

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

N_ens = 500
N_r = 5e4
T0 = 0
Tf = 6
dT = .1
days = np.arange(0, Tf, 1)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_pr = .1 # hour^-1
k_pr = k_pr*24 #days^-1
qs = [1, 2]
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e3
E_ms = -27

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

beta_r = lambdas[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
Kd_r = np.exp(E_r)
#----------------------------------------------------------------


#lambda_A = 1
lambda_Bs = np.array([np.flip([.5])*lambda_A, np.flip([.5])*lambda_A], dtype=object)

d=20
energy_model = 'MJ'
colors_gm = np.array([plt.cm.Oranges(np.linspace(0,1,len(lambda_Bs[0])+2)),plt.cm.Reds(np.linspace(0,1,len(lambda_Bs[1])+2)) ], dtype=object)

for q in [1, 2]:
	beta_q = lambdas[lambdas>q][-1]
	E_q = Es[lambdas>q][-1]
	K_d_q = np.exp(E_q)
	Kd_act = np.max([K_d_q, Kd_r])
	for j, gm in enumerate(growth_models):
		fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		for n_lambda_B, lambda_B in enumerate(lambda_Bs[j]):

			parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_q-%d_linear-%d_'%(lambda_A, 0.5, k_pr/24, q, j)+energy_model
			data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
			data2 = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/summary_ensemble.txt', sep = '\t', header=None)
			N_clones = np.array(data2[0])
			data_active = data.loc[data[1]==1]
			data_plasma = data_active.loc[data_active[2]==1]
			data_GC = data_active.loc[data_active[2]==0]
			activations_times = np.array(data_plasma[3])
			min_e_data = np.min(data[0])
			max_e_data = np.max(data[0])
			clone_sizes = np.exp(lambda_B*(Tf - activations_times))

			N_clones_active = np.array([], dtype = int)
			N_ens = len(N_clones)
			counter = 0
			for n in np.arange(N_ens):
				temp_data = data[counter:counter+int(N_clones[n])]
				N_clones_active = np.append(N_clones_active, int(len(np.array(temp_data.loc[temp_data[1]==1][0]))))
				counter = counter + int(N_clones[n])


			Kds0 = np.exp(data[0])
			data_Kds0 = ax0.hist(Kds0, bins = np.logspace(np.log10(np.exp(min_e_data)), np.log10(np.exp(max_e_data)), 30), density = False, color = colors_gm[j][n_lambda_B+2],histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
			counts0 = data_Kds0[0][np.where(data_Kds0[0]!=0)]
			Kds_array_data0 = (data_Kds0[1][np.where(data_Kds0[0]!=0)])

			Kds = np.exp(data_plasma[0])
			data_Kds = ax0.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data)), np.log10(np.exp(max_e_data)), 30), density = False, color = colors_gm[j][n_lambda_B+2],histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
			counts = data_Kds[0][np.where(data_Kds[0]!=0)]
			Kds_array_data = (data_Kds[1][np.where(data_Kds[0]!=0)])

			clone_sizes_binned = np.zeros([len(counts)])
			clone_sizes_binned_2= np.zeros([len(counts)])
			var_clone_sizes_binned = np.zeros([len(counts)])
			max_clone_sizes_binned = np.zeros([len(counts)])
			for n in np.arange(N_ens-1):
				counter = 0
				for i in np.arange(int(len(counts))):
					clone_sizes_binned[i] += np.mean( np.concatenate((clone_sizes[counter:counter + N_clones_active[n]][(Kds[counter:counter + N_clones_active[n]]>=data_Kds[1][i]) & (Kds[counter:counter + N_clones_active[n]]<data_Kds[1][i+1])], np.array([0]) )) )/N_ens
					clone_sizes_binned_2[i] += np.mean( np.concatenate(((clone_sizes[counter:counter + N_clones_active[n]][(Kds[counter:counter + N_clones_active[n]]>=data_Kds[1][i]) & (Kds[counter:counter + N_clones_active[n]]<data_Kds[1][i+1])])**2, np.array([0]) )) )/N_ens
					max_clone_sizes_binned[i] += np.max(clone_sizes[counter:counter + N_clones_active[n]][(Kds[counter:counter + N_clones_active[n]]>=data_Kds[1][i]) & (Kds[counter:counter + N_clones_active[n]]<data_Kds[1][i+1]) ], initial=1)/N_ens
				counter = counter + N_clones_active[n]

			ax.plot(Kds_array_data0, counts0/N_ens, color = colors_gm[j][n_lambda_B+2], alpha = .8, marker = 'o', ms = 10, linestyle = '')
			ax.plot(Kds_array_data, counts/N_ens, color = colors_gm[j][n_lambda_B+2], alpha = .8, marker = '^', ms = 10, linestyle = '')
			
			#clone_sizes_binned = clone_sizes_binned[Kds_array_data>Kd_act]
			#clone_sizes_binned_2 = clone_sizes_binned_2[Kds_array_data>Kd_act]
			#counts = counts[Kds_array_data>Kd_act]
			#Kds_array_data = Kds_array_data[Kds_array_data>Kd_act]

			#counts0 = counts0[Kds_array_data0>Kd_act]
			#Kds_array_data0 = Kds_array_data0[Kds_array_data0>Kd_act]

			popt, pcov = curve_fit(f = my_linear_func , xdata = np.log(Kds_array_data0[0:4]), ydata= np.log(counts0)[0:4] )
			print('beta_act = %.2f'%(popt[1]))
			beta_act2 = popt[1]
			beta_act = np.min([q, beta_r])
			#popt2, pcov2 = curve_fit(f = my_linear_func , xdata = np.log(Kds_array_data0[0:-2]), ydata= np.log(counts0)[0:-2] )
			#print('beta = %.2f'%(popt2[1]))
			#lambd2 = popt2[1]

			u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf)/N_A, Es, q, lambda_A, N_c, dE)

			ax.plot(np.exp(Es[:-1]), Q0*N_r, color = colors_gm[j][n_lambda_B+2])
			ax.plot(np.exp(Es[:-1]), QR*N_r, color = colors_gm[j][n_lambda_B+2])
			#ax0.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data-1)), np.log10(np.exp(max_e_data)), 10), density = False, color = colors_gm[j][n_lambda_B+2], alpha = .2, histtype = 'step', weights = clone_sizes, zorder=0, align = 'left', linewidth = 2)
			ax.plot(Kds_array_data0[0:5], (counts0[0]/N_ens)*(Kds_array_data0[0:5]/Kds_array_data0[0])**(beta_act), color = colors_gm[j][n_lambda_B+2], linewidth =5, linestyle = ':', marker = '', ms = 15, alpha = .8)
			ax.plot(Kds_array_data0[0:5], (counts0[0]/N_ens)*(Kds_array_data0[0:5]/Kds_array_data0[0])**(beta_act2), color = colors_gm[j][n_lambda_B+2], linewidth =5, linestyle = ':', marker = '', ms = 15, alpha = .4)
			#ax.vlines([Kd_r, Kds[lambdas[:-1]<q][0], Kds[np.where(lambdas[:-1]<beta_r)][0]], ax2.get_ylim()[0], ax2.get_ylim()[1], colors = ['black','grey','darkred'], linestyles = ['--', ':', ':'], linewidths = [2, 2, 2])

			#ax.plot(Kds_array_data[4:], counts0[4]/N_ens*(Kds_array_data[4:]/np.min(Kds_array_data[4:]))**(lambd2), color = colors_gm[j][n_lambda_B+2], linewidth =2, linestyle = ':', marker = '', ms = 15, alpha = .8)
			#ax.vlines(np.exp(Es)[lambdas < q][0], ax2.get_ylim()[0], 1e4, color = 'black', linestyle = ':', linewidth = 2)
			#ax.vlines(k_pr/k_on, 0, 1e4, color = 'grey', linestyle = ':')
			ax2.plot(Kds_array_data[clone_sizes_binned!=0][:-1], clone_sizes_binned[clone_sizes_binned!=0][:-1], color = 'orange', linewidth =5, linestyle = '', marker = 's', ms = 8)
			#ax2.plot(Kds_array_data[clone_sizes_binned!=1], max_clone_sizes_binned[max_clone_sizes_binned!=1], color = 'orange', linewidth =5, linestyle = '', marker = '*', ms = 8)
			#ax2.plot(Kds_array_data[:-4], clone_sizes_binned[-6]*(Kds_array_data[:-4]/Kds_array_data[-6])**((lambda_B/alpha)*(-q+lambd)), color = 'orange', linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)                
			ax2.plot(Kds_array_data[:-1], clone_sizes_binned[6]*(Kds_array_data[:-1]/Kds_array_data[6])**((lambda_B/lambda_A)*(-q)), color = 'orange', linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)
			#ax2.plot(Kds_array_data[:-1], clone_sizes_binned[4]*(Kds_array_data[:-1]/Kds_array_data[4])**((lambda_B/lambda_A)*(-beta_act2/q)), color = 'orange', linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .4)                
			#ax2.plot(Kds_array_data[:-4], clone_sizes_binned[-6]*(Kds_array_data[:-4]/Kds_array_data[-6])**((lambda_B/lambda_A)*(-q+lambd2)), color = 'orange', linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)                
			#ax2.plot(Kds_array_data[:-1], clone_sizes_binned[-4]*(Kds_array_data[:-1]/Kds_array_data[-4])**((lambda_B/lambda_A)*(-q*lambd2)), color = 'orange', linewidth =3, linestyle = ':', marker = '', ms = 15, alpha = .8)                
			#ax2.errorbar(x = Kds_array_data[:], y = clone_sizes_binned[:], yerr=1.9*np.sqrt(clone_sizes_binned_2 - clone_sizes_binned**2 ), xerr=None, linestyle = '', linewidth = 2, capsize = 8, ecolor = 'orange')
			#ax2.plot(Kds_array_data[:5], max_clone_sizes_binned[0]*(Kds_array_data[:5]/np.min(Kds_array_data))**((lambda_B/lambda_A)*(lambd)), color = colors_gm[j][n_lambda_B+2], linewidth =2, linestyle = '--', marker = '', ms = 15, alpha = .8)
			#ax2.vlines(np.exp(Es)[lambdas < q][0], ax2.get_ylim()[0], 10, color = 'black', linestyle = ':', linewidth = 2)
			#ax2.vlines(k_pr/k_on, 0, 10, color = 'grey', linestyle = ':')
			#ax2.vlines([Kd_r, Kds[lambdas[:-1]<q][0], Kds[np.where(lambdas[:-1]<beta_r)][0]], ax2.get_ylim()[0], ax2.get_ylim()[1], colors = ['black','grey','darkred'], linestyles = ['--', ':', ':'], linewidths = [2, 2, 2])

		my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		#ax.legend(title=r'$\lambda_A/\lambda_B$', fontsize = 30, title_fontsize = 35)
		ax.set_ylim(bottom = 1e-6)
		fig.savefig('../../Figures/1_Dynamics/Ensemble/Q_K_q-%d.pdf'%q)

		my_plot_layout(ax = ax2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		#ax2.legend(title=r'$\lambda_A/\lambda_B$', fontsize = 30, title_fontsize = 35)
		ax2.set_ylim(bottom = 0.5)
		fig2.savefig('../../Figures/1_Dynamics/Ensemble/N_vs_K_q-%d.pdf'%q)





