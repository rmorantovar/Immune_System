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
antigen = 'TACNSEYPNTTKCGRWYC'
#antigen = 'FTSENAYCGR'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'

L=len(antigen)

N_ens = 500
N_r = 5e4
N_r = 1e5
#N_r = 1e6
T0 = 0
Tf = 6
Tf = 8
dT = .1
days = np.arange(0, Tf, 1)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_pr = .1 # hour^-1
k_pr = k_pr*24 #days^-1
qs = [1, 2, 3]
colors_q = ['darkred', 'olive', 'navy']
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e3
E_ms = -28

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

beta_r = lambdas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
Kd_r = np.exp(E_r)

E_pr = Es[:-1][Kds<(k_pr/k_on)][-1]
Kd_pr = np.exp(E_pr)
#----------------------------------------------------------------

lambda_Bs = np.array([np.flip([.5])*lambda_A, np.flip([.5])*lambda_A], dtype=object)

d=20
energy_model = 'MJ'
colors_gm = np.array([plt.cm.Oranges(np.linspace(0,1,len(lambda_Bs[0])+2)),plt.cm.Reds(np.linspace(0,1,len(lambda_Bs[1])+2)) ], dtype=object)
FIG, AX = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
for q in qs:
	beta_q = lambdas[lambdas>q][-1]
	E_q = Es[lambdas>q][-1]
	Kd_q = np.exp(E_q)
	Kd_act = np.max([Kd_q, Kd_r])
	for j, gm in enumerate(growth_models):
		fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		for n_lambda_B, lambda_B in enumerate(lambda_Bs[j]):

			parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_q-%d_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, q, j, N_ens)+energy_model
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

			N_active_clones_plasma = np.ones_like(N_clones)
			counter = 0
			for n, n_clones in enumerate(N_clones):
				temp_data = data[counter:counter+n_clones]
				temp_data2 = temp_data.loc[temp_data[1]==1]
				N_active_clones_plasma[n] = len(temp_data2.loc[temp_data2[2]==1])
				counter += n_clones

			Kds0 = np.exp(data[0])
			data_Kds0 = ax0.hist(Kds0, bins = np.logspace(np.log10(np.exp(min_e_data)), np.log10(np.exp(max_e_data)), 40), density = False, color = colors_gm[j][n_lambda_B+2],histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
			counts0 = data_Kds0[0][np.where(data_Kds0[0]!=0)]
			Kds_array_data0 = (data_Kds0[1][np.where(data_Kds0[0]!=0)])

			Kds = np.exp(data_plasma[0])
			data_Kds = ax0.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data)), np.log10(np.exp(max_e_data)), 40), density = False, color = colors_gm[j][n_lambda_B+2],histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
			counts = data_Kds[0][np.where(data_Kds[0]!=0)]
			Kds_array_data = (data_Kds[1][np.where(data_Kds[0]!=0)])

			clone_sizes_binned = np.zeros([len(counts)])
			clone_sizes_binned_2= np.zeros([len(counts)])
			var_clone_sizes_binned = np.zeros([len(counts)])
			max_clone_sizes_binned = np.zeros([len(counts)])

			# counter = 0
			# for n in np.arange(N_ens):
			# 	for i in np.arange(int(len(counts))):
			# 		start = counter
			# 		stop = counter + N_active_clones_plasma[n]
			# 		clone_sizes_temp = clone_sizes[start:stop]
			# 		Kds_temp = Kds[start:stop]
			# 		clone_sizes_binned[i] += np.mean( np.concatenate((clone_sizes_temp[(Kds_temp>=data_Kds[1][i]) & (Kds_temp<data_Kds[1][i+1])], np.array([0]) )) )#/N_ens
			# 		clone_sizes_binned_2[i] += np.mean( np.concatenate(((clone_sizes_temp[(Kds_temp>=data_Kds[1][i]) & (Kds_temp<data_Kds[1][i+1])])**2, np.array([0]) )) )#/N_ens
			# 		max_clone_sizes_binned[i] += np.max(clone_sizes_temp[(Kds_temp>=data_Kds[1][i]) & (Kds_temp<data_Kds[1][i+1]) ], initial=1)#/N_ens
			# 		#if(np.mean( np.concatenate((clone_sizes_temp[(Kds_temp>=data_Kds[1][i]) & (Kds_temp<data_Kds[1][i+1])], np.array([0]))))!=0):

			# 	counter = counter + N_active_clones_plasma[n]

			for i in np.arange(int(len(counts))):
		 		clone_sizes_binned[i] += np.mean( np.concatenate((clone_sizes[(Kds>=data_Kds[1][i]) & (Kds<data_Kds[1][i+1])], np.array([0]) )) )#/N_ens
		 		clone_sizes_binned_2[i] += np.mean( np.concatenate(((clone_sizes[(Kds>=data_Kds[1][i]) & (Kds<data_Kds[1][i+1])])**2, np.array([0]) )) )#/N_ens
		 		max_clone_sizes_binned[i] += np.max(clone_sizes[(Kds>=data_Kds[1][i]) & (Kds<data_Kds[1][i+1]) ], initial=1)#/N_ens

			max_clone_size = np.max(clone_sizes_binned)
			ax.plot(Kds_array_data0, counts0/N_ens, color = colors_gm[j][n_lambda_B+2], alpha = .8, marker = 'o', ms = 8, linestyle = '')
			ax.plot(Kds_array_data, counts/N_ens, color = colors_gm[j][n_lambda_B+2], alpha = .8, marker = '^', ms = 8, linestyle = '')

			popt, pcov = curve_fit(f = my_linear_func , xdata = np.log(Kds_array_data0[0:4]), ydata= np.log(counts0)[0:4] )
			beta_act2 = popt[1]
			beta_act = np.min([q, beta_r])

			print('beta_act = %.2f'%(beta_act))
			print(np.max([beta_r, beta_q]))

			u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf)/N_A, Es, q, lambda_A, N_c, dE)

			ax.plot(np.exp(Es[:-1]), Q0*N_r, color = colors_gm[j][n_lambda_B+2])
			ax.plot(np.exp(Es[:-1]), QR*N_r, color = colors_gm[j][n_lambda_B+2])
			ax.plot(Kds_array_data0[0:10], (counts0[0]/N_ens)*(Kds_array_data0[0:10]/Kds_array_data0[0])**(beta_act), color = colors_gm[j][n_lambda_B+2], linewidth = 5, linestyle = ':', marker = '', ms = 15, alpha = .8)
			ax.plot(Kds_array_data0[0:10], (counts0[0]/N_ens)*(Kds_array_data0[0:10]/Kds_array_data0[0])**(beta_act2), color = colors_gm[j][n_lambda_B+2], linewidth = 5, linestyle = ':', marker = '', ms = 15, alpha = .4)

			Kds_array_data = Kds_array_data[clone_sizes_binned!=0]
			clone_sizes_binned = clone_sizes_binned[clone_sizes_binned!=0]
			ax2.plot(Kds_array_data[:], clone_sizes_binned[:]/max_clone_size, color = 'orange', linewidth =5, linestyle = '', marker = 's', ms = 8)

			AX.plot(Kds_array_data[:], clone_sizes_binned[:]/max_clone_size, linewidth =5, linestyle = '', marker = 's', ms = 8, label = '%d'%q, color = colors_q[q-1])

			cross_over = np.where(clone_sizes_binned==max_clone_size)[0][0]
			ax2.plot(Kds_array_data[0:cross_over+1], (clone_sizes_binned[cross_over]/max_clone_size)*(Kds_array_data[0:cross_over+1]/Kds_array_data[cross_over])**((lambda_B/lambda_A)), color = 'orange', linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)
			ax2.plot(Kds_array_data[cross_over:], (clone_sizes_binned[cross_over+1]/max_clone_size)*(Kds_array_data[cross_over:]/Kds_array_data[cross_over+1])**((lambda_B/lambda_A)*(-q)), color = 'orange', linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)
			#ax2.vlines([Kd_pr, Kd_q, Kd_r], 4e-3, 1.5, linestyles = ['-',':', '--'], color = 'grey')

			#AX.plot(Kds_array_data[0:cross_over+1], (clone_sizes_binned[cross_over]/max_clone_size)*(Kds_array_data[0:cross_over+1]/Kds_array_data[cross_over])**((lambda_B/lambda_A)), linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)
			AX.plot(Kds_array_data[cross_over:], (clone_sizes_binned[cross_over]/max_clone_size)*(Kds_array_data[cross_over:]/Kds_array_data[cross_over])**((lambda_B/lambda_A)*(-q)), linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8, color = colors_q[q-1])

	my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	#ax.legend(title=r'$\lambda_A/\lambda_B$', fontsize = 30, title_fontsize = 35)
	ax.set_xlim(left = np.exp(E_ms), right = np.exp(E_ms+25))
	ax.set_ylim(bottom = 1e-6)
	fig.savefig('../../Figures/1_Dynamics/Ensemble/Q_K_q-%d.pdf'%q)

	my_plot_layout(ax = ax2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	#ax2.legend(fontsize = 30, title_fontsize = 35)
	ax2.set_xlim(left = np.exp(E_ms), right = np.exp(E_ms+25))
	ax2.set_ylim(bottom = 4e-3, top = 2e0)
	ax2.set_yticks([1, 0.1, 0.01])
	ax2.set_yticklabels([1, 0.1, 0.01])
	fig2.savefig('../../Figures/1_Dynamics/Ensemble/N_vs_K_q-%d.pdf'%q)


my_plot_layout(ax = AX, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
AX.legend(fontsize = 30, title_fontsize = 35, title = r'$q$')
AX.set_xlim(left = np.exp(E_ms), right = np.exp(E_ms+25))
AX.set_ylim(bottom = 4e-3, top = 2e0)
AX.set_yticks([1, 0.1, 0.01])
AX.set_yticklabels([1, 0.1, 0.01])
FIG.savefig('../../Figures/1_Dynamics/Ensemble/N_vs_K.pdf')





