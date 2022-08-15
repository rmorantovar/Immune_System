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

N_ens = 10
N_r = 5e4
N_r = 1e6
#N_r = 1e6
T0 = 0
Tf = 10
#Tf = 8
dT = .01
days = np.arange(0, Tf, 1)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_pr = .1 # hour^-1
k_pr = k_pr*24 #days^-1
thetas = [2, 1.5, 1.0]
thetas = [2, 1.5]
colors_theta = ['darkred', 'olive', 'navy']
colors_theta = ['darkred', 'olive', 'darkblue']
colors_theta = ['tab:blue','darkblue', 'olive', 'orange', 'darkred']
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
E_ms = -28
C = 1e4

print('k_pr/k_on = %.1e'%(k_on/k_pr)**(-1))

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

Es, dE, Q0, betas = calculate_Q0(0.01, 50, 200000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
Kd_r = np.exp(E_r)

print('beta_r = %.2f'%beta_r)

E_pr = Es[:-1][Kds<(k_pr/k_on)][-1]
Kd_pr = np.exp(E_pr)
#----------------------------------------------------------------

lambda_Bs = np.array([np.flip([.5])*lambda_A, np.flip([.5])*lambda_A], dtype=object)
delta_Nb = lambda t, tb, Nb, N, lambda_B, C: lambda_B*Nb*(1-(N/C))*np.heaviside(t-tb, 1)

d=20
energy_model = 'MJ'
colors_gm = np.array([plt.cm.Oranges(np.linspace(0,1,len(lambda_Bs[0])+2)),plt.cm.Reds(np.linspace(0,1,len(lambda_Bs[1])+2)) ], dtype=object)
FIG, AX = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
FIG2, AX2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
FIG3, AX3 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
for i_theta, theta in enumerate(thetas):
	print('theta = %.1f...'%theta)
	beta_q = betas[betas>theta][-1]
	E_q = Es[betas>theta][-1]
	Kd_q = np.exp(E_q)
	Kd_act = np.max([Kd_q, Kd_r])
	for j, gm in enumerate(growth_models):
		fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		for n_lambda_B, lambda_B in enumerate(lambda_Bs[j]):

			#--------------------------m(t)---------------------------
			u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, theta, lambda_A, N_c, dE)
			M_r = N_r*N_c*np.sum(Q0*p_a*dE)
			m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])

			t_act = time[m_bar>1][0]
			t_C = t_act+1.2
			#--------------------------------------------------------

			parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, theta, j, N_ens)+energy_model
			data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
			data2 = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/summary_ensemble.txt', sep = '\t', header=None)
			N_clones = np.array(data2[0])

			min_e_data = np.min(data[0])
			max_e_data = np.max(data[0])

			Kds0 = np.exp(data[0])	
			NC = np.zeros_like(time)
			for i_ens in np.arange(N_ens):

				data_i = data.loc[data[4]==i_ens]
				data_active = data_i.loc[data_i[1]==1]
				#data_active = data_active.loc[data_active[3]<t_C]
				energies = data_active[0]
				data_plasma = data_active.loc[data_active[2]==1]
				data_GC = data_active.loc[data_active[2]==0]
				activation_times = np.array(data_active[3]) #Can be changed for data_plasma
				
				clone_sizes = np.ones((len(activation_times), len(time)))

				for i_t, t in enumerate(time[:-1]):
					for i in np.arange(0, len(activation_times)):
						#clone_sizes[i, int(activation_times[i]/dT):] = np.exp(lambda_B*(time[int(activation_times[i]/dT):] - activation_times[i] ))
						tb = activation_times[i]
						Nb = clone_sizes[i, i_t]# * np.heaviside(tb-t)
						N = np.sum(clone_sizes[:, i_t]) - np.sum(clone_sizes[:, 0])
						clone_sizes[i, i_t+1] = Nb + delta_Nb(t, tb, Nb, N, lambda_B, C)*dT

				clone_sizes_f = clone_sizes[:,-1]

				#-----t_C filter-------
				# filter_C = activation_times<t_C
				filter_C = clone_sizes_f > 1

				clone_sizes_C = clone_sizes[filter_C, :]
				activation_times_C = activation_times[filter_C]
				energies_C = energies[filter_C]

				#-------Simulations-------
				Kds_C = np.exp(energies_C)
				NC_i = 1-np.array([np.product(1-1/(1+(Kds_C/((1e5*(clone_sizes_C[:,t]-1))/N_A)))) for t in np.arange(len(time))])
				NC += NC_i/N_ens
				if(i_ens%2==0):
					AX.plot(time, NC_i, color = colors_theta[i_theta], alpha = .1)

			AX.plot(time, NC, color = colors_theta[i_theta], alpha = 1)
			#-------Theory-------




my_plot_layout(ax = AX, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#AX.legend(fontsize = 30, title_fontsize = 35, title = r'$\theta$')
#AX.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#AX.set_ylim(bottom = 5e-5, top = 2e0)
#AX.set_yticks([1, 0.1, 0.01, 0.001])
#AX.set_yticklabels([1, 0.1, 0.01])
FIG.savefig('../../Figures/1_Dynamics/Ensemble/NC.pdf')





