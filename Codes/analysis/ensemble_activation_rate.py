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

N_ens = 500
N_r = 5e4
N_r = 1e5
T0 = 0
Tf = 6
Tf = 8
dT = .01
days = np.arange(0, Tf, 1)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
time_bins = time[::5]
lambda_A = 6 #days^-1
k_pr = .1 # hour^-1
k_pr = k_pr*24 #days^-1
qs = [1, 2, 3]
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e3
E_ms = -28

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'TACNSEYPNTTK'
antigen = 'TACNSEYPNTTKCGRWYC'
#antigen = 'FTSENAYCGR'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'

L=len(antigen)
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
S = np.cumsum(lambdas[:-1]*dE)
Omega = np.sum(np.exp(S)*dE)
Ks = np.exp(Es[:-1])

beta_pr = lambdas[:-1][Es[:-1]<np.log(k_pr/k_on)][-1]
E_pr = Es[:-1][Es[:-1]<np.log(k_pr/k_on)][-1]
print('E_pr:%.2f'%E_pr, 'beta_pr:%.2f'%beta_pr)

beta_r = lambdas[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
print('E_r:%.2f'%E_r, 'beta_r:%.2f'%beta_r)

#----------------------------------------------------------------
t_act = [3.8, 4.5]

energy_model = 'MJ'
for q in qs:
	beta_q = lambdas[lambdas>q][-1]
	E_q = Es[lambdas>q][-1]
	print('E^(%d):%.2f'%(q, E_q), 'beta_q:%.2f'%beta_q)
	for j, gm in enumerate(growth_models):
		fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})

		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_q-%d_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, q, j, N_ens)+energy_model
		data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		data2 = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/summary_ensemble.txt', sep = '\t', header=None)
		N_clones = np.array(data2[0])
		N_ens = len(data2[0])

		data_active = data.loc[data[1]==1]	
		activations_times = np.array(data_active[3])
		ar1, ar2 = np.histogram(activations_times, bins = time_bins)
		data_N_active_linages = np.cumsum(ar1/N_ens)

		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf)/N_A, Es, q, lambda_A, N_c, dE)
		psi = ((k_on*N_c)/(N_A))*p_a
		M_r = N_r*N_c*np.sum(Q0*p_a*dE)

		m_bar = np.array([N_r*(1-np.sum(np.exp(-(psi/lambda_A)*(np.exp(lambda_A*t)))*Q0*dE)) for t in time])

		m_bar2 = np.array([N_r*(np.sum((1-np.exp(-(np.exp(lambda_A*t)/N_A)*k_on*p_a*N_c/lambda_A))*Q0*dE)) for t in time])

		t_act = (1/lambda_A)*(np.log((lambda_A*N_A*Omega*beta_q)/(N_c*k_on*N_r)) - beta_q*E_q)
		print(t_act)
		t_act2 = (1/lambda_A)*np.log((lambda_A*N_A)/(N_c*k_on))
		

		delta_t = 0
		delta_t2 = 0
		if(E_r>E_q):
			delta_t =  q*(E_r-E_pr)/lambda_A
			t_act2 =  q*(E_r-E_q-(1/q)*np.log((N_c*k_on)/(lambda_A*N_A)))/lambda_A
		print(t_act2)
		

		t_act_data = time_bins[:-1][data_N_active_linages>1][0]
		t_act_theory = time[m_bar2>1][0]

		#Simulation
		ax.plot(time_bins[:-1], data_N_active_linages, linestyle = '-', marker = 'o', ms = 4, linewidth = 1, label = 'simulation', color = colors[j])
		#Exact analytics
		#ax.plot(time, m_bar, color = colors[j], label = growth_models[j] + ' growth', linestyle = '--', linewidth = 3, alpha = .6)
		ax.plot(time, m_bar2, color = colors[j], label = growth_models[j] + ' growth', linestyle = '-', linewidth = 3, alpha = .6)
		#approximations
		#ax.plot(time, data_N_active_linages[-10]*np.exp(lambda_A*(time))/np.exp(lambda_A*(time_bins[:-1][-10])), color = 'black', linestyle = ':', linewidth = 1, alpha = .6)
		#ax.plot(time, data_N_active_linages[-1]*np.exp((lambda_A*(np.min([beta_q, beta_r]))/q)*(time))/np.exp((lambda_A*(np.min([beta_q, beta_r]))/q)*(time_bins[:-1][-1])), color = 'black', linestyle = '--', linewidth = 1, alpha = .6)
		
		ax.vlines([t_act, t_act2, t_act_data, t_act_theory], 1e-6, N_r, color = ['black', 'green', 'darkblue', 'darkblue'], linestyle = ['-', ':', '--', '-'], linewidth = [2, 2, 2, 2], alpha = .6)
		ax.hlines(1, 0, Tf, linestyle = 'dashed', color = 'black', linewidth = 1, alpha = .6)
		
		my_plot_layout(ax = ax, xscale='linear', yscale= 'log', xlabel=r'time', ylabel = r'$\bar m(t)$', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		ax.set_xlim(left = 3, right = Tf)
		ax.set_ylim(bottom = 1e-5, top = 1e5)
		#ax.legend(title=r'$\lambda_A/\lambda_B$', fontsize = 24, title_fontsize = 24)
		fig.savefig('../../Figures/1_Dynamics/Ensemble/activation_rate_q-%d.pdf'%q)

	#ax.text(x=text_pos[i][0], y=text_pos[i][1], s = text[i], fontsize=44, color = colors_fit[i])




