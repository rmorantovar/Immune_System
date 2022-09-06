import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.special as sc

import pickle
import json
import warnings

from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
from tqdm import tqdm
from scipy.optimize import curve_fit
from matplotlib import style
from matplotlib.collections import LineCollection

# ----- CONSTANTS ------
N_A = 6.02214076e23
k_BT = 1.380649e-23*293


# --- functions ----

delta_Nb = lambda t, tb, Nb, N, lambda_B, C: lambda_B*Nb*(1-(N/C))*np.heaviside(t-tb, 1)

E_t = lambda t, theta:lambda_A*t/theta - np.log((lambda_A*N_A)/(k_on*N_c))/theta + np.log(k_pr/k_on) 

def get_clones_sizes_C(n_act, time, activation_times, lambda_B, C, dT):
	clone_sizes = np.ones((n_act, len(time)))
	for i_t, t in enumerate(time[:-1]):
		for i in np.arange(0, n_act):
		    #clone_sizes[i, int(activation_times[i]/dT):] = np.exp(lambda_B*(time[int(activation_times[i]/dT):] - activation_times[i] ))
			tb = activation_times[i]
			Nb = clone_sizes[i, i_t]# * np.heaviside(tb-t)
			N = np.sum(clone_sizes[:, i_t]) - np.sum(clone_sizes[:, 0])
			clone_sizes[i, i_t+1] = Nb + delta_Nb(t, tb, Nb, N, lambda_B, C)*dT
	return clone_sizes

def get_motif(antigen, Matrix, Text_files_path):

	if(Matrix == 'TCRen'):
		TCRen = pd.read_csv('../Input_files/TCRen_potential.csv')
		Alphabet = ["C", "S", "T", "P", "A", "G", "N", "D", "E", "Q", "H", "R", "K", "M", "I", "L", "V", "F", "Y", "W"]
		Alphabet_list = Alphabet
		TCRen_dict = dict()
		alphabet_from  = np.unique(TCRen['residue.aa.from'])
		alphabet_to  = np.unique(TCRen['residue.aa.to'])
		for i in range(len(alphabet_from)):
			TCRen_dict[alphabet_from[i]] = dict()
		for i in range(len(TCRen)):
			residue_from = TCRen['residue.aa.from'][i]
			residue_to = TCRen['residue.aa.to'][i]
			TCRen_dict[residue_from][residue_to] = TCRen['TCRen'][i]
		TCR_matrix = np.zeros((20, 20))
		for i, aa1 in enumerate(Alphabet_list):
			for j, aa2 in enumerate(Alphabet_list):
				TCR_matrix[i, j] = TCRen_dict[aa1][aa2]
		M = TCR_matrix
	if(Matrix == 'MJ2'):
		M = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
		Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
		Alphabet_list = Alphabet.tolist()

	antigen_list = [i for i in antigen]
	antigen_seq = np.array([], dtype = int)
	for i, aa in enumerate(antigen_list):
	    index = Alphabet_list.index(aa)
	    antigen_seq = np.append(antigen_seq, int(index))

	return M[:,antigen_seq]

def get_repertoire_properties(betas, Q0, Es, dE, N_r):
	beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
	E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
	Kd_r = np.exp(E_r)

	return beta_r, E_r, Kd_r

def get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on):
	E_pr = Es[:-1][Es[:-1]<np.log(k_pr/k_on)][-1]
	Kd_pr = np.exp(E_pr)
	beta_pr = betas[:-1][Es[:-1]<E_pr][-1]

	return beta_pr, E_pr, Kd_pr

def get_n_properties(betas, Q0, Es, dE, theta):

	beta_theta = betas[betas>theta][-1]
	E_theta = Es[betas>theta][-1]
	Kd_theta = np.exp(E_theta)

	return beta_theta, E_theta, Kd_theta

def hamming_distance(chaine1, chaine2):

    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_complementary_seq(sequence, Energy_Matrix):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],min(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def find_complementary_seq_2(sequence, Energy_Matrix):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],max(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def calculate_energy(Energy_Matrix, seq1, seq2):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq1 = list(seq1)
	list_seq2 = list(seq2)
	Energy = 0
	for i in range(9):
		pos_i = np.where(np.isin(Alphabet,list_seq1[i]))[0][0]
		pos_j = np.where(np.isin(Alphabet,list_seq2[i]))[0][0]
		Energy += M[pos_i][pos_j]
	return Energy

def my_linear_func(x, a, b):

    return a + b*x

def my_quadratic_func(x, a, b, c):

    return np.log(a)+np.log(np.sqrt(-b)) + b*(x-c)**2

def binding_affinity(x, a, b):
	
    return 1/(1+((np.exp(a+b*x))))

def Z_PWM(PWM, T):
    Z = 1
    for i in range(len(PWM[0,:])):
        Z_i = 0
        for j in range(len(PWM[:,0])):
            Z_i = Z_i + np.exp((-PWM[j, i]/T))
        Z = Z*Z_i
    return Z

def Z_PWMj(PWM, beta):
    Z = 1
    for i in range(len(PWM[0,:])):
        Z_i = 0.j
        for j in range(len(PWM[:,0])):
            Z_i = Z_i + np.exp(-PWM[j, i]*beta)
        Z = Z*Z_i
    return Z

def Z_PWM_integral(T, lamda, k):
    return np.exp(-lamda*min_E)*(np.exp(max_E*(lamda-(1/T)))-np.exp(min_E*(lamda-(1/T))))/(lamda-(1/T))*k

def Z_PWM_integral2(T, lamda):
    return 2*np.exp(-lamda*min_E)*(np.exp(avg_E*(lamda-(1/T)))-np.exp(min_E*(lamda-(1/T))))/(lamda-(1/T))

def calculate_Q0(Tmin, Tmax, n_T, E_matrix, E_ms, L):

	Ts = np.linspace(Tmin, Tmax, n_T)
	lambdas = 1/Ts[:-1]
	F_PWM = -Ts*np.log(Z_PWM(E_matrix, Ts))
	Es = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts)) + E_ms
	dE = np.diff(Es)
	
	S = np.cumsum(lambdas[:-1]*dE)

	Omega = 20**L
	Omega = np.sum(np.exp(S)*dE)
	Omega = 2*np.sum(np.exp(S)*dE)
	
	Q0 = np.exp(S)/Omega

	return Es, dE, Q0, lambdas

def calculate_QR(Q0, k_on, k_act, rho_A, Es, q, lambda_A, N_c, dE):

	p_a = (1/(1 + (k_on*np.exp(Es[:-1])/k_act)**q) )
	u_on = rho_A*k_on
	R = 1-np.exp(-u_on*p_a*N_c/lambda_A)
	QR = Q0*R

	return u_on, p_a, R, QR

def get_t_act(time, N_r, Q0, k_on, k_pr, lambda_A, Es, dE, theta, N_c):
	u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, theta, lambda_A, N_c, dE)
	m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])
	t_act = time[m_bar>1][0] # Activation time	
	return t_act

def P_e_gaussian(avg_E, var_E, Es):

    return (2*np.pi*var_E)**(-0.5)*np.exp(-(Es-avg_E)**2/(2*var_E))

def P_min_e(N, avg_E, var_E, Es, dE):

    return (N*(1-np.cumsum(P_e_gaussian(avg_E, var_E, Es)*dE))**(N-1)*(P_e_gaussian(avg_E, var_E, Es)))

def apply_filter_C(clone_sizes, activation_times, energies, lim_size):
	filter_C = clone_sizes[:, -1] > lim_size
	n_C = np.sum(filter_C)
	clone_sizes_C = clone_sizes[filter_C, :]
	activation_times_C = activation_times[filter_C]
	energies_C = energies[filter_C]

	return clone_sizes_C, activation_times_C, energies_C, filter_C, n_C

def my_plot_layout(ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 24, xlabel = '', ylabel = '', title = '', x_fontsize=24, y_fontsize = 24, t_fontsize = 24):
    ax.tick_params(labelsize = ticks_labelsize)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize = x_fontsize)
    ax.set_ylabel(ylabel, fontsize = y_fontsize)
    ax.set_title(title, fontsize = t_fontsize)


#----------------- Plots -----------------

def plot_energy_matrix(Energy_Matrix, Alphabet, title, ax):

	M = Energy_Matrix
	
	Alphabet = Alphabet

	sns.heatmap(np.flip(M, axis = 0), ax = ax, cmap=plt.cm.viridis, center = 0, cbar = True)
	ax.set_title(title, fontsize = 22)
	ax.tick_params(labelsize = 20)
	ax.set_xticklabels(Alphabet)
	ax.set_yticklabels(np.flip(Alphabet));
	cbar = ax.collections[0].colorbar
	cbar.ax.tick_params(labelsize=18)

def plot_PWM(PWM, Alphabet, sequence, title, ax):

	M = PWM
	
	Alphabet = Alphabet

	sns.heatmap(np.flip(M, axis = 0), ax = ax, cmap=plt.cm.viridis, center = 0, cbar = True)
	ax.set_title(title, fontsize = 22)
	ax.tick_params(labelsize = 20)
	ax.set_xticklabels(sequence)
	ax.set_yticklabels(np.flip(Alphabet));
	cbar = ax.collections[0].colorbar
	cbar.ax.tick_params(labelsize=18)

def plot_histogram_hamming_distance(Sequences, ax):

	distances = np.array([i.hamming_distance for i in Sequences])
	data_distances = np.histogram(distances, bins=range(int(max(distances))))

	#ax.plot(data_distances[1][0:-1], scipy.special.comb(9, data_distances[1][0:-1]), linewidth = 4 , label = 'Binary')
	ax.plot(data_distances[1][0:-1], sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), color = 'steelblue', linewidth = 4 , label = '20-Alphabet')
	#ax.plot(data_distances[1][0:-1], scipy.special.comb(9, data_distances[1][0:-1])*((4-1)**data_distances[1][0:-1]), linewidth = 4 , label = '4-Alphabet')

	#ax.plot(data_distances[1][0:-1], np.exp(4*data_distances[1][0:-1]), linewidth = 4, label = r'$e^{\lambda r}$')
	ax.plot(data_distances[1][0:-1], data_distances[0], linewidth = 4, label = 'Data', linestyle = '', marker = 'o')

	ax.set_yscale('log')
	#ax.set_ylim(1,1e8)
	ax.set_xlabel(r'Hamming Distance $r$', fontsize = 20)
	ax.set_ylabel(r'', fontsize = 20)
	ax.tick_params(labelsize = 20)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(np.concatenate(([],handles)),np.concatenate(([],labels)), loc = 2, fontsize = 20)

	return distances

def plot_histogram_energy(Sequences, normalization, bins, color, density, n_seq, ax):

	energies = np.array([i.energy for i in Sequences])
	data_energies = np.histogram(energies, bins=bins, density = density)

	ax.plot(data_energies[1][0:-1], data_energies[0]/normalization, linewidth = 4, color = color, linestyle = '', marker = 'o', label = 'N=%.e'%(n_seq))
	ax.set_yscale('log')
	#ax.set_ylim(1,1e10)
	ax.set_xlabel(r'Energy $r$', fontsize = 20)
	ax.set_ylabel(r'', fontsize = 20)
	ax.tick_params(labelsize = 20)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(np.concatenate(([],handles)),np.concatenate(([],labels)), loc = 0, fontsize = 20)

	return energies, data_energies

def plot_scatter_hamming_distance_energy(distances, energies, color, ax):

	ax.scatter(distances, energies, color = color, s = 15, marker = '^')
	ax.hlines(energies[0],0,9, linestyle = 'dashed', label = r'Master Seq. energy $\epsilon_0$')
	#ax.set_yscale('log')
	#ax.set_ylim(1,1e10)
	ax.set_xlabel(r'Hamming distance $d$', fontsize = 20)
	ax.set_ylabel(r'Energy $\epsilon$', fontsize = 20)
	ax.tick_params(labelsize = 20)
	ax.legend(loc = 0, fontsize = 20)

def plot_histogram_energy_subsampling(bins, n_linages, Sequences, sub_energies, fig, ax):

	energies, data_energies = plot_histogram_energy(Sequences = Sequences, bins = bins, ax = ax)

	#____________ fit and plot gaussian function
	popt, pcov = curve_fit(my_quadratic_func, data_energies[1][0:-1][data_energies[0]!=0], np.log(data_energies[0][data_energies[0]!=0]) , p0 = (5e5/np.sqrt(np.pi), -0.04, -30))
	print(r'integral a*$\pi$:',popt[0]*np.pi)
	r_array = np.linspace(np.min(energies)-5, np.max(energies)+5, 1000)
	ax.plot(r_array, np.exp(my_quadratic_func(r_array, *popt)), linestyle = '--', linewidth = 4, color = 'steelblue', alpha = 0.4)
	r_array2 = np.linspace(-80, 50, 10000)


	#____________ Plot histogram of sub_energies
	data_energies = np.histogram(sub_energies, bins=20, density = False)
	ax.plot(data_energies[1][0:-1], data_energies[0]/(2000), linewidth = 4, color = 'indigo', label = 'After samplig', linestyle = '', marker = 'o')
	ax.plot(r_array, (2e2/5e5)*np.exp(my_quadratic_func(r_array, *popt)), linestyle = '--', linewidth = 4, color = 'indigo', alpha = 0.4)
	ax.plot(r_array, (2e2/5e5)*popt[0]*np.sqrt(-popt[1])*(1+popt[1]*(r_array-popt[2])**2) )
	ax.set_ylim(1,8e4)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(np.concatenate(([],handles)),np.concatenate(([],labels)), loc = 2, fontsize = 20)

	#____________ Create inset with the integral of the gaussian function  
	left, bottom, width, height = [0.65, 0.67, 0.25, 0.2]
	ax2 = fig.add_axes([left, bottom, width, height])
	ax2.plot(r_array, np.cumsum((2e2/5e5)*np.exp(my_quadratic_func(r_array, *popt))*(np.max(energies)-np.min(energies))/1000), linewidth = 2, color = 'indigo')
	ax2.set_xlabel(r'Energy $r$', fontsize = 15)
	ax2.set_ylabel(r'', fontsize = 15)
	ax2.tick_params(labelsize = 15)
	ax2.yaxis.tick_right()
	ax2.set_yscale('log')

	return popt, pcov

#----------------- Generate files -----------------

		
def print_raw_file(Sequences, filename):

	file = open(filename, 'w+')
	for i in range(len(Sequences)):
		np.savetxt(file, np.array([Sequences[i].parent_id, Sequences[i].id]), fmt = '%d', delimiter = ' ', newline = ' ')
		file.write("\n")
	file.close()

def generate_newick_format(filename):

	file  = np.loadtxt(filename, dtype = int)
	n_f = '0()'
	print(file)


#----------------- Plots for ensemble averages -----------------

def plot_activation_rate_ensemble_deterministic(beta, b, nu, gamma, T, initial_time, eo, n_linages, n_left_tail, rho_min, rho_max, dt, popt, energies, comment, gaussian, exponential, ax):

	N_A = 6.02214076e23
	to = (eo+np.log(N_A))/beta
	#____________ Read and plot the activation of linages as function of antigen concetration
	t_new = np.linspace(initial_time, T, int((T-initial_time)/dt))
	activation_time_series = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System//Text_files/ensemble_deterministic_model_activation_time_series_"+comment+".pkl", "rb" ) )
	activation_time_series_var = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System//Text_files/ensemble_deterministic_model_activation_time_series_var_"+comment+".pkl", "rb" ) )
	ax.plot(np.exp(beta*t_new[1:][::50][np.where(activation_time_series[1:][::50]!=0)])/N_A, activation_time_series[1:][::50][np.where(activation_time_series[1:][::50]!=0)], linestyle = '', marker = '.', ms = 15, linewidth = 4, color = 'indigo', label = 'simulation')
	ax.fill_between(np.exp(beta*t_new[1:][::10][np.where(activation_time_series[1:][::10]!=0)])/N_A , activation_time_series[1:][::10][np.where(activation_time_series[1:][::10]!=0)] - np.sqrt(activation_time_series_var[1:][::10][np.where(activation_time_series[1:][::10]!=0)]), activation_time_series[1:][::10][np.where(activation_time_series[1:][::10]!=0)] + np.sqrt(activation_time_series_var[1:][::10][np.where(activation_time_series[1:][::10]!=0)]), linewidth = 4, color = 'indigo', alpha = 0.2)

	#____________ Plot the gaussian integral
	if(gaussian):
		r_array = np.linspace(np.min(energies), np.max(np.log(np.exp(t_new[1:][::10][np.where(activation_time_series[1:][::10]!=0)])/N_A)), 5000)
		ax.plot(np.exp(r_array), np.cumsum((2e2/5e5)*np.exp(my_quadratic_func(r_array, *popt))*(np.max(energies)-np.min(energies))/5000), linestyle = '--', ms = 15, linewidth = 4, color = 'violet', label = 'Gaussian integral')

	#____________ Plot the exponential integral
	if(exponential):
		rho_new = np.logspace(np.log10(rho_min), np.log10(rho_max), 50)
		alpha = n_linages/n_left_tail
		ax.plot(rho_new, (alpha/(b*beta))*(np.exp(-beta*b*to)*(rho_new*N_A)**(b)-1), linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = 'exponential model', color = 'indigo')


	ax.set_xlabel(r'Antigen concentration $[M]$', fontsize = 20)
	ax.set_ylabel(r'Activated linages', fontsize = 20)
	ax.tick_params(labelsize = 22)
	ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_ylim(0.1, max(activation_time_series[1:]*1.5))
	ax.legend(loc = 0, fontsize = 20)

def plot_size_distribution_ensemble_deterministic(beta, b, nu, gamma, T, eo, dt, n_bins, density, popt, comment, gaussian, exponential, ax):

	N_A = 6.02214076e23
	to = (eo+np.log(N_A))/beta

	#____________ Read and plot the distribution of clone sizes
	activated_linages_size = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_linage_sizes_"+comment+".pkl", "rb" ) )
	bins = np.logspace(0, np.log10(np.max(activated_linages_size)), n_bins)
	data_activated_linages_log = np.histogram(activated_linages_size, bins = bins, density = True)
	
	#Distribution
	ax.plot(data_activated_linages_log[1][:-1], 1 - np.cumsum(data_activated_linages_log[0]*np.diff(data_activated_linages_log[1])), marker = '.', ms = 15, linestyle = '', linewidth = 3, color = 'indigo', label = 'Simulation')
	
	n_array = np.logspace(0,np.log10(np.max(activated_linages_size)), 50)
	#____________ Plot the gaussian integral
	if(gaussian):
		ax.plot(n_array, 1-((len(activated_linages_size)/2e2)*(2e2/5e5)*np.exp(my_quadratic_func(np.log((np.exp(T)/N_A)/(n_array)), *popt))), linestyle = '--', ms = 20, linewidth = 4, color = 'violet', label = 'Gaussian model')

	#____________ Plot the exponential integral
	if(exponential):
		ax.plot(n_array, ((1/(beta*b))*(np.exp(beta*b*(T-to))*n_array**(-(b*beta)/nu)-1))/((((1/(beta*b))*(np.exp(beta*b*(T-to))*n_array**(-b*beta)-1)))[0]), linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = 'exponential model', color = 'indigo')

	ax.set_xlabel(r'Clone size $n_i$', fontsize = 20)
	ax.set_ylabel(r'counts', fontsize = 20)
	ax.tick_params(labelsize = 22)
	ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_xlim(.9,1000)
	ax.legend(loc = 0, fontsize = 20)

def plot_N_total_ensemble_deterministic(beta, b, nu, gamma, T, initial_time, eo, n_linages, n_left_tail, dt, popt, comment, gaussian, exponential, ax):

	N_A = 6.02214076e23
	to = (eo+np.log(N_A))/beta
	t_new = np.linspace(initial_time, T, int((T-initial_time)/dt))

	#____________ Read and plot the distribution of clone sizes
	N_total_sim = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_N_total_"+comment+".pkl", "rb" ) )

	#Total linage size
	ax.plot(t_new[::50], N_total_sim[::50], marker = '.', ms = 20, linestyle = '', linewidth = 3, color = 'indigo', label = 'Simulation')

	#____________ Plot the gaussian integral
	if(gaussian):
		N_total = np.array([np.cumsum(np.exp(t-np.linspace(0,t, 100))*(2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt))*(t/100))[-1] + 200 - np.cumsum((2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt))*(t/100))[-1] for t in t_new])
		ax.plot(t_new, N_total, linestyle = '--', linewidth = 4, color = 'violet', alpha = 0.4, label = 'gaussian model')

	#____________ Plot the exponential integral
	if(exponential):
		tau = np.exp(t_new-to)
		alpha = n_linages/n_left_tail
		N_total_exp = alpha*((tau**(beta*b)-tau**(nu))/(beta*b-nu)) + n_linages - (alpha/(beta*b))*(tau**(beta*b)-1)
		ax.plot(t_new, N_total_exp, linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = 'exponential model', color = 'indigo')

	ax.set_xlabel(r'Time $t$', fontsize = 20)
	ax.set_ylabel(r'size $N_{total}$', fontsize = 20)
	ax.tick_params(labelsize = 22)
	#ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_xlim(T/2,T)
	ax.legend(loc = 0, fontsize = 20)

def plot_entropy_ensemble_deterministic(beta, b, nu, gamma, T, initial_time, eo, n_linages, n_left_tail, dt, popt, comment, gaussian, exponential, ax):

	N_A = 6.02214076e23

	to = (eo+np.log(N_A))/beta
	t_new = np.linspace(initial_time, T, int((T-initial_time)/dt))
	#____________ Read and plot the distribution of clone sizes
	entropy = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_entropy_"+comment+".pkl", "rb" ) )

	#Total linage size
	ax.plot(t_new[::50], entropy[::50], marker = '.', ms = 20, linestyle = '', linewidth = 3, color = 'indigo', label = 'Simulation')

	#____________ Plot the gaussian integral
	if(gaussian):
		N_total = np.array([np.cumsum(np.exp(t-np.linspace(0,t, 100))*(2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt))*(t/100))[-1] + 200 - np.cumsum((2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt))*(t/100))[-1] for t in t_new])
		N_alpha_log_N_alpha = np.array([np.cumsum((np.exp(t-np.linspace(0,t, 100))*(t-np.linspace(0,t, 100))*(2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt)))*(t/100))[-1] for t in t_new])
		ax.plot(t_new, (-1/N_total)*(N_alpha_log_N_alpha) + np.log(N_total), linestyle = '--', linewidth = 4, color = 'violet', alpha = 0.4, label = 'gaussian model')

	#____________ Plot the exponential integral
	if(exponential):
		tau = np.exp(t_new-to)
		alpha = n_linages/n_left_tail
		N_total_exp = alpha*((tau**(beta*b)-tau**(nu))/(beta*b-nu)) + n_linages - (alpha/(beta*b))*(tau**(beta*b)-1)
		Entropy_exp = ((nu*alpha)/(N_total_exp*(beta*b-nu))) * (tau**(nu)*np.log(tau)  - ((tau**(beta*b)-tau**(nu))/(beta*b-nu))) + ((np.log(N_total_exp)*1)/(1))
		ax.plot(t_new, Entropy_exp, linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = 'exponential model', color = 'indigo')
		ax.hlines(1/(1-beta*b), initial_time, T, linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = r'$S_{\infty}$')

	ax.set_xlabel(r'Time $t$', fontsize = 20)
	ax.set_ylabel(r'Entropy $S_{Linages}$', fontsize = 20)
	ax.tick_params(labelsize = 22)
	#ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_xlim(T/2,T)
	ax.legend(loc = 0, fontsize = 20)
