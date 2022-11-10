import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")
my_red = np.array((228,75,41))/256.
my_purple = np.array((125,64,119))/256.
my_purple2 = np.array((116,97,164))/256.
my_green = np.array((125,165,38))/256.
my_blue = np.array((76,109,166))/256.
my_gold = np.array((215,139,45))/256.
my_brown = np.array((182,90,36))/256.
my_blue2 = np.array((80,141,188))/256.
my_yellow = np.array((246,181,56))/256.
my_green2 = np.array((158,248,72))/256.
my_olive = np.array((108,114,44))/256
my_cyan = 'tab:cyan'
Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
#N_ensss = [[200], [501, 502, 503, 504, 505, 506, 507, 508, 509, 400, 300, 200, 100, 50], [200, 150, 100], [200, 100], [200]]
#N_ensss = [[400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 501, 502, 503, 504, 505, 506, 507, 508, 509]] #for p=2.5
N_ensss = [[50, 100, 150, 200] + [400+i for i in range(1, 10)] + [500+i for i in range(0, 91)] + [1000+i for i in range(1, 51)]] #for p=3

N_r = 2e8

T0 = 3
Tf = 12
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1
lambda_B = lambda_A/2
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5
#N_c = 1e5
E_ms = -27.63
C = 3e4
AA = 1

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

kappas = [2.5]
kappas = [3.0]

antigen_color = my_yellow

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_red, my_green, my_blue2, my_gold, my_purple])
color_list = np.array([my_green])

colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'
L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path)
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_K_distribution, ax_K_distribution = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_distribution2, ax_K_distribution2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

K_final_best_renorm = []
counter_total = 0
for i_kappa, kappa in enumerate(kappas):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0] 
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	beta_act = np.min([beta_r, beta_kappa])

	N_enss = N_ensss[i_kappa]

	K_all = np.zeros_like(time)
	K_final_all = []
	Counter_all = 0

	for N_ens in N_enss:
		print('N_ens = %d'%N_ens)

		#-----------------Loading data----------------------------
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
		#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		data, return_data_type = get_data_ensemble_K_elite(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
		
		if(return_data_type):
			K_final_all = data[0]
			Counter_all = data[1]
		else:

			for i_ens in tqdm(np.arange(N_ens)):
				data_i = data.loc[data[4]==i_ens]
				data_active = data_i.loc[data_i[1]==1]
				t_act_data = np.min(data_active[3])
				data_active = data_active.loc[data_active[3]<(t_act_data+1.0+0.1*(kappa-1))]
				
				data_active_all = data_active#.loc[data_active[3]<(t_act_theory)]
				activation_times_all = np.array(data_active_all[3])
				energies_all = np.array(data_active_all[0])

				#---------------------------- B cell linages ----------------------
				clone_sizes_all = get_clones_sizes_C(len(activation_times_all), time, activation_times_all, lambda_B, C, dT)
				#--------------------------t_C filter-------------------------
				lim_size = 2
				clone_sizes_C_all, activation_times_C_all, energies_C_all, filter_C_all, n_C_all = apply_filter_C(clone_sizes_all, activation_times_all, energies_all, lim_size)
				#-------Simulations-------
				if(len(energies_C_all)>0):
					Kds_C_all = np.exp(energies_C_all)
					Avidities = np.divide(((clone_sizes_C_all-1).T)/1, Kds_C_all).T
					K_i_all = np.sum(Avidities, axis = 0)
										
					K_final_best_renorm.append(((C/1)/np.min(Kds_C_all)))

					if(np.sum(K_i_all)!=0):
						#K_all += K_i_all
						K_final_all.append(K_i_all[-1])
						Counter_all+=1

			f = open(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/processed_data_K_elite.pkl', 'wb')
			pickle.dump([K_final_all, Counter_all], f, pickle.HIGHEST_PROTOCOL)	

	counter_total+=Counter_all
	Counter_all = 1
	normalization_all = 1
	print('%.2e'%np.max(K_final_all))
	#K_data_all = np.histogram(np.log10(np.array(K_final_all)), bins = np.linspace(10.5, 14, 34), density = False)
	K_data_all = np.histogram(np.log10(np.array(K_final_all)), bins = 'auto', density = True)
	
	ax_K_distribution.plot(10**(K_data_all[1][:-1]), K_data_all[0]/Counter_all, color = 'limegreen', linestyle='', marker = 'D', linewidth = 2, ms = 5, label = r'$p=%.1f$'%(kappa))

	ax_K_distribution2.plot(10**(K_data_all[1][:-1]), 1-np.cumsum(K_data_all[0]*np.diff(K_data_all[1])/Counter_all), color = 'limegreen', linestyle='-', marker = '', linewidth = 5, ms = 8, label = r'$p=%.1f$'%(kappa), zorder = 20)
	
	ax_K_distribution2.vlines(8e11, 1e-5, 1, color = my_green, linestyle=':', linewidth = 4)		
print(counter_total)

# Printing K from Gumbel
Nb = C
#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
K_array = ((Nb/1)/Kds)
p_K = P_min_e_Q0(N_r, Q0, dE)#/K_array**2*(Nb/1)
p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
ax_K_distribution.plot(((np.flip(K_array[:-1]))/1), np.flip(p_K[:-1]), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 5, alpha = .8, label = 'Gumbel')
#ax_K_distribution2.plot(((np.flip(K_array[:-1]))/1), 1-np.cumsum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array)))), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 4, alpha = .4, label = 'Gumbel')

K_array_tail = 10**np.linspace(9, 14.5, 50)
exponent_tail  = beta_r+1
#fit_tail = np.exp(-exponent_tail*(K_array_tail))/np.exp(-exponent_tail*(12.3))
fit_tail =(K_array_tail**(-exponent_tail))/((10**13.2)**(-exponent_tail))
#fit_tail *= (1-np.cumsum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array)))))[np.log10(np.flip(K_array)[:-1]/1)<13.2][-1]
fit_tail *= (1-np.cumsum(K_data_all[0]*np.diff(K_data_all[1])))[((K_data_all[1][:-1]))<13.2][-1]
#print((1-np.cumsum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array)))))[(np.flip(K_array)[:-1]/normalization_all)<27.5][-1])
ax_K_distribution2.plot(K_array_tail, fit_tail, linewidth = 4, color = 'black', linestyle = 'dashed')

my_plot_layout(ax = ax_K_distribution, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K_distribution.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_K_distribution.set_ylim(bottom = 5e-6, top = 3)
ax_K_distribution.set_xlim(left = 5*10**(10), right = 10**(14.5))
#ax_K_distribution.set_xticks([])
#ax_K_distribution.set_yticks([])
#ax_K_distribution.set_yticklabels([1, 0.1, 0.01])
fig_K_distribution.savefig('../../Figures/1_Dynamics/Ensemble/K_elite_P_'+energy_model+'_p-%.1f.pdf'%(kappa))

my_plot_layout(ax = ax_K_distribution2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K_distribution2.legend(fontsize = 28, title_fontsize = 30, loc = 0)
ax_K_distribution2.set_ylim(bottom = 5e-5, top = 1)
ax_K_distribution2.set_xlim(left = 6*10**(11), right = 6*10**(13))
#ax_K_distribution2.set_xticks([])
#ax_K_distribution2.set_yticks([])
#ax_K_distribution2.set_yticklabels([1, 0.1, 0.01])
fig_K_distribution2.savefig('../../Figures/1_Dynamics/Ensemble/K_elite_F_'+energy_model+'_p-%.1f.pdf'%(kappa))



