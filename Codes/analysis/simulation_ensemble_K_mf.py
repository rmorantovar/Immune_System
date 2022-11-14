import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 500
N_r = 2e8
transparencies_p = [.8, 1, .8, .8, .8]

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


kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2, 3, 4]#, 5]
#kappas = [3]

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
my_cyan = 'tab:cyan'

antigen_color = my_yellow/256.

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_red, my_blue2, my_green, my_gold, my_brown])
#color_list = np.array([my_green, my_blue2, my_gold])

colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

# antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
# antigen = 'FMLFMAVFVMTSWYC'
# antigen = 'FTSENAYCGR'
# antigen = 'TACNSEYPNTTK'
antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'
#antigen = 'TACNSEYPNTTKCGRWYC'
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
print('N_r = %.e'%N_r)
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_K_mf, ax_K_mf = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
print('--------')
print('--- Processing all clones ---')
print('--------')
max_potency_simulations = dict()
max_potency_simulations_std = dict()
max_potency_theory = dict()

for i_kappa, kappa in enumerate(kappas):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0] 
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)

	if kappa==1:
		t_act_1 = t_act_theory
	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
	#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
	#data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
	data, return_data_type = get_data_ensemble_K_mf(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)

	if(return_data_type):
		K_final = data[0]
		Counter = data[1]
	else:
		K_final = []
		Counter = 0
		for i_ens in tqdm(np.arange(N_ens)):
			data_i = data.loc[data[4]==i_ens]
			data_active = data_i.loc[data_i[1]==1]
			t_act_data = np.min(data_active[3])
			data_active = data_active.loc[data_active[3]<(t_act_data+1.0+0.1*(kappa-1))]
			activation_times = np.array(data_active[3])
			energies  = np.array(data_active[0])

			#---------------------------- B cell linages ----------------------
			clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

			#--------------------------t_C filter-------------------------
			lim_size = 2
			clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
			#-------Simulations-------
			if(len(energies_C)>0):
				Kds_C = np.exp(energies_C)

				K_i = [np.sum(((clone_sizes_C[:,t]-1)/1)/Kds_C) for t in np.arange(len(time))]

				if(np.sum(K_i)!=0):
					K_final.append(K_i[-1])
					Counter+=1

				#if(i_ens%1==0):
				#	ax_K.plot(time, K_i, color = colors_kappa[i_kappa], alpha = .1, linewidth = 1)

		f = open(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/processed_data_K_mf.pkl', 'wb')
		pickle.dump([K_final, Counter], f, pickle.HIGHEST_PROTOCOL)	


	normalization = 1
	if(kappa==1):
		normalization = 1# K[-1]

	if(kappa==1):
		# Printing K from Gumbel
		Nb = C
		#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
		K_array = ((Nb/1)/Kds)
		p_K = P_min_e_Q0(N_r, Q0, dE)#/K_array**2*(Nb/1)
		p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
		Kd_r_renorm = Kds[(P_min_e_Q0(N_r, Q0, dE)/Kds)==np.max(P_min_e_Q0(N_r, Q0, dE)/Kds)]
		ax_K_mf.hlines(0, .9, 4.1, linewidth = 2, color = 'black', linestyle = 'dashed')
		QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_1))/N_A, Es, kappa, lambda_A, N_c, dE)[3]
		Kd_1_renorm = Kds[(QR/1)==np.max(QR/1)]
		ax_K_mf.hlines(np.log10((C*1)/Kd_1_renorm/(C/Kd_r_renorm)), .9, 4.1, linewidth = 2, color = 'black', linestyle = 'dashed')
		Kd_1_renorm = Kds[(QR/Kds)==np.max(QR/Kds)]
		ax_K_mf.hlines(np.log10((C*1)/Kd_1_renorm/(C/Kd_r_renorm)), .9, 4.1, linewidth = 2, color = 'black', linestyle = 'dashed')
		

	max_potency_simulations[kappa] = np.mean(np.array(K_final)/(C/Kd_r_renorm))
	max_potency_simulations_std[kappa] = np.sqrt(np.var((np.array(K_final)/(C/Kd_r_renorm))))

ax_K_mf.plot(kappas, np.log10(np.array(list(max_potency_simulations.values()))), color = my_purple2, linestyle = '', marker = 'D', linewidth = 3, label = 'Total', ms = 10, alpha = 1)
ax_K_mf.errorbar(x=kappas, y=np.log10(np.array(list(max_potency_simulations.values()))), yerr = 1.8*np.log10(np.array(list(max_potency_simulations_std.values()))), ls = 'none', color = my_purple2, linewidth = 2, alpha = .8)


#-------------------------# 
print('--------')
print('--- Processing largest clone ---')
print('--------')
N_enss = [500, 500, 500, 500]
max_potency_simulations2 = dict()
max_potency_simulations_std2 = dict()
max_potency_theory2 = dict()

for i_kappa, kappa in enumerate(kappas):
	N_ens = N_enss[i_kappa]
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0] 
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	beta_act = np.min([beta_r, beta_kappa])

	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
	#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
	#data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
	data, return_data_type = get_data_ensemble_K_largest(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)

	if(return_data_type):
		final_biggest = data[0]
		final_biggest_affinity = data[1]
	else:	
		final_biggest = []
		final_biggest_affinity = []

		for i_ens in tqdm(np.arange(N_ens)):
			data_i = data.loc[data[4]==i_ens]
			data_active = data_i.loc[data_i[1]==1]
			t_act_data = np.min(data_active[3])
			data_active = data_active.loc[data_active[3]<(t_act_data+1.0+.1*(kappa-1))]
			
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
				
				final_biggest_affinity.append(energies_C_all[clone_sizes_C_all[:,-1]==np.max(clone_sizes_C_all[:,-1])][0])
				final_biggest.append(np.max(clone_sizes_C_all[:,-1]))

		f = open(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/processed_data_K_largest.pkl', 'wb')
		pickle.dump([final_biggest, final_biggest_affinity], f, pickle.HIGHEST_PROTOCOL)	
		
	max_potency_simulations2[kappa] = np.mean((np.array(final_biggest)/np.exp(final_biggest_affinity))/(C/Kd_r_renorm))
	max_potency_simulations_std2[kappa] = np.sqrt(np.var(((np.array(final_biggest)/np.exp(final_biggest_affinity))/(C/Kd_r_renorm))))

ax_K_mf.plot(kappas, np.log10(np.array(list(max_potency_simulations2.values()))), color = my_purple, linestyle = '', marker = '*', linewidth = 3, label = 'Largest', ms = 14, alpha = 1)
#ax_K_mf.errorbar(x=kappas, y=np.log(np.array(list(max_potency_simulations2.values()))), yerr = 1.8*np.log(np.array(list(max_potency_simulations_std2.values()))), ls = 'none', color = my_purple2)


print('--------')
print('--- Theory ---')
print('--------')
E_0_integral = np.log(Kd_r)
E_0_integral = E_ms
Nb_func = lambda C, tf, tb, lambda_B : C*np.exp(lambda_B*(tf-tb))/(C-1+np.exp(lambda_B*(tf-tb)))

kappas_theory = np.linspace(1, 4.1, 10)
for kappa in tqdm(kappas_theory):
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	E_0_integral = E_kappa
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>=1][0] 
	QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_theory))/N_A, Es, kappa, lambda_A, N_c, dE)[3]
	#numerator = np.sum(dE[:]*QR[:]*np.exp(-lambda_B*kappa/lambda_A*Es[:-1][:] - Es[:-1][:]))
	#denominator = np.sum(dE[:]*QR[:]*np.exp(-lambda_B*kappa/lambda_A*Es[:-1][:]))
	#numerator = np.sum(dE[Es[:-1]>np.log(Kd_r)]*QR[Es[:-1]>np.log(Kd_r)]*np.exp(-lambda_B*kappa/lambda_A*Es[:-1][Es[:-1]>np.log(Kd_r)] - Es[:-1][Es[:-1]>np.log(Kd_r)]))
	#denominator = np.sum(dE[Es[:-1]>np.log(Kd_r)]*QR[Es[:-1]>np.log(Kd_r)]*np.exp(-lambda_B*kappa/lambda_A*Es[:-1][Es[:-1]>np.log(Kd_r)]))
	if (kappa<beta_r):
		numerator = np.sum(dE[Es[:-1]>E_0_integral]*QR[Es[:-1]>E_0_integral]*Nb_func(C, Tf, t_E(Es[:-1][Es[:-1]>E_0_integral], kappa, lambda_A, k_on, N_c, k_pr), lambda_B)*np.exp(-Es[:-1][Es[:-1]>E_0_integral]))
		denominator = np.sum(dE[Es[:-1]>E_0_integral]*QR[Es[:-1]>E_0_integral]*Nb_func(C, Tf, t_E(Es[:-1][Es[:-1]>E_0_integral], kappa, lambda_A, k_on, N_c, k_pr), lambda_B))
		#Es[:-1]>E_r
		K = (C*(numerator/denominator))
		#print(K)
		if(kappa == 1):
			normalization_theory = 1
		max_potency_theory[kappa] = K - normalization_theory
	if (kappa>=beta_r):
		numerator = np.sum(dE[Es[:-1]>E_0_integral]*QR[Es[:-1]>E_0_integral]*Nb_func(C, Tf, t_E(Es[:-1][Es[:-1]>E_0_integral], kappa, lambda_A, k_on, N_c, k_pr), lambda_B)*np.exp(-Es[:-1][Es[:-1]>E_0_integral]))
		denominator = np.sum(dE[Es[:-1]>E_0_integral]*QR[Es[:-1]>E_0_integral]*Nb_func(C, Tf, t_E(Es[:-1][Es[:-1]>E_0_integral], kappa, lambda_A, k_on, N_c, k_pr), lambda_B))
		#Es[:-1]>E_r
		K = (C*(numerator/denominator))
		#print(K)
		if(kappa == 1):
			normalization_theory = 1
		max_potency_theory[kappa] = K - normalization_theory

ax_K_mf.plot(kappas_theory, np.log10(np.array(list(max_potency_theory.values()))/(C/Kd_r_renorm)), color = my_purple2, linestyle = '-', marker = '', linewidth = 3, label = 'theory total', alpha = .8)

print('--------')
for kappa in tqdm(kappas_theory):
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	E_0_integral = E_kappa
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>=1][0] 
	QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_theory))/N_A, Es, kappa, lambda_A, N_c, dE)[3]
	#numerator = np.sum(dE[:]*QR[:]*np.exp(-lambda_B*kappa/lambda_A*Es[:-1][:] - Es[:-1][:]))
	#denominator = np.sum(dE[:]*QR[:]*np.exp(-lambda_B*kappa/lambda_A*Es[:-1][:]))
	numerator = Nb_func(C, Tf, t_E(Es[:-1][QR==np.max(QR)], kappa, lambda_A, k_on, N_c, k_pr), lambda_B)*np.exp(-Es[:-1][QR==np.max(QR)])
	denominator = Nb_func(C, Tf, t_E(Es[:-1][QR==np.max(QR)], kappa, lambda_A, k_on, N_c, k_pr), lambda_B)
	denominator = np.sum(dE[Es[:-1]>E_0_integral]*QR[Es[:-1]>E_0_integral]*Nb_func(C, Tf, t_E(Es[:-1][Es[:-1]>E_0_integral], kappa, lambda_A, k_on, N_c, k_pr), lambda_B))
	#Es[:-1]>E_r
	K = (1e-4*(numerator/denominator))
	#print(K)
	if(kappa == 1):
		normalization_theory = 1
	max_potency_theory[kappa] = K - normalization_theory

ax_K_mf.plot(kappas_theory, np.log10(np.array(list(max_potency_theory.values()))/(C/Kd_r_renorm))-0.4, color = my_purple, linestyle = '-', marker = '', linewidth = 3, label = 'theory largest', alpha = .8)

print(Kds[betas[:-1]>1][-1])
t_growth = (1/lambda_B)*np.log(C/50)

print('--------')
my_plot_layout(ax = ax_K_mf, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K_mf.legend(fontsize = 26, title_fontsize = 30, loc = 4)
ax_K_mf.set_xlim(left = 0.8, right = 4.2)
#ax_K_mf.set_ylim(bottom = -2.5, top = 0.15)
#ax_K_mf.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K_mf.set_yticklabels([1, 0.1, 0.01])
fig_K_mf.savefig('../../Figures/1_Dynamics/Ensemble/K_max_'+energy_model+'.pdf')



