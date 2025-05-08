import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------
N_ens = 400
N_enss = [N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens]

L_0 = 1e9
transparencies_p = [.8, 1, .8, .8, .8]

T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
p = 4
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*10
#N_c = 1e5
#E_m = -27.63
E_m = -24
C = 1e4
AA = 1

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

k_steps = np.array([1, 3, 6, 12, 30, 60, 120, 300, 600, 1200])#]) # hour^-1
E_lim = -12
t_lims = [9.0, 8.5, 8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5]#]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])
color_list = np.array([my_blue2, my_blue2, my_green, my_green, my_red, my_red, my_gold, my_red, my_red, my_gold, my_gold])


colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(k_steps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18'

L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path + "primary_immune_response/in/")
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_m, L)
Kds = np.exp(Es[:-1])

#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
print('L_0 = %.e'%L_0)
print('beta_r = %.1f'%beta_r)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Look_steps...')

print('p = %.2f...'%p)
beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
#--------------------------Look_steps--------------------------
#fig_mf, ax_mf = plt.subplots(figsize=(11,8), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_mf, ax_mf = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
#ax_mf_2 = ax_mf.twinx()
print('--------')
print('--- Processing all clones ---')
print('--------')
max_potency_simulations = dict()
max_potency_simulations_std = dict()
max_potency_theory = dict()

# Printing K from Gumbel
Nb = C
#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
K_array = ((Nb/1)/Kds)
p_K = P_min_e_Q0(L_0, Q0, dE)#/K_array**2*(Nb/1)
p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
Kd_r_renorm = Kds[(P_min_e_Q0(L_0, Q0, dE))==np.max(P_min_e_Q0(L_0, Q0, dE))]
ax_mf.hlines(0, .9, 4.1, linewidth = 2, color = 'indigo', linestyle = 'dotted')


for i_k_step, k_step in enumerate(k_steps):
	k_step = k_step*24 #days^-1
	t_lim = t_lims[i_k_step]
	# QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t_act_1))/N_A, Es, p, lambda_A, N_c, dE)[3]
	# Kd_1_renorm = Kds[(QR/1)==np.max(QR/1)]
	# #ax_mf.hlines(np.log10((C*1)/Kd_1_renorm/(C/Kd_r_renorm)), .9, 4.1, linewidth = 2, color = 'black', linestyle = 'dotted')
	# Kd_1_renorm = Kds[(QR/Kds)==np.max(QR/Kds)]
	#ax_mf.hlines(np.log10((C*1)/Kd_1_renorm/(C/Kd_r_renorm)), .9, 4.1, linewidth = 2, color = 'black', linestyle = 'dashed')
	#--------------------------Proofreading properties--------------------------
	beta_step, E_step, Kd_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
	print('k_step = %.2f'%(k_step/24), 'beta_step = %.2f'%beta_step)
	N_ens = N_enss[i_k_step]
	m_bar_theory = np.array([np.sum(L_0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
	t_act_theory = time_array[m_bar_theory>1][0] 

	if p==1:
		t_act_1 = t_act_theory
	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3, k_step/24, p, N_c, linear, N_ens)+energy_model
	return_data_type = 0
	data, return_data_type = get_data_b(folder_path = '../../out/primary_immune_response', data_type = 'K_mf_2_p-%.1f_k_step-%.d'%(p, k_step))

	if(return_data_type==1):
		K_final = data[0]
		Counter = data[1]
	elif(return_data_type==0):
		data = pd.read_csv(Text_files_path + 'primary_immune_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
		K_final = []
		Counter = 0
		for i_ens in tqdm(np.arange(N_ens)):
			data_active = data.loc[data['ens_id']==i_ens]
			#data_active = data_i.loc[data_i['active']==1]
			t_act_data = np.min(data_active['time'])
			data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
			activation_times = np.array(data_active['time'])
			energies  = np.array(data_active['E'])

			#---------------------------- B cell linages ----------------------
			clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
			#--------------------------t_C filter-------------------------
			lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
			#lim_size = 2
			clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
			#-------Simulations-------
			if(len(energies_C)>0):
				Kds_C = np.exp(energies_C)

				K_i = [np.sum(((clone_sizes_C[:,t]-1)/1)/Kds_C) for t in np.arange(len(time_array))]

				if(np.sum(K_i)!=0):
					K_final.append(K_i[-1])
					Counter+=1

				#if(i_ens%1==0):
				#	ax_K.plot(time_array, K_i, color = colors_p[i_k_step], alpha = .1, linewidth = 1)

		f = open('../../out/primary_immune_response/processed_data_K_mf_2_p-%.1f_k_step-%.d.pkl'%(p, k_step), 'wb')
		pickle.dump([K_final, Counter], f, pickle.HIGHEST_PROTOCOL)	
		

	max_potency_simulations[k_step] = np.mean(np.array(K_final)/(C/Kd_r_renorm))
	max_potency_simulations_std[k_step] = np.sqrt(np.var((np.array(K_final)/(C/Kd_r_renorm))))
	print('--------')

ax_mf.plot(k_steps/(60), (np.array(list(max_potency_simulations.values()))), color = 'indigo', linestyle = '-', marker = '^', linewidth = 3, label = 'Total', ms = 10, alpha = 1)
#ax_mf.errorbar(x=k_steps, y=(np.array(list(max_potency_simulations.values()))), yerr = (np.array(list(max_potency_simulations_std.values()))), ls = 'none', color = 'indigo', linewidth = 2, alpha = .8)

#-------------------------# 
print('--------')
print('--- Processing largest clone ---')
print('--------')
max_potency_simulations2 = dict()
max_potency_simulations_std2 = dict()
max_potency_theory2 = dict()

max_potency_simulations3 = dict()
max_potency_simulations_std3 = dict()
max_potency_theory3 = dict()

for i_k_step, k_step in enumerate(k_steps):
	k_step = k_step*24 #days^-1
	t_lim = t_lims[i_k_step]
	N_ens = N_enss[i_k_step]
	m_bar_theory = np.array([np.sum(L_0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
	t_act_theory = time_array[m_bar_theory>1][0] 
	print('--------')
	#--------------------------Proofreading properties--------------------------
	beta_step, E_step, Kd_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
	print('k_step = %.2f'%(k_step/24), 'beta_step = %.2f'%beta_step)

	beta_act = np.min([beta_r, beta_p])

	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
	return_data_type = 0
	data, return_data_type = get_data_b(folder_path = '../../out/primary_immune_response', data_type = 'K_largest_2_k_step-%d'%k_step)
		

	if(return_data_type):
		final_biggest = data[0]
		final_biggest_affinity = data[1]
		final_highest_potency = data[2]
	else:	
		data = pd.read_csv(Text_files_path + 'primary_immune_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
		final_biggest = []
		final_biggest_affinity = []
		final_highest_potency = []

		for i_ens in tqdm(np.arange(N_ens)):
			data_active = data.loc[data['ens_id']==i_ens]
			#data_active = data_i.loc[data_i['active']==1]
			t_act_data = np.min(data_active['time'])
			data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
			activation_times = np.array(data_active['time'])
			energies  = np.array(data_active['E'])

			#---------------------------- B cell linages ----------------------
			clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
			#--------------------------t_C filter-------------------------
			lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
			clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
			#-------Simulations-------
			potencies_C = (clone_sizes_C.T/np.exp(energies_C)).T

			if(len(energies_C)>0):
				
				final_biggest_affinity.append(energies_C[clone_sizes_C[:,-1]==np.max(clone_sizes_C[:,-1])][0])
				final_biggest.append(np.max(clone_sizes_C[:,-1]))
				final_highest_potency.append(np.max(potencies_C[:,-1]))

		f = open('../../out/primary_immune_response/processed_data_K_largest_2_k_step-%d.pkl'%k_step, 'wb')
		pickle.dump([final_biggest, final_biggest_affinity, final_highest_potency], f, pickle.HIGHEST_PROTOCOL)	
		
	max_potency_simulations2[k_step] = np.mean((np.array(final_biggest)/np.exp(final_biggest_affinity))/(C/Kd_r_renorm))
	max_potency_simulations_std2[k_step] = np.sqrt(np.var(((np.array(final_biggest)/np.exp(final_biggest_affinity))/(C/Kd_r_renorm))))

	max_potency_simulations3[k_step] = np.mean((np.array(final_highest_potency))/(C/Kd_r_renorm))
	max_potency_simulations_std3[k_step] = np.sqrt(np.var((np.array(final_highest_potency))/(C/Kd_r_renorm)))

ax_mf.plot(k_steps/(60), (np.array(list(max_potency_simulations2.values()))), color = 'indigo', linestyle = '--', marker = '*', linewidth = 3, label = 'Largest', ms = 14, alpha = .8)
#ax_mf.plot(k_steps, np.log10(np.array(list(max_potency_simulations3.values()))), color = 'indigo', linestyle = '', marker = '*', linewidth = 3, label = 'Largest', ms = 14, alpha = .6)
#ax_mf.errorbar(x=k_steps, y=(np.array(list(max_potency_simulations2.values()))), yerr = (np.array(list(max_potency_simulations_std2.values()))), ls = 'none', color = 'indigo', alpha = .6)

print('--------')

#ax_mf.spines['top'].set_color('indigo')
#ax_mf.xaxis.label.set_color('indigo')
#ax_mf.tick_params(axis='y', colors='indigo')

my_plot_layout(ax = ax_mf, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_mf.legend(fontsize = 26, title_fontsize = 30, loc = 8)
ax_mf.set_xlim(left = .5/(60), right = 1800/(60) )
#ax_mf.set_ylim(bottom = 1e-3, top = 1.1)
#ax_mf.set_yticks([1, 0.1, 0.01, 0.001])
#ax_mf.set_yticklabels([1, 0.1, 0.01])

print('--------')
print('--- Entropy ---')
print('--------')
# Delta_S = dict()
# N_acts = dict()
# fig_entropy, ax_entropy = plt.subplots(figsize=(4.5*2,3.0*1.805), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
# fig_N_act, ax_N_act = plt.subplots(figsize=(4.5*2,3.0*1.805), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
# for freq in np.array([0.0100])*10:
# 	# ENTROPYYYYYYYY
# 	N_acts[freq] = dict()
# 	max_entropy_simulations = dict()
# 	max_entropy_simulations_std = dict()
# 	max_entropy_theory = dict()
# 	S1 = 1
# 	S4 = 1
# 	for i_k_step, k_step in enumerate(k_steps):
# 		k_step = k_step*24 #days^-1
# 		N_ens = N_enss[i_k_step]
# 		#--------------------------Proofreading properties--------------------------
# 		beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
# 		print('beta_pr = %.2f'%beta_pr)
# 		#-----------------Loading data----------------------------
# 		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
# 		data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'S_mf_2_%.3f'%freq)

# 		if(return_data_type):
# 			S_final = data[0]
# 			N_final = data[1]
# 			N_act = data[2]
# 			Counter = data[3]
# 		else:
# 			S_final = []
# 			N_final = []
# 			N_act = []
# 			Counter = 0
# 			for i_ens in tqdm(np.arange(N_ens)):
# 				data_active = data.loc[data['ens_id']==i_ens]
# 				#data_active = data_i.loc[data_i['active']==1]
# 				t_act_data = np.min(data_active['time'])
# 				data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
# 				activation_times = np.array(data_active['time'])
# 				energies  = np.array(data_active['E'])
# 				#---------------------------- B cell linages ----------------------
# 				clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
# 				#--------------------------t_C filter-------------------------
# 				#lim_size = np.max([int(C*freq), 2])
# 				lim_size = np.max([int(np.max(clone_sizes[:, -1])*freq), 2])
# 				clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
# 				N_act.append(len(activation_times_C))
# 				#-------Simulations-------
# 				#unactivated_clones = np.array([np.size(activation_times_C[activation_times_C>time_array[t]]) for t in range(len(time))])
# 				#total_size = np.sum(clone_sizes_C, axis = 0)
# 				#S_i = -np.array([np.sum(((clone_sizes_C[:, t])/total_size[t])*np.log((clone_sizes_C[:, t]-np.size(activation_times_C[activation_times_C>time_array[t]]))/total_size[t])) for t in range(len(time_array))])
# 				#S_i = np.nan_to_num(S_i, nan = 0, posinf = 0, neginf = 0)
# 				sort_inds = clone_sizes_C[:, -1].argsort()
# 				clone_sizes_C_sorted = clone_sizes_C[sort_inds, :][:,:]
# 				N_final_i = np.sum(clone_sizes_C_sorted[:, -1])
# 				S_final_i = -np.sum((clone_sizes_C_sorted[:, -1]/N_final_i)*np.log(clone_sizes_C_sorted[:, -1]/N_final_i))
# 				if(S_final_i!=0):
# 					S_final.append(S_final_i)
# 					N_final.append(N_final_i)#/np.log(N_final))
# 					Counter+=1
# 				#if(i_ens%1==0):
# 				#	ax_mf.plot(time_array, entropy_i, color = colors_p[  i_k_step], alpha = .1, linewidth = 1)

# 			f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_S_mf_2_%.3f.pkl'%freq, 'wb')
# 			pickle.dump([S_final, N_final, N_act, Counter], f, pickle.HIGHEST_PROTOCOL)	
		
# 		if(p==1):
# 			S1 = np.mean(S_final)

# 		if(p==4):
# 			S4 = np.mean(S_final)

# 		S_final = np.array(S_final)
# 		N_acts[freq][k_step] = N_act
# 		max_entropy_simulations[k_step] = np.mean(S_final)
# 		max_entropy_simulations_std[k_step] = np.std(S_final)

# 	Delta_S[freq] = S1 - S4
# 	k_steps_theory = np.linspace(.5, 4, 400)
# 	y_interp2 = np.interp(k_steps_theory, k_steps, np.array(list(max_entropy_simulations.values())))

# 	ax_entropy.plot(k_steps, np.array(list(max_entropy_simulations.values())), linestyle = '-', marker = '^', linewidth = 3, label = r'$%.3f$'%freq, ms = 10, alpha = 1)
# 	ax_N_act.plot(k_steps, [np.mean(N_acts[freq][k_step]) for k_step in k_steps*24], linestyle = '-', marker = '^', linewidth = 3, label = r'$%.3f$'%freq, ms = 10, alpha = 1)
# 	#ax_entropy.plot(k_steps_theory, y_interp2, linestyle = '-', marker = '', linewidth = 3, ms = 10, alpha = 1)

# 	if(freq==0.1):
# 		ax_mf_2.plot(k_steps, np.array(list(max_entropy_simulations.values())), color = 'olive', linestyle = '-', marker = '^', linewidth = 3, label = 'Total', ms = 10, alpha = 1)
# 		#ax_mf_2.errorbar(x=k_steps, y=(np.array(list(max_entropy_simulations.values()))), yerr = (np.array(list(max_entropy_simulations_std.values()))), ls = 'none', color = 'olive', alpha = .6)

# 		if np.any(y_interp2>4.01):
# 			#ax_mf_2.scatter(k_steps_theory[y_interp2>4.01][-1], 4.01, s = 180, facecolors='none', edgecolors='olive', marker = 'o', lw=2)
# 			#ax_mf_2.errorbar(x = k_steps_theory[y_interp2>4.01][-1], y = 4.01, yerr = 0.57, xerr = np.array([[k_steps_theory[y_interp2>4.01][-1] - k_steps_theory[y_interp2>(4.01+0.57)][-1], k_steps_theory[y_interp2>(4.01-0.57)][-1] - k_steps_theory[y_interp2>4.01][-1]]]).T, color = 'olive', alpha = .6, fmt='ko')
# 			print(k_steps_theory[y_interp2>4.01][-1], k_steps_theory[y_interp2>(4.01+0.57)][-1], k_steps_theory[y_interp2>(4.01-0.57)][-1])

# ax_mf_2.spines['top'].set_color('olive')
# ax_mf_2.xaxis.label.set_color('olive')
# ax_mf_2.tick_params(axis='y', colors='olive')

# my_plot_layout(ax = ax_mf_2, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
# #ax_mf_2.legend(fontsize = 28, title_fontsize = 30, loc = 8)
# #ax_mf_2.set_xlim(left = 0.8, right = 4.2)
# #ax_mf_2.set_ylim(bottom = 0)
# #ax_mf_2.set_yticks([1, 0.1, 0.01, 0.001])
# #ax_mf_2.set_yticklabels([1, 0.1, 0.01])

# my_plot_layout(ax = ax_entropy, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
# ax_entropy.legend(fontsize = 26, title_fontsize = 24, loc = 0, title = r'$f_{\rm{min}}$')
# #ax_entropy.set_xlim(left = 0.8, right = 4.2)
# #ax_entropy.set_ylim(bottom = 0)
# #ax_entropy.set_yticks([1, 0.1, 0.01, 0.001])
# #ax_entropy.set_yticklabels([1, 0.1, 0.01])
# fig_entropy.savefig('../../Figures/1_Dynamics/Ensemble/L%d/MF_2_entropy_'%(L)+energy_model+'.pdf')

# my_plot_layout(ax = ax_N_act, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
# ax_N_act.legend(fontsize = 26, title_fontsize = 24, loc = 0, title = r'$f_{\rm{min}}$')
# #ax_N_act.set_xlim(left = 0.8, right = 4.2)
# #ax_N_act.set_ylim(bottom = 0)
# #ax_N_act.set_yticks([1, 0.1, 0.01, 0.001])
# #ax_N_act.set_yticklabels([1, 0.1, 0.01])
# fig_N_act.savefig('../../Figures/1_Dynamics/Ensemble/L%d/MF_2_N_act_'%(L)+energy_model+'.pdf')


fig_mf.savefig('../../../Figures/primary_immune_response/1_Dynamics/CSV/MF_step_'+energy_model+'.pdf')

#print(Delta_S)
