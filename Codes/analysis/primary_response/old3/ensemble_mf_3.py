import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 500
N_enss = [N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens]

transparencies_p = [.8, 1, .8, .8, .8]

T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
k_pr = 120 # hour^-1
k_pr = k_pr*24 #days^-1
p = 2
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*1000
#N_c = 1e5
#E_ms = -27.63
E_ms = -25
C = 1e4
AA = 1

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

N_rs = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])
color_list = np.array([my_blue2, my_blue2, my_green, my_green, my_red, my_red, my_gold, my_red, my_red, my_gold, my_gold])


colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(N_rs)):
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
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path)
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Look_prs...')

print('p = %.2f...'%p)
beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
#--------------------------Look_prs--------------------------
#fig_mf, ax_mf = plt.subplots(figsize=(11,8), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_mf, ax_mf = plt.subplots(figsize=(4.5*2,3.0*1.805), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
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
p_K = P_min_e_Q0(1e8, Q0, dE)#/K_array**2*(Nb/1)
p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
Kd_r_renorm = Kds[(P_min_e_Q0(1e8, Q0, dE)/Kds)==np.max(P_min_e_Q0(1e8, Q0, dE)/Kds)]
ax_mf.hlines(0, .9, 4.1, linewidth = 2, color = 'indigo', linestyle = 'dotted')


for i_N_r, N_r in enumerate(N_rs):

	# QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_1))/N_A, Es, p, lambda_A, N_c, dE)[3]
	# Kd_1_renorm = Kds[(QR/1)==np.max(QR/1)]
	# #ax_mf.hlines(np.log10((C*1)/Kd_1_renorm/(C/Kd_r_renorm)), .9, 4.1, linewidth = 2, color = 'black', linestyle = 'dotted')
	# Kd_1_renorm = Kds[(QR/Kds)==np.max(QR/Kds)]
	#ax_mf.hlines(np.log10((C*1)/Kd_1_renorm/(C/Kd_r_renorm)), .9, 4.1, linewidth = 2, color = 'black', linestyle = 'dashed')
	#--------------------------Repertoire properties--------------------------
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
	print('N_r = %.e'%N_r)
	print('beta_r = %.1f'%beta_r)
	N_ens = N_enss[i_N_r]
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
	t_act_theory = time_array[m_bar_theory>1][0] 
	print('--------')

	if p==1:
		t_act_1 = t_act_theory
	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_Linear-%d_N_ens-%d_'%(lambda_A, 3, k_pr/24, p, N_c, linear, N_ens)+energy_model
	data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'K_mf_3')
	if(return_data_type==1):
		K_final = data[0]
		Counter = data[1]
	elif(return_data_type==0):
		K_final = []
		Counter = 0
		for i_ens in tqdm(np.arange(N_ens)):
			data_i = data.loc[data['i_ens']==i_ens]
			data_active = data_i.loc[data_i['active']==1]
			t_act_data = np.min(data_active['act_time'])
			data_active = data_active.loc[data_active['act_time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
			activation_times = np.array(data_active['act_time'])
			energies  = np.array(data_active['energy'])

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
				#	ax_K.plot(time_array, K_i, color = colors_p[i_N_r], alpha = .1, linewidth = 1)

		f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_K_mf_3.pkl', 'wb')
		pickle.dump([K_final, Counter], f, pickle.HIGHEST_PROTOCOL)	
		

	max_potency_simulations[N_r] = np.mean(np.array(K_final)/(C/Kd_r_renorm))
	max_potency_simulations_std[N_r] = np.sqrt(np.var((np.array(K_final)/(C/Kd_r_renorm))))

ax_mf.plot(N_rs, (np.array(list(max_potency_simulations.values()))), color = 'indigo', linestyle = '-', marker = '^', linewidth = 3, label = 'Total', ms = 10, alpha = 1)
#ax_mf.errorbar(x=k_prs, y=(np.array(list(max_potency_simulations.values()))), yerr = (np.array(list(max_potency_simulations_std.values()))), ls = 'none', color = 'indigo', linewidth = 2, alpha = .8)

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

for i_N_r, N_r in enumerate(N_rs):

	N_ens = N_enss[i_N_r]
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
	t_act_theory = time_array[m_bar_theory>1][0] 
	print('--------')
	#--------------------------Repertoire properties--------------------------
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
	print('N_r = %.e'%N_r)
	print('beta_r = %.1f'%beta_r)

	beta_act = np.min([beta_r, beta_p])

	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, p, N_c, linear, N_ens)+energy_model
	data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'K_largest_3')

	if(return_data_type):
		final_biggest = data[0]
		final_biggest_affinity = data[1]
		final_highest_potency = data[2]
	else:	
		final_biggest = []
		final_biggest_affinity = []
		final_highest_potency = []

		for i_ens in tqdm(np.arange(N_ens)):
			data_i = data.loc[data['i_ens']==i_ens]
			data_active = data_i.loc[data_i['active']==1]
			t_act_data = np.min(data_active['act_time'])
			data_active = data_active.loc[data_active['act_time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
			activation_times = np.array(data_active['act_time'])
			energies  = np.array(data_active['energy'])

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

		f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_K_largest_3.pkl', 'wb')
		pickle.dump([final_biggest, final_biggest_affinity, final_highest_potency], f, pickle.HIGHEST_PROTOCOL)	
		
	max_potency_simulations2[N_r] = np.mean((np.array(final_biggest)/np.exp(final_biggest_affinity))/(C/Kd_r_renorm))
	max_potency_simulations_std2[N_r] = np.sqrt(np.var(((np.array(final_biggest)/np.exp(final_biggest_affinity))/(C/Kd_r_renorm))))

	max_potency_simulations3[N_r] = np.mean((np.array(final_highest_potency))/(C/Kd_r_renorm))
	max_potency_simulations_std3[N_r] = np.sqrt(np.var((np.array(final_highest_potency))/(C/Kd_r_renorm)))

ax_mf.plot(N_rs, (np.array(list(max_potency_simulations2.values()))), color = 'indigo', linestyle = '--', marker = '*', linewidth = 3, label = 'Largest', ms = 14, alpha = .8)
#ax_mf.plot(k_prs, np.log10(np.array(list(max_potency_simulations3.values()))), color = 'indigo', linestyle = '', marker = '*', linewidth = 3, label = 'Largest', ms = 14, alpha = .6)
#ax_mf.errorbar(x=k_prs, y=(np.array(list(max_potency_simulations2.values()))), yerr = (np.array(list(max_potency_simulations_std2.values()))), ls = 'none', color = 'indigo', alpha = .6)

print('--------')

ax_mf.spines['top'].set_color('indigo')
ax_mf.xaxis.label.set_color('indigo')
ax_mf.tick_params(axis='y', colors='indigo')

my_plot_layout(ax = ax_mf, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_mf.legend(fontsize = 26, title_fontsize = 30, loc = 8)
ax_mf.set_xlim(left = 5e2, right = 5e9 )
#ax_mf.set_ylim(bottom = 1e-3, top = 1.1)
ax_mf.set_ylim(bottom = 7e-4, top = 1e-1)
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
# 	for i_N_r, k_pr in enumerate(k_prs):
# 		k_pr = k_pr*24 #days^-1
# 		N_ens = N_enss[i_N_r]
# 		#--------------------------Proofreading properties--------------------------
# 		beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
# 		print('beta_pr = %.2f'%beta_pr)
# 		#-----------------Loading data----------------------------
# 		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, p, N_c, linear, N_ens)+energy_model
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
# 				data_active = data.loc[data['i_ens']==i_ens]
# 				#data_active = data_i.loc[data_i['active']==1]
# 				t_act_data = np.min(data_active['act_time'])
# 				data_active = data_active.loc[data_active['act_time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
# 				activation_times = np.array(data_active['act_time'])
# 				energies  = np.array(data_active['energy'])
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
# 				#	ax_mf.plot(time_array, entropy_i, color = colors_p[  i_N_r], alpha = .1, linewidth = 1)

# 			f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_S_mf_2_%.3f.pkl'%freq, 'wb')
# 			pickle.dump([S_final, N_final, N_act, Counter], f, pickle.HIGHEST_PROTOCOL)	
		
# 		if(p==1):
# 			S1 = np.mean(S_final)

# 		if(p==4):
# 			S4 = np.mean(S_final)

# 		S_final = np.array(S_final)
# 		N_acts[freq][k_pr] = N_act
# 		max_entropy_simulations[k_pr] = np.mean(S_final)
# 		max_entropy_simulations_std[k_pr] = np.std(S_final)

# 	Delta_S[freq] = S1 - S4
# 	k_prs_theory = np.linspace(.5, 4, 400)
# 	y_interp2 = np.interp(k_prs_theory, k_prs, np.array(list(max_entropy_simulations.values())))

# 	ax_entropy.plot(k_prs, np.array(list(max_entropy_simulations.values())), linestyle = '-', marker = '^', linewidth = 3, label = r'$%.3f$'%freq, ms = 10, alpha = 1)
# 	ax_N_act.plot(k_prs, [np.mean(N_acts[freq][k_pr]) for k_pr in k_prs*24], linestyle = '-', marker = '^', linewidth = 3, label = r'$%.3f$'%freq, ms = 10, alpha = 1)
# 	#ax_entropy.plot(k_prs_theory, y_interp2, linestyle = '-', marker = '', linewidth = 3, ms = 10, alpha = 1)

# 	if(freq==0.1):
# 		ax_mf_2.plot(k_prs, np.array(list(max_entropy_simulations.values())), color = 'olive', linestyle = '-', marker = '^', linewidth = 3, label = 'Total', ms = 10, alpha = 1)
# 		#ax_mf_2.errorbar(x=k_prs, y=(np.array(list(max_entropy_simulations.values()))), yerr = (np.array(list(max_entropy_simulations_std.values()))), ls = 'none', color = 'olive', alpha = .6)

# 		if np.any(y_interp2>4.01):
# 			#ax_mf_2.scatter(k_prs_theory[y_interp2>4.01][-1], 4.01, s = 180, facecolors='none', edgecolors='olive', marker = 'o', lw=2)
# 			#ax_mf_2.errorbar(x = k_prs_theory[y_interp2>4.01][-1], y = 4.01, yerr = 0.57, xerr = np.array([[k_prs_theory[y_interp2>4.01][-1] - k_prs_theory[y_interp2>(4.01+0.57)][-1], k_prs_theory[y_interp2>(4.01-0.57)][-1] - k_prs_theory[y_interp2>4.01][-1]]]).T, color = 'olive', alpha = .6, fmt='ko')
# 			print(k_prs_theory[y_interp2>4.01][-1], k_prs_theory[y_interp2>(4.01+0.57)][-1], k_prs_theory[y_interp2>(4.01-0.57)][-1])

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


fig_mf.savefig('../../Figures/1_Dynamics/Ensemble/L%d/MF_3_'%(L)+energy_model+'.pdf')

#print(Delta_S)
