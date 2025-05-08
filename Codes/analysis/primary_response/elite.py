import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------

L0s = [1e9, 1e7]#, 1e8]


T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1
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

E_lims = [-12]

antigen_color = my_yellow

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_red, my_green, my_blue2, my_gold, my_purple])
color_list = np.array([my_red, my_green, my_purple])

colors_L0 = []
for i in range(len(color_list)):
        colors_L0.append(np.array(color_list[i]))

# #colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
# colors_R = []
# for i in range(len(ps)):
#     colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18'

L=len(antigen)
print('--------')
print('l=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path + "primary_response/in/")
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))

#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

min_E = np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
avg_E = (np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]))
print(min_E, avg_E, E_m)
var_E = (np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]))
#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.001, 80, 1000000, PWM_data, E_m, L)
Kds = np.exp(Es[:-1])


#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_K_distribution2_total, ax_K_distribution2_total = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

Kd_r_renorm = Kds[(P_min_e_Q0(1e8, Q0, dE))==np.max(P_min_e_Q0(1e8, Q0, dE))][0]

for i_L0, L0 in enumerate(L0s):
	if L0 == 1e7:
		# print('-------- L0 = %.0e --------'%L0)
		# N_ensss = [[1000+i for i in range(0, 1, 1)] + [10000+i for i in range(0, 1, 1)] + [20000+i for i in range(0, 3, 1)] + [40000+i for i in range(0, 1, 1)]] # L0 1e7
		# L0 = 1e7
		# t_lims = [9.]
		# ps = [5.0]

		print('-------- L0 = %.0e --------'%L0)
		N_ensss = [[40000+i for i in range(0, 1, 1)] + [80000+i for i in range(0, 1, 1)] + [100000+i for i in range(0, 1, 1)]] # L0 1e7
		L0 = 1e7
		t_lims = [9.]
		ps = [4.0]

	if L0 == 1e8:
		print('--------L0 = %.0e --------'%L0)
		N_ensss = [[1000+i for i in range(0, 1, 1)] + [20000+i for i in range(0, 1, 1)] + [8000+i for i in range(0, 1, 1)] + [4000+i for i in range(0, 1, 1)]] # L0 1e7=8
		L0 = 1e8
		t_lims = [8.]
		ps = [5.0]

	if L0 == 1e9:
		print('--------L0 = %.0e --------'%L0)
		N_ensss = [[1000+i for i in range(0, 6, 1)] + [2000+i for i in range(0, 1, 1)] + [4000+i for i in range(0, 1, 1)] + [8000+i for i in range(0, 1, 1)]] # L0 1e9
		L0 = 1e9
		t_lims = [7.]
		ps = [4.0]

		# N_ensss = [[1000+i for i in range(0, 1, 1)] ] # L0 1e7
		# L0 = 1e9
		# t_lims = [8.]
		# ps = [5.0]
	#--------------------------Repertoire properties--------------------------
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
	print('beta_r = %.1f'%beta_r)
	fig_K_distribution, ax_K_distribution = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	fig_K_distribution2, ax_K_distribution2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

	for i_p, p in enumerate(ps):
		E_lim = E_lims[i_p]
		t_lim = t_lims[i_p]
		m_bar_theory = np.array([np.sum(L0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
		t_act_theory = time_array[m_bar_theory>1][0] 
		print('--------')
		print('p = %.2f...'%p)
		beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
		beta_act = np.min([beta_r, beta_p])

		Kd_r_renorm = Kds[(P_min_e_Q0(L0, Q0, dE))==np.max(P_min_e_Q0(L0, Q0, dE))][0]
		print('L0:%.0e and Z0:%.2e'%(L0, C/Kd_r_renorm))
		print('L0:%.0e and Z0:%.2e'%(L0, C/Kd_r))

		N_enss = N_ensss[i_p]

		#K = np.zeros_like(time_array)
		Z_final = []
		Counter = 0
		Z_final_largest = []
		Z_final_largest_2 = []
		#K_final_best_renorm = []
		counter_total = 0
		for N_ens in N_enss:
			print('N_ens = %d'%N_ens)

			#-----------------Loading data----------------------------
			parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
			#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
			return_data_type = 0
			data, return_data_type = get_data(folder_path = '../../out/primary_response', data_type = 'K_elite_Nens-%d_p-%.1f_L0-%.0e'%(N_ens, p, L0))
			
			if(return_data_type==-1):
				continue
			elif(return_data_type==1):
				Z_final_n = data[0]
				Counter_n = data[1]
				Z_final_largest_n = data[2]
				Z_final_largest_2_n = data[3]
			else:
				data = pd.read_csv(Text_files_path + 'primary_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L0, p, k_step, E_lim, t_lim, E_m))    
				Z_final_n = []
				Counter_n = 0
				Z_final_largest_n = []
				Z_final_largest_2_n = []

				for i_ens in tqdm(np.arange(0, N_ens, 1)):
					data_active = data.loc[data['ens_id']==i_ens]
					#data_active = data_i.loc[data_i['active']==1]
					t_act_data = np.min(data_active['time'])
					data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.2*(p-1))] # it was 1.0 + 0.1*...
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
						#Potencies = np.divide(((clone_sizes_C-1).T)/1, Kds_C).T
						#Z_i = np.sum(Potencies, axis = 0)
						Potencies = (clone_sizes_C[:, -1]-1)/Kds_C
						Z_i = np.sum(Potencies)	
						#if(Z_i[-1]>1e12):
						# if(Z_i>1e12):
						# 	print('----------- FOUND ELITE:', i_ens, 'at', N_ens, '-----------')
						# 	print(Z_i)
						#K_final_best_renorm.append(((C/1)/np.min(Kds_C)))
						#if(np.sum(Z_i)!=0):
						#Z_final_n.append(Z_i[-1])
						#Z_final_largest_n.append(np.max(Potencies[:, -1]))
						Z_final_n.append(Z_i)
						Z_final_largest_n.append(np.max(Potencies))
						Z_final_largest_2_n.append(C/np.min(Kds_C))
						Counter_n+=1

				f = open('../../out/primary_response/processed_data_K_elite_Nens-%d_p-%.1f_L0-%.0e.pkl'%(N_ens, p, L0), 'wb')
				pickle.dump([Z_final_n, Counter_n, Z_final_largest_n, Z_final_largest_2_n], f, pickle.HIGHEST_PROTOCOL)	

			Z_final+=Z_final_n
			Z_final_largest+=Z_final_largest_n
			Z_final_largest_2+=Z_final_largest_2_n
			Counter+=Counter_n
			print('counts in n = ',Counter_n)

		print('total counts = ', Counter)
		#Counter = 1
		normalization = 1
		print('%.2e'%np.max(Z_final))
		
		bins = np.logspace(np.min(np.log10(Z_final)), np.max(np.log10(Z_final)), 50)
		bins2 = np.logspace(np.min(np.log10(Z_final_largest)), np.max(np.log10(Z_final_largest)), 50)
		bins3 = np.logspace(np.min(np.log10(Z_final_largest_2)), np.max(np.log10(Z_final_largest_2)), 50)
		#bins = 'auto'
		#bins2 = 'auto'
		#bins = 500
		#bins2 = 500

		P_data_sim = plt.hist((np.array(Z_final)), bins = 'auto', density = False, cumulative  = True, alpha = 0)
		#P_data_sim = np.histogram(np.log10(np.array(Z_final)), bins = bins, density = True)

		P_data_sim2 = plt.hist((np.array(Z_final_largest)), bins = bins2, density = False, cumulative  = True, alpha = 0)
		#P_data_sim2 = np.histogram(np.log10(np.array(Z_final_largest)), bins = bins, density = True)

		P_data_sim3 = plt.hist((np.array(Z_final_largest_2)), bins = bins3, density = False, cumulative  = True, alpha = 0)
		#P_data_sim3 = np.histogram(np.log10(np.array(Z_final_largest_2)), bins = bins, density = True)

		Zs_sim = P_data_sim[1]
		P_Z_sim = P_data_sim[0]

		Zs_sim2 = P_data_sim2[1]
		P_Z_sim2 = P_data_sim2[0]

		Zs_sim3 = P_data_sim3[1]
		P_Z_sim3 = P_data_sim3[0]

		#P_data_sim = np.histogram(np.array(Z_final), bins = np.logspace(np.log10(np.min(Z_final)), np.log10(np.max(Z_final)), 200), density = False)
		# density = True
		#ax_K_distribution.plot((Zs_sim[:-1])/(C/Kd_r_renorm), P_Z_sim, color = 'darkorange', linestyle='', marker = 'D', linewidth = 2, ms = 4, label = r'$p=%.1f$'%(p))
		#ax_K_distribution.plot((Zs_sim2[:-1])/(C/Kd_r_renorm), P_Z_sim2, color = 'tab:cyan', linestyle='', marker = 'o', linewidth = 2, ms = 4, alpha = .6)
		# density = False
		#ax_K_distribution.plot((Zs_sim[:-1])/(C/Kd_r_renorm), P_Z_sim/Counter, color = 'darkorange', linestyle='', marker = 'D', linewidth = 2, ms = 5, label = r'$p=%.1f$'%(p))
		#ax_K_distribution.plot((Zs_sim2[:-1])/(C/Kd_r_renorm), P_Z_sim2/Counter, color = 'tab:cyan', linestyle='', marker = 'o', linewidth = 2, ms = 5)
		# cumulative  = True
		
		#ax_K_distribution.plot((Zs_sim2[:-2])/(C/Kd_r_renorm), np.diff(P_Z_sim2/Counter)/np.diff(Zs_sim2[:-1])*1e9, color = 'tab:cyan', linestyle='', marker = '*', linewidth = 2, ms = 8)
		ax_K_distribution.plot((Zs_sim3[:-2])/(C/Kd_r_renorm), np.diff(P_Z_sim3/Counter)/np.diff(Zs_sim3[:-1])*1e9, color = colors_L0[i_L0], linestyle='', marker = '*', linewidth = 2, ms = 8)
		ax_K_distribution.plot((Zs_sim[:-2])/(C/Kd_r_renorm), np.diff(P_Z_sim/Counter)/np.diff(Zs_sim[:-1])*1e9, color = colors_L0[i_L0], linestyle='', marker = '^', linewidth = 2, ms = 6, label = r'$p=%d$'%(p))

		# density = True	
		#ax_K_distribution2.plot(np.exp(Zs_sim[:-1])/(C/Kd_r_renorm), 1-np.cumsum(P_Z_sim*np.diff(Zs_sim)), color = 'darkorange', linestyle='-', marker = '', linewidth = 5, ms = 8, label = r'$p=%.1f$'%(p), zorder = 20)
		#ax_K_distribution2.plot(np.exp(Zs_sim2[:-1])/(C/Kd_r_renorm), 1-np.cumsum(P_Z_sim2*np.diff(Zs_sim2)), color = 'tab:cyan', linestyle='-', marker = '', linewidth = 5, ms = 8, zorder = 20, alpha = .6)
		# density = False
		#ax_K_distribution2.plot((np.flip(Zs_sim[:-1]))/(C/Kd_r_renorm), np.cumsum(np.flip(P_Z_sim))/Counter, color = 'darkorange', linestyle='-', marker = '', linewidth = 5, ms = 8, label = r'$p=%.1f$'%(p), zorder = 20)
		#ax_K_distribution2.plot((np.flip(Zs_sim2[:-1]))/(C/Kd_r_renorm), np.cumsum(np.flip(P_Z_sim2))/Counter, color = 'tab:cyan', linestyle='-', marker = '', linewidth = 5, ms = 8, zorder = 20)
		# cumulative  = True
		#ax_K_distribution2.plot(((Zs_sim2[:-1]))/(C/Kd_r_renorm), 1-P_Z_sim2/Counter, color = 'tab:cyan', linestyle='', marker = '*', linewidth = 5, ms = 8, zorder = 20)
		ax_K_distribution2.plot(((Zs_sim3[:-1]))/(C/Kd_r_renorm), 1-P_Z_sim3/Counter, color = colors_L0[i_L0], linestyle='-', marker = '', linewidth = 3, ms = 8, zorder = 20)
		ax_K_distribution2.plot(((Zs_sim[:-1]))/(C/Kd_r_renorm), 1-P_Z_sim/Counter, color = colors_L0[i_L0], linestyle='-', marker = '', linewidth = 5, ms = 6, label = r'$p=%d'%(p), zorder = 20)

		#ax_K_distribution2_total.plot(((Zs_sim3[:-1]))/(C/Kd_r_renorm), 1-P_Z_sim3/Counter, color = 'darkred', linestyle='', marker = '*', linewidth = 3, ms = 8, zorder = 20)
		if L0 == 1e7:
			ax_K_distribution2_total.plot(((Zs_sim[:-1]))/(C/Kd_r_renorm), 1-P_Z_sim/Counter, color = colors_L0[i_L0], linestyle='-', marker = '', linewidth = 3, ms = 6, label = r'$p=%d'%(p), zorder = 20)


		#print(np.sum(P_Z_sim*np.diff(Zs_sim)))
		#print(np.sum(P_Z_sim2*np.diff(Zs_sim2)))

		print(np.sum(P_Z_sim), np.sum(P_Z_sim/Counter), np.sum(P_Z_sim[-1]/Counter))
		print(np.sum(P_Z_sim2), np.sum(P_Z_sim2/Counter), np.sum(P_Z_sim2[-1]/Counter))

		mean_Q = np.sum(10**Zs_sim[:-1]*P_Z_sim*np.diff(Zs_sim)/Counter)

		#ax_K_distribution.vlines(mean_Q/(C/Kd_r_renorm), 1e-6, 1, color = my_red, linestyle='-', linewidth = 4)	
		#ax_K_distribution2.vlines(mean_Q/(C/Kd_r_renorm), 1e-6, 1, color = my_red, linestyle='-', linewidth = 4)

		#ax_K_distribution.vlines(2.3e-1, 1e-6, 1, color = my_red, linestyle='--', linewidth = 4)	
		ax_K_distribution.vlines(np.mean(Z_final)/(C/Kd_r_renorm), 1e-6, 1, color = colors_L0[i_L0], linestyle='-', linewidth = 4)	
		#ax_K_distribution2.vlines(2.3e-1, 1e-6, 1, color = my_red, linestyle='--', linewidth = 4)		
		ax_K_distribution2.vlines(np.mean(Z_final)/(C/Kd_r_renorm), 1e-6, 1, color = colors_L0[i_L0], linestyle='-', linewidth = 2)		
		if L0 != 1e8:
			ax_K_distribution2_total.vlines(np.mean(Z_final)/(C/Kd_r_renorm), 1e-6, 1, color = colors_L0[i_L0], linestyle='-', linewidth = 3, alpha = .6)		
		print('mean_Z = %.3e (%.3e)'%(np.mean(Z_final), np.mean(Z_final)/(C/Kd_r)))

		# Printing K from Gumbel
		Nb = C
		#Z_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
		Z_array = np.flip((Nb/1)/Kds)
		p_E = P_min_e_Q0(L0, Q0, dE)
		#p_E = P_min_e(L0, avg_E, var_E, Es[:-1], dE)
		print('N=', sum(p_E*dE))
		p_E/=np.sum(p_E*dE)
		p_Z = np.flip(p_E)/Z_array
		#p_Z /= np.sum(p_Z[:-1]*np.diff(Z_array))
		#p_Z = p_Z/np.sum(np.flip(p_Z[:-1])/np.flip(Z_array[:-1])*abs(np.diff(np.flip(Z_array))))
		print('N = ',np.sum((p_Z[:-1])*(np.diff((Z_array)))))

		p_Z /= np.sum(p_Z[:-1]*np.diff(Z_array))
		ax_K_distribution.plot((((Z_array[:-1]))/(C/Kd_r_renorm)), (p_Z[:-1]), linestyle = '--', marker = '',  color = colors_L0[i_L0], ms = 2, linewidth = 3, alpha = .5, label = 'Gumbel')
		ax_K_distribution2.plot((((Z_array[:-1]))/(C/Kd_r_renorm)), 1-np.cumsum((abs(p_Z[:-1]))*(np.diff((Z_array)))), linestyle = '--', marker = '',  color = colors_L0[i_L0], ms = 2, linewidth = 3, alpha = .5)
		
		if L0 != 1e8:
			ax_K_distribution2_total.plot((((Z_array[:-1]))/(C/Kd_r_renorm)), 1-np.cumsum((abs(p_Z[:-1]))*(np.diff((Z_array)))), linestyle = '--', marker = '',  color = colors_L0[i_L0], ms = 2, linewidth = 2, alpha = .7)
		#ax_K_distribution2.plot((((Z_array[:]))/(C/Kd_r_renorm)), 1-np.cumsum(np.flip(p_E[:]*dE)), linestyle = '-', marker = '',  color = 'darkred', ms = 2, linewidth = 3, alpha = .5)

		Z_array_tail = 10**np.linspace(9, 14.5, 100)

		alpha = (lambda_A*beta_r)/(lambda_B*p + lambda_A)
		zeta = (lambda_B*p)/(lambda_A*beta_r)

		print('exponent = ', -1/zeta)

		fit2p = (Z_array_tail**(-1/zeta-1))/((1e11)**(-1/zeta-1))
		fit2p *= (1e9*np.diff(P_Z_sim/Counter)/np.diff(Zs_sim[:-1]))[((Zs_sim[:-2]))<1e11][-1]

		fit2p2 = (Z_array_tail**(-1*alpha-1))/((1e11)**(-1*alpha-1))
		fit2p2 *= (1e9*np.diff(P_Z_sim2/Counter)/np.diff(Zs_sim2[:-1]))[((Zs_sim2[:-2]))<1e11][-1]

		ax_K_distribution.plot(Z_array_tail/(C/Kd_r_renorm), fit2p, linewidth = 3, color = colors_L0[i_L0], linestyle = ':', alpha = .8)
		#ax_K_distribution.plot(Z_array_tail/(C/Kd_r_renorm), fit2p2, linewidth = 3, color = 'tab:cyan', linestyle = '-')

		fit = Z_array_tail**(-(beta_r+1))/((Zs_sim[-2])**(-(beta_r+1)))
		#fit2 = Z_array_tail**(-2)/((1e11)**(-2))
		fit *= (1-P_Z_sim/Counter)[((Zs_sim[:-1]))<Zs_sim[-2]][-1]

		fit2 = Z_array_tail**(-(1/(zeta)))/((1e11)**(-(1/(zeta))))
		#fit2 = Z_array_tail**(-2)/((1e11)**(-2))
		fit2 *= (1-P_Z_sim/Counter)[((Zs_sim[:-1]))<1e11][-1]

		exp_potency = (lambda_A*beta_r*1.5/p)/(lambda_B + lambda_A/p)
		fit22 = (Z_array_tail**((-exp_potency)))/((1e11)**((-exp_potency)))
		#fit22 *= (1-np.cumsum(P_Z_sim*np.diff(Zs_sim2)))[((Zs_sim2[:-1]))<10**10.8][-1]
		fit22 *= (1-P_Z_sim2/Counter)[((Zs_sim2[:-1]))<1e11][-1]

		print('EXPONENTS:', 1/zeta, exp_potency)

		ax_K_distribution2.plot(Z_array_tail/(C/Kd_r_renorm), fit2, linewidth = 3, color = colors_L0[i_L0], linestyle = ':', alpha = .6)
		ax_K_distribution2.plot(Z_array_tail/(C/Kd_r_renorm), fit22, linewidth = 3, color = 'k', linestyle = ':', alpha = .6)
		
		if L0 != 1e8:
			#ax_K_distribution2_total.plot(Z_array_tail/(C/Kd_r_renorm), fit, linewidth = 4, color = colors_L0[i_L0], linestyle = '--', alpha = .7)
			ax_K_distribution2_total.plot(Z_array_tail/(C/Kd_r_renorm), fit2, linewidth = 3, color = colors_L0[i_L0], linestyle = ':', alpha = .7)
			# ax_K_distribution2_total.plot(Z_array_tail/(C/Kd_r_renorm), fit22, linewidth = 4, color = colors_L0[i_L0], linestyle = ':', alpha = .7)
		# ax_K_distribution2.plot(Z_array_tail/(C/Kd_r_renorm), fit22, linewidth = 3, color = 'tab:cyan', linestyle = '-')


		my_plot_layout(ax = ax_K_distribution, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		#ax_K_distribution.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
		ax_K_distribution.set_ylim(bottom = 1e-8, top = 2)
		ax_K_distribution.set_xlim(left = 1e-2, right = 10**(13.5)/(C/Kd_r_renorm))
		#ax_K_distribution.set_xticks([])
		#ax_K_distribution.set_yticks([])
		#ax_K_distribution.set_yticklabels([1, 0.1, 0.01])
		fig_K_distribution.savefig('../../../Figures/primary_response/1_Dynamics/CSV/Z_elite_P_'+energy_model+'_p-%.1f_L0-%.0e.pdf'%(p, L0))

		my_plot_layout(ax = ax_K_distribution2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		#ax_K_distribution2.legend(fontsize = 28, title_fontsize = 30, loc = 0)
		ax_K_distribution2.set_ylim(bottom = 8e-6, top = 1.05)
		ax_K_distribution2.set_xlim(left = (np.mean(Z_final)/(C/Kd_r_renorm))/3, right = 7e1)
		#ax_K_distribution2.set_xticks([])
		#ax_K_distribution2.set_yticks([])
		#ax_K_distribution2.set_yticklabels([1, 0.1, 0.01])
		fig_K_distribution2.savefig('../../../Figures/primary_response/1_Dynamics/CSV/Z_elite_F_'+energy_model+'_p-%.1f_L0-%.0e.pdf'%(p, L0))

		print((C/Kd_r_renorm))

my_plot_layout(ax = ax_K_distribution2_total, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K_distribution2_total.legend(fontsize = 28, title_fontsize = 30, loc = 0)
ax_K_distribution2_total.set_yticks([1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
ax_K_distribution2_total.set_ylim(bottom = 2e-6, top = 1.2)
ax_K_distribution2_total.set_xlim(left = 1e-1, right = 1e2)
#ax_K_distribution2_total.set_xticks([])
#ax_K_distribution2_total.set_yticklabels([1, 0.1, 0.01])
fig_K_distribution2_total.savefig('../../../Figures/primary_response/1_Dynamics/CSV/Z_elite_F_'+energy_model+'.pdf')
