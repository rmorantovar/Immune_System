import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
#N_ensss = [[400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 501, 502, 503, 504, 505, 506, 507, 508, 509]] #for p=2.5
N_ensss = [[400+i for i in range(1, 51, 1)] + [600+i for i in range(1, 100, 1)] + [700+i for i in range(1, 100, 1)] + [800+i for i in range(1, 100, 1)] + [900+i for i in range(1, 100, 1)] + [1000+i for i in range(1, 100, 1)]] #for p=3
#N_ensss = [[400+i for i in range(1, 100, 1)]]
N_r = 1e8

T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
k_pr = k_pr*24 #days^-1
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

ps = [3.0]
#ps = [4.0]
#ps = [2.5]

antigen_color = my_yellow

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_red, my_green, my_blue2, my_gold, my_purple])
color_list = np.array([my_red])

colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(ps)):
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

min_E = np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
avg_E = (np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]))
print(min_E, avg_E, E_ms)
var_E = (np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]))
#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.001, 80, 1000000, PWM_data, E_ms, L)
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
#fig_K_distribution, ax_K_distribution = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
#fig_K_distribution2, ax_K_distribution2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

fig_K_distribution, ax_K_distribution = plt.subplots(figsize=(4.5*2,3.0*2), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_distribution2, ax_K_distribution2 = plt.subplots(figsize=(4.5*2,3.0*2), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

for i_p, p in enumerate(ps):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
	t_act_theory = time_array[m_bar_theory>1][0] 
	print('--------')
	print('p = %.2f...'%p)
	beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
	beta_act = np.min([beta_r, beta_p])

	Kd_r_renorm = Kds[(P_min_e_Q0(N_r, Q0, dE))==np.max(P_min_e_Q0(N_r, Q0, dE))]
	
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
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, p, N_c, linear, N_ens)+energy_model
		#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d_elite/'%L+parameters_path, data_type = 'K_elite')
		if(return_data_type==-1):
			continue
		elif(return_data_type==1):
			Z_final_n = data[0]
			Counter_n = data[1]
			Z_final_largest_n = data[2]
			Z_final_largest_2_n = data[3]
		else:
			Z_final_n = []
			Counter_n = 0
			Z_final_largest_n = []
			Z_final_largest_2_n = []

			for i_ens in tqdm(np.arange(N_ens/1)):
				data_active = data.loc[data['i_ens']==i_ens]
				#data_active = data_i.loc[data_i['active']==1]
				t_act_data = np.min(data_active['act_time'])
				data_active = data_active.loc[data_active['act_time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
				activation_times = np.array(data_active['act_time'])
				energies  = np.array(data_active['energy'])

				#---------------------------- B cell linages ----------------------
				clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
				#--------------------------t_C filter-------------------------
				#lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
				lim_size = 2
				clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
				#-------Simulations-------
				if(len(energies_C)>0):
					Kds_C = np.exp(energies_C)
					#Potencies = np.divide(((clone_sizes_C-1).T)/1, Kds_C).T
					#Z_i = np.sum(Potencies, axis = 0)
					Potencies = (clone_sizes_C[:, -1]-1)/Kds_C
					Z_i = np.sum(Potencies)	
					#if(Z_i[-1]>1e12):
					if(Z_i>1e12):
						print('----------- FOUND ELITE:', i_ens, 'at', N_ens, '-----------')
						print(Z_i)
					#K_final_best_renorm.append(((C/1)/np.min(Kds_C)))

					#if(np.sum(Z_i)!=0):
					
					#Z_final_n.append(Z_i[-1])
					#Z_final_largest_n.append(np.max(Potencies[:, -1]))
					Z_final_n.append(Z_i)
					Z_final_largest_n.append(np.max(Potencies))
					Z_final_largest_2_n.append(C/np.min(Kds_C))
					Counter_n+=1

			f = open(Text_files_path + 'Dynamics/Ensemble/L%d_elite/'%L+parameters_path+'/processed_data_K_elite.pkl', 'wb')
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
	
	bins = np.logspace(np.min(np.log10(Z_final)), np.max(np.log10(Z_final)), 100)
	bins2 = np.logspace(np.min(np.log10(Z_final_largest)), np.max(np.log10(Z_final_largest)), 100)
	bins3 = np.logspace(np.min(np.log10(Z_final_largest_2)), np.max(np.log10(Z_final_largest_2)), 100)
	#bins = 'auto'
	#bins2 = 'auto'
	#bins = 500
	#bins2 = 500

	P_data_sim = plt.hist((np.array(Z_final)), bins = bins, density = False, cumulative  = True, alpha = 0)
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
	ax_K_distribution.plot((Zs_sim3[:-2])/(C/Kd_r_renorm), np.diff(P_Z_sim3/Counter)/np.diff(Zs_sim3[:-1])*1e9, color = 'darkred', linestyle='', marker = '*', linewidth = 2, ms = 8)
	ax_K_distribution.plot((Zs_sim[:-2])/(C/Kd_r_renorm), np.diff(P_Z_sim/Counter)/np.diff(Zs_sim[:-1])*1e9, color = my_red, linestyle='', marker = '^', linewidth = 2, ms = 6, label = r'$p=%d$'%(p))

	# density = True	
	#ax_K_distribution2.plot(np.exp(Zs_sim[:-1])/(C/Kd_r_renorm), 1-np.cumsum(P_Z_sim*np.diff(Zs_sim)), color = 'darkorange', linestyle='-', marker = '', linewidth = 5, ms = 8, label = r'$p=%.1f$'%(p), zorder = 20)
	#ax_K_distribution2.plot(np.exp(Zs_sim2[:-1])/(C/Kd_r_renorm), 1-np.cumsum(P_Z_sim2*np.diff(Zs_sim2)), color = 'tab:cyan', linestyle='-', marker = '', linewidth = 5, ms = 8, zorder = 20, alpha = .6)
	# density = False
	#ax_K_distribution2.plot((np.flip(Zs_sim[:-1]))/(C/Kd_r_renorm), np.cumsum(np.flip(P_Z_sim))/Counter, color = 'darkorange', linestyle='-', marker = '', linewidth = 5, ms = 8, label = r'$p=%.1f$'%(p), zorder = 20)
	#ax_K_distribution2.plot((np.flip(Zs_sim2[:-1]))/(C/Kd_r_renorm), np.cumsum(np.flip(P_Z_sim2))/Counter, color = 'tab:cyan', linestyle='-', marker = '', linewidth = 5, ms = 8, zorder = 20)
	# cumulative  = True
	#ax_K_distribution2.plot(((Zs_sim2[:-1]))/(C/Kd_r_renorm), 1-P_Z_sim2/Counter, color = 'tab:cyan', linestyle='', marker = '*', linewidth = 5, ms = 8, zorder = 20)
	#ax_K_distribution2.plot(((Zs_sim3[:-1]))/(C/Kd_r_renorm), 1-P_Z_sim3/Counter, color = 'darkred', linestyle='', marker = '*', linewidth = 3, ms = 8, zorder = 20)
	ax_K_distribution2.plot(((Zs_sim[:-1]))/(C/Kd_r_renorm), 1-P_Z_sim/Counter, color = my_red, linestyle='', marker = '^', linewidth = 5, ms = 6, label = r'$p=%d'%(p), zorder = 20)

	#print(np.sum(P_Z_sim*np.diff(Zs_sim)))
	#print(np.sum(P_Z_sim2*np.diff(Zs_sim2)))

	print(np.sum(P_Z_sim/Counter))
	print(np.sum(P_Z_sim2/Counter))

	mean_Q = np.sum(10**Zs_sim[:-1]*P_Z_sim*np.diff(Zs_sim)/Counter)

	#ax_K_distribution.vlines(mean_Q/(C/Kd_r_renorm), 1e-6, 1, color = my_red, linestyle='-', linewidth = 4)	
	#ax_K_distribution2.vlines(mean_Q/(C/Kd_r_renorm), 1e-6, 1, color = my_red, linestyle='-', linewidth = 4)

	#ax_K_distribution.vlines(2.3e-1, 1e-6, 1, color = my_red, linestyle='--', linewidth = 4)	
	ax_K_distribution.vlines(np.mean(Z_final)/(C/Kd_r_renorm), 1e-6, 1, color = my_red, linestyle='-', linewidth = 4)	
	#ax_K_distribution2.vlines(2.3e-1, 1e-6, 1, color = my_red, linestyle='--', linewidth = 4)		
	ax_K_distribution2.vlines(np.mean(Z_final)/(C/Kd_r_renorm), 1e-6, 1, color = my_red, linestyle=':', linewidth = 2)		
	print('mean_Z = ', '%.3e'%np.mean(Z_final))

# Printing K from Gumbel
Nb = C
#Z_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
Z_array = np.flip((Nb/1)/Kds)
p_E = P_min_e_Q0(N_r, Q0, dE)
#p_E = P_min_e(N_r, avg_E, var_E, Es[:-1], dE)
print('N=', sum(p_E*dE))
p_E/=np.sum(p_E*dE)
p_Z = np.flip(p_E)/Z_array
#p_Z /= np.sum(p_Z[:-1]*np.diff(Z_array))
#p_Z = p_Z/np.sum(np.flip(p_Z[:-1])/np.flip(Z_array[:-1])*abs(np.diff(np.flip(Z_array))))
print('N = ',np.sum((p_Z[:-1])*(np.diff((Z_array)))))

p_Z /= np.sum(p_Z[:-1]*np.diff(Z_array))
ax_K_distribution.plot((((Z_array[:-1]))/(C/Kd_r_renorm)), (p_Z[:-1]), linestyle = '-', marker = '',  color = 'grey', ms = 2, linewidth = 3, alpha = .5, label = 'Gumbel')
ax_K_distribution2.plot((((Z_array[:-1]))/(C/Kd_r_renorm)), 1-np.cumsum((abs(p_Z[:-1]))*(np.diff((Z_array)))), linestyle = '--', marker = '',  color = 'grey', ms = 2, linewidth = 3, alpha = .5)
#ax_K_distribution2.plot((((Z_array[:]))/(C/Kd_r_renorm)), 1-np.cumsum(np.flip(p_E[:]*dE)), linestyle = '-', marker = '',  color = 'darkred', ms = 2, linewidth = 3, alpha = .5)

Z_array_tail = 10**np.linspace(9, 14.5, 50)

alpha = (lambda_A*beta_r)/(lambda_B*p + lambda_A)
zeta = (lambda_B*p)/(lambda_A*beta_r)

print('exponent = ', -1/zeta)

fit2p = (Z_array_tail**(-1/zeta-1))/((1e11)**(-1/zeta-1))
fit2p *= (1e9*np.diff(P_Z_sim/Counter)/np.diff(Zs_sim[:-1]))[((Zs_sim[:-2]))<1e11][-1]

fit2p2 = (Z_array_tail**(-1*alpha-1))/((1e11)**(-1*alpha-1))
fit2p2 *= (1e9*np.diff(P_Z_sim2/Counter)/np.diff(Zs_sim2[:-1]))[((Zs_sim2[:-2]))<1e11][-1]

ax_K_distribution.plot(Z_array_tail/(C/Kd_r_renorm), fit2p, linewidth = 3, color = 'darkorange', linestyle = '-', alpha = .6)
ax_K_distribution.plot(Z_array_tail/(C/Kd_r_renorm), fit2p2, linewidth = 3, color = 'tab:cyan', linestyle = '-')

fit2 = Z_array_tail**(-(1/(zeta)))/((1e11)**(-(1/(zeta))))
fit2 = Z_array_tail**(-2)/((1e11)**(-2))
fit2 *= (1-P_Z_sim/Counter)[((Zs_sim[:-1]))<1e11][-1]

fit22 = (Z_array_tail**(-1*alpha))/((1e11)**(-1*alpha))
#fit22 *= (1-np.cumsum(P_Z_sim*np.diff(Zs_sim2)))[((Zs_sim2[:-1]))<10**10.8][-1]
fit22 *= (1-P_Z_sim2/Counter)[((Zs_sim2[:-1]))<1e11][-1]

ax_K_distribution2.plot(Z_array_tail/(C/Kd_r_renorm), fit2, linewidth = 3, color = 'grey', linestyle = '-', alpha = .6)
#ax_K_distribution2.plot(Z_array_tail/(C/Kd_r_renorm), fit22, linewidth = 3, color = 'tab:cyan', linestyle = '-')


my_plot_layout(ax = ax_K_distribution, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K_distribution.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_K_distribution.set_ylim(bottom = 1e-8, top = 2)
ax_K_distribution.set_xlim(left = 1e-2, right = 10**(13.5)/(C/Kd_r_renorm))
#ax_K_distribution.set_xticks([])
#ax_K_distribution.set_yticks([])
#ax_K_distribution.set_yticklabels([1, 0.1, 0.01])
fig_K_distribution.savefig('../../Figures/1_Dynamics/Ensemble/L%d/Z_elite_P_'%L+energy_model+'_p-%.1f.pdf'%(p))

my_plot_layout(ax = ax_K_distribution2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K_distribution2.legend(fontsize = 28, title_fontsize = 30, loc = 0)
ax_K_distribution2.set_ylim(bottom = 1.5e-6, top = 1.05)
ax_K_distribution2.set_xlim(left = 1e-1, right = 7e1)
#ax_K_distribution2.set_xticks([])
#ax_K_distribution2.set_yticks([])
#ax_K_distribution2.set_yticklabels([1, 0.1, 0.01])
fig_K_distribution2.savefig('../../Figures/1_Dynamics/Ensemble/L%d/Z_elite_F_'%L+energy_model+'_p-%.1f.pdf'%(p))


print((C/Kd_r_renorm))
