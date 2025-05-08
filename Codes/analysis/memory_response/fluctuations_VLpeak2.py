import sys
sys.path.append('../../lib/')
from funcs import*

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

project = 'memory_response'
subproject = 'multi-epitope'
subproject = 'Panel'
subproject = 'fluctuations_potency'
experiment = 1
output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
os.makedirs(output_plot, exist_ok=True)

#--------------- PARAMETERS ---------------------

L0s = [1e7]#, 1e8]


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
antigen_seq = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18'
L=len(antigen_seq)

print('--------')
print('l=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
antigen = from_aa_to_i(antigen_seq, energy_model, '../../')
motif = get_motif(antigen, energy_model, '../../')
print('min_e_PWM=%.2f'%(np.sum([np.min(motif[:,i]) for i in range(len(motif[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(motif[:,i]) for i in range(len(motif[0,:]))])))

#Change values by the minimum
for i in np.arange(L):
    motif[:,i]-=np.min(motif[:,i], axis=0)

min_E = np.sum([np.min(motif[:,i]) for i in range(len(motif[0,:]))])
avg_E = (np.sum([np.mean(motif[:,i]) for i in range(len(motif[0,:]))]))
print(min_E, avg_E, E_m)
var_E = (np.sum([np.var(motif[:,i]) for i in range(len(motif[0,:]))]))
#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.001, 80, 1000000, motif, E_m, L)
Kds = np.exp(Es[:-1])


#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_K_distribution2_total, ax_N_distribution2_total = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

Kd_r_renorm = Kds[(P_min_e_Q0(1e8, Q0, dE))==np.max(P_min_e_Q0(1e8, Q0, dE))][0]
my_markers = ['o', 'D']
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
		N_ensss = [[1000+i for i in range(0, 6, 1)]]# + [2000+i for i in range(0, 1, 1)] + [4000+i for i in range(0, 1, 1)] + [8000+i for i in range(0, 1, 1)]] # L0 1e9
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
	fig_K_distribution, ax_N_distribution = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
	fig_K_distribution2, ax_N_distribution2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})

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
		# print('L0:%.0e and Z0:%.2e'%(L0, C/Kd_r))

		N_enss = N_ensss[i_p]

		#K = np.zeros_like(time_array)

		zs = []
		zs_largest_N = []
		zs_smallest_K = []
		zs_largest_z = []
		Zs = []
		counter = 0
		for N_ens in N_enss:
			print('N_ens = %d'%N_ens)

			#-----------------Loading data----------------------------
			return_data_type = 0
			data, return_data_type = get_data(folder_path = '../../out/memory_response', data_type = 'Z_Nens-%d_p-%.1f_L0-%.0e'%(N_ens, p, L0))
			
			if(return_data_type==1):
				Zs_i = data[4]
				counter_i = data[5]
			else:
				data = pd.read_csv(Text_files_path + 'primary_response/simulations/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L0, p, k_step, E_lim, t_lim, E_m))    
				Zs_i = []
				counter_i = 0

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
						Potencies = (clone_sizes_C[:, -1]-1)/Kds_C
						Zs_i += [np.sum(Potencies)]
						zs_i += list(Potencies)
						zs_largest_N_i += [Potencies[np.argmax(clone_sizes_C[:, -1])]]
						zs_smallest_K_i += [Potencies[np.argmin(Kds_C)]]
						zs_largest_z_i += [Potencies[np.argmax(Potencies)]]
						counter_i+=1

				f = open('../../out/memory_response/processed_data_Z_Nens-%d_p-%.1f_L0-%.0e.pkl'%(N_ens, p, L0), 'wb')
				pickle.dump([zs_i, zs_largest_N_i, zs_smallest_K_i, zs_largest_z_i, Zs_i, counter_i], f, pickle.HIGHEST_PROTOCOL)	

			Zs += Zs_i
			counter += counter_i
		
		ab = 1e-11
		zeta = (lambda_B*p)/(lambda_A*beta_r)
		lams = [1, 3.0, 5]
		lam_c = 1/zeta - 1
		lam_c2 = beta_r - 1
		my_colors = [my_red, my_green, my_brown, my_purple2, my_yellow, my_cyan, my_purple]
		for i_lam, lam in enumerate(lams):
			Ns = (ab*np.array(Zs))**(-lam)
			bins0 = np.logspace(np.min(np.log10(np.array(Ns)/(1))), np.max(np.log10(np.array(Ns)/(1))), 100)

			P_data_sim0 = plt.hist((np.array(Ns)/(1)), bins = bins0, density = False, cumulative  = True, alpha = 0)
			# ax_N_distribution.hist((np.array(Zs)/(C/Kd_r_renorm)), bins = bins0, density = True, cumulative  = False, alpha = .2)
			#P_data_sim = np.histogram(np.log10(np.array(Z_final)), bins = bins, density = True)

			Ns_sim0 = P_data_sim0[1]
			F_N_sim0 = P_data_sim0[0]
			normalization0 = len(Ns)

			pN = np.diff(F_N_sim0/normalization0)/np.diff(Ns_sim0[:-1])

			print('Normalization p(N):', '%.2f'%(np.sum(pN*np.diff(Ns_sim0[:-1]))))
			
			ax_N_distribution.plot(Ns_sim0[:-2], pN, color = my_red, markeredgecolor = my_colors[i_lam], markerfacecolor="None", linestyle='', marker = 'o', linewidth = 3, ms = 3, alpha = .9, label = r'$%.1f$'%(lam), zorder = 20)

			Z0_array = (Ns_sim0[:-1])**(-1/lam)/ab

			ax_N_distribution2.plot(Z0_array/	(C/Kd_r_renorm), (F_N_sim0/normalization0), color = my_green, markeredgecolor = my_colors[i_lam], markerfacecolor="None", linestyle='', marker = 'o', linewidth = 5, ms = 3, alpha = .9, label = r'$\lambda_A/\lambda_Z = %.1f$'%(lam), zorder = 20)

			# ax_N_distribution2.plot(Z0_array, 1 - 1/(1+(Z0_array/(1.7e10))**-2.5), color = my_green, markeredgecolor = my_colors[i_lam], markerfacecolor="None", linestyle='--', marker = '', linewidth = 3, ms = 3, alpha = .9, label = r'$\textrm{Hill}\, \textrm{fit}$', zorder = 20)
			# ax_N_distribution2.plot(Z0_array, 1 - 1/(1+(Z0_array*ab)**-1.2), color = my_blue2, markeredgecolor = my_colors[i_lam], markerfacecolor="None", linestyle='--', marker = '', linewidth = 3, ms = 3, alpha = .9, label = r'$in\, vitro$', zorder = 20)

			ax_N_distribution2_total.plot(((Ns_sim0[:-1])), 1-F_N_sim0/normalization0, markeredgecolor = my_colors[i_lam], markerfacecolor="None", linestyle='', marker = my_markers[i_L0], linewidth = 4, ms = 3, alpha = .8, label = r'$p=%d'%(p), zorder = 20)

		#  -------------------------------- Printing K from Gumbel --------------------------------
		Nb = C
		dE
		Z_array = np.flip((Nb/1)/Kds/(C/Kd_r_renorm))[-30000:]
		p_E = P_min_e_Q0(L0, Q0, dE)
		#p_E = P_min_e(L0, avg_E, var_E, Es[:-1], dE)
		print('N =', sum(p_E*dE))
		p_E/=np.sum(p_E*dE)
		p_Z = np.flip(p_E)[-30000:]/Z_array
		#p_Z /= np.sum(p_Z[:-1]*np.diff(Z_array))
		#p_Z = p_Z/np.sum(np.flip(p_Z[:-1])/np.flip(Z_array[:-1])*abs(np.diff(np.flip(Z_array))))
		print('N =',np.sum((p_Z[:-1])*(np.diff((Z_array)))))
		print('mean_Z_gumbel = %.3e'%(np.sum((Z_array[:-1]*p_Z[:-1])*(np.diff((Z_array))))))
		

		p_Z /= np.sum(p_Z[:-1]*np.diff(Z_array))
		# ax_N_distribution.plot((((Z_array[-20000:-1]))), (p_Z[-20000:-1]), linestyle = '--', marker = '',  color = 'grey', ms = 2, linewidth = 2, alpha = 1, label = r'$\textrm{Gumbel}$')
		# ax_N_distribution2.plot((((Z_array[-20000:-1]))), 1-np.cumsum((abs(p_Z[-20000:-1]))*(np.diff((Z_array[-20000:])))), linestyle = '--', marker = '',  color = 'grey', ms = 2, linewidth = 2, alpha = 1, label = r'$\textrm{Gumbel}$')
		
		# if L0 == 1e7:
			# ax_N_distribution2_total.plot((((Z_array[-20000:-1]))), 1-np.cumsum((abs(p_Z[-20000:-1]))*(np.diff((Z_array[-20000:])))), linestyle = '--', marker = '',  color = 'grey', ms = 2, linewidth = 2, alpha = 1, label = r'$\textrm{Gumbel}$')
		#ax_N_distribution2.plot((((Z_array[:]))/(C/Kd_r_renorm)), 1-np.cumsum(np.flip(p_E[:]*dE)), linestyle = '-', marker = '',  color = 'darkred', ms = 2, linewidth = 3, alpha = .5)
		#  -------------------------------- ---------------------- --------------------------------

		N_array_fit = 10**(np.linspace(-3.8, 4, 50))

		# -------- fits for p --------
		p_fit = (N_array_fit**(-1 + 1/(lam*zeta)))/((1)**(-1 + 1/(lam*zeta)))
		p_fit *= pN[((Ns_sim0[:-2]))<1][-1]

		ax_N_distribution.plot(N_array_fit, p_fit, ls = '--', color = 'k', label = r'${N_A^c}^{-1 + 1/\lambda\zeta}$')

		p_fit = (N_array_fit**(-1 + beta_r/(lam)))/((1)**(-1 - 1/lam + beta_r/(lam)))
		p_fit *= pN[((Ns_sim0[:-2]))<1][-1]

		# ax_N_distribution.plot(N_array_fit, p_fit, ls = '--', color = 'grey', label = r'$N^{-1 + \beta^*/\lambda}$')

		alpha = (lambda_A*beta_r)/(lambda_B*p + lambda_A)
		
		exp_potency = (lambda_A*beta_r/p)/(lambda_B + lambda_A/p)

		print('EXPONENTS:', 1/zeta, exp_potency)


		# -------- Gaussian fits --------
		if L0 == 1e7:
			varZ = 0.35
			norm_gauss = 2e3
			shift = 1
		elif L0 == 1e9:
			varZ = 1.80
			norm_gauss = 9e3
			shift = 1

		# ax_N_distribution.plot((Z_array), 1e-3*norm_gauss*np.exp(-(np.log((Z_array)) - np.log(shift*meanZ))**2/(2*np.sqrt(varZ))), color = my_blue, linestyle='-', marker = '', linewidth = 1, ms = 3, alpha = .8, label = r'$\textrm{gaussian fit}$')


		my_plot_layout(ax = ax_N_distribution, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		# ax_N_distribution.set_ylim(bottom = 1e-7, top = np.max(pN)*1.2)
		# ax_N_distribution.set_xlim(left = 1e-2, right = 1e2)
		#ax_N_distribution.set_xticks([])
		#ax_N_distribution.set_yticks([])
		#ax_N_distribution.set_yticklabels([1, 0.1, 0.01])
		ax_N_distribution.legend(fontsize = 10, title = r'$\lambda/\lambda_c$')
		fig_K_distribution.savefig(output_plot + '/VLpeak_p_'+energy_model+'_p-%.1f_L0-%.0e_linear_2.pdf'%(p, L0))

		my_plot_layout(ax = ax_N_distribution, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		# ax_N_distribution.set_ylim(bottom = np.min(pN[pN!=0])/5, top = np.max(pN)*2)
		# ax_N_distribution.set_xlim(left = 1e-2, right = 1e2)
		#ax_N_distribution.set_xticks([])
		#ax_N_distribution.set_yticks([])
		#ax_N_distribution.set_yticklabels([1, 0.1, 0.01])
		ax_N_distribution.legend(fontsize = 18, title = r'$\lambda_A/\lambda_Z$', loc = 0, title_fontsize = 20)
		fig_K_distribution.savefig(output_plot + '/VLpeak_p_'+energy_model+'_p-%.1f_L0-%.0e_log_2.pdf'%(p, L0))

		my_plot_layout(ax = ax_N_distribution2, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30)
		#ax_N_distribution2.legend(fontsize = 28, title_fontsize = 30, loc = 0)
		ax_N_distribution2.set_ylim(bottom = 0, top = 1.1)
		ax_N_distribution2.set_xlim(left = 5e-2, right = 1e2)
		#ax_N_distribution2.set_xticks([])
		#ax_N_distribution2.set_yticks([])
		#ax_N_distribution2.set_yticklabels([1, 0.1, 0.01])
		ax_N_distribution2.legend(fontsize = 20, title_fontsize = 16, loc = 0)#, title = r'$\lambda_A/\lambda_Z$')
		fig_K_distribution2.savefig(output_plot + '/VLpeak_F_'+energy_model+'_p-%.1f_L0-%.0e_linear_2.pdf'%(p, L0))

		my_plot_layout(ax = ax_N_distribution2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		#ax_N_distribution2.legend(fontsize = 28, title_fontsize = 30, loc = 0)
		# ax_N_distribution2.set_ylim(bottom = 1e-5, top = 1.8)
		# ax_N_distribution2.set_xlim(left = 1e-2, right = 1e2)
		#ax_N_distribution2.set_xticks([])
		#ax_N_distribution2.set_yticks([])
		#ax_N_distribution2.set_yticklabels([1, 0.1, 0.01])
		ax_N_distribution2.legend(fontsize = 20, title_fontsize = 16, loc = 3)
		fig_K_distribution2.savefig(output_plot + '/VLpeak_F_'+energy_model+'_p-%.1f_L0-%.0e_log_2.pdf'%(p, L0))

my_plot_layout(ax = ax_N_distribution2_total, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_N_distribution2_total.legend(fontsize = 28, title_fontsize = 30, loc = 0)
# ax_N_distribution2_total.set_yticks([1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
ax_N_distribution2_total.set_ylim(bottom = 2e-6, top = 1.8)
# ax_N_distribution2_total.set_xlim(left = 1e-4, right = 2e2)
#ax_N_distribution2_total.set_xticks([])
#ax_N_distribution2_total.set_yticklabels([1, 0.1, 0.01])
fig_K_distribution2_total.savefig(output_plot + '/VLpeak_F_'+energy_model+'_2.pdf')
