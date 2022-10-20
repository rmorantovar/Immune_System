import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ensss = [[200], [400, 300, 200, 100, 50], [200, 150, 100], [200, 100]]
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

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2.5, 3.0, 4.0]
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
color_list = np.array([my_red, my_green, my_blue2, my_gold])
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
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_NC, ax_NC = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_distribution, ax_NC_distribution = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_distribution2, ax_NC_distribution2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_scatter_affinity, ax_NC_scatter_affinity = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_scatter_clone_size, ax_NC_scatter_clone_size = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_scatter_avidity, ax_NC_scatter_avidity = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_avidity_order, ax_NC_avidity_order = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_avidity_cumulative, ax_NC_avidity_cumulative = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

NC_final_best_renorm = []
counter_total = 0
for i_kappa, kappa in enumerate(kappas):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0] 
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	beta_act = np.min([beta_r, beta_kappa])

	N_enss = N_ensss[i_kappa]

	NC_all = np.zeros_like(time)
	NC_final_all = []
	Counter_all = 0

	NC_best = np.zeros_like(time)
	NC_final_best = []
	Counter_best = 0

	NC_biggest = np.zeros_like(time)
	NC_final_biggest = []
	Counter_biggest = 0

	NC_avidity = np.zeros_like(time)
	NC_final_avidity = []
	Counter_avidity = 0

	Avidity_orders = []

	for N_ens in N_enss:
		print('N_ens = %d'%N_ens)

		#-----------------Loading data----------------------------
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
		#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
		

		for i_ens in tqdm(np.arange(N_ens)):
			data_i = data.loc[data[4]==i_ens]
			data_active = data_i.loc[data_i[1]==1]
			t_act_data = np.min(data_active[3])
			data_active = data_active.loc[data_active[3]<(t_act_data+1.0+0.1*(kappa-1))]
			
			data_active_all = data_active#.loc[data_active[3]<(t_act_theory)]
			activation_times_all = np.array(data_active_all[3])
			energies_all = np.array(data_active_all[0])

			#data_active_common = data_active.loc[data_active[3]>(t_act_theory)]
			#activation_times_common = np.array(data_active_common[3])
			#energies_common = np.array(data_active_common[0])

			#---------------------------- B cell linages ----------------------
			clone_sizes_all = get_clones_sizes_C(len(activation_times_all), time, activation_times_all, lambda_B, C, dT)
			#clone_sizes_common = get_clones_sizes_C(len(activation_times_common), time, activation_times_common, lambda_B, C, dT)

			#--------------------------t_C filter-------------------------
			lim_size = 2
			clone_sizes_C_all, activation_times_C_all, energies_C_all, filter_C_all, n_C_all = apply_filter_C(clone_sizes_all, activation_times_all, energies_all, lim_size)
			#clone_sizes_C_common, activation_times_C_common, energies_C_common, filter_C_common, n_C_common = apply_filter_C(clone_sizes_common, activation_times_common, energies_common, lim_size)
			#-------Simulations-------
			Kds_C_all = np.exp(energies_C_all)
			Avidities = np.divide(((clone_sizes_C_all-1).T)/N_A, Kds_C_all).T
			final_potencies = (Avidities[:,-1])
			final_potency = (np.sum(Avidities[:,-1]))
			sort_inds_avidity = np.flip(final_potencies.argsort())

			cum_avidity_freq = np.cumsum(final_potencies[sort_inds_avidity])/final_potency

			ax_NC_avidity_cumulative.plot(np.arange(len(final_potencies))+1, cum_avidity_freq, color = colors_kappa[i_kappa], alpha = .5 )

			order = 0
			K_i = final_potencies[sort_inds_avidity[order]]
			while (((K_i)/(final_potency))<.9):
				order+=1
				K_i+=final_potencies[sort_inds_avidity[order]]

			Avidity_orders.append(order+1)

			#Kds_C_common = np.exp(energies_C_common)

			NC_i_all = np.log(np.sum(Avidities, axis = 0))
			NC_i_best = np.log(np.sum(Avidities[energies_C_all==np.min(energies_C_all),:], axis = 0))
			NC_i_biggest = np.log(np.sum(Avidities[clone_sizes_C_all[:,-1]==np.max(clone_sizes_C_all[:,-1]),:], axis = 0))
			NC_i_avidity = np.log(np.sum(Avidities[Avidities[:,-1]==np.max(Avidities[:,-1]),:], axis = 0))
			
			NC_final_best_renorm.append(np.log((C/N_A)/np.min(Kds_C_all)))

			if(np.sum(~np.isinf(NC_i_all))!=0):
				NC_all += NC_i_all
				NC_final_all.append(NC_i_all[-1])
				Counter_all+=1

				NC_best += NC_i_best
				NC_final_best.append(NC_i_best[-1])
				Counter_best+=1

				NC_biggest += NC_i_biggest
				NC_final_biggest.append(NC_i_biggest[-1])
				Counter_biggest+=1

				NC_avidity += NC_i_avidity
				NC_final_avidity.append(NC_i_avidity[-1])
				Counter_avidity+=1

			
			#if(i_ens%1==0):
			#	ax_NC.plot(time, NC_i, color = colors_kappa[i_kappa], alpha = .1, linewidth = 1)
			
	counter_total+=Counter_all

	ax_NC_avidity_order.hist(Avidity_orders, color = colors_kappa[i_kappa], bins = np.logspace(0, np.log10(np.max(2*Avidity_orders)), 20))

	print(Counter_all, Counter_best, Counter_biggest, Counter_avidity)

	NC_all = (NC_all/Counter_all)
	NC_best = (NC_best/Counter_best)
	NC_biggest = (NC_biggest/Counter_biggest)
	NC_avidity = (NC_avidity/Counter_avidity)

	if(kappa==1):
		normalization_all = NC_all[-1]
		normalization_best = NC_all[-1]
		normalization_biggest = NC_all[-1]
		normalization_avidity = NC_all[-1]
	
	else:
		ax_NC.plot(time, NC_all - normalization_all, color = colors_kappa[i_kappa], alpha = 1, linewidth = 5, linestyle = '-', label = 'all')
		#ax_NC.plot(time, NC_best - normalization_best, color = colors_kappa[i_kappa], alpha = 1, linewidth = 5, linestyle = '--', label = 'affinity')
		#ax_NC.plot(time, NC_biggest - normalization_biggest, color = colors_kappa[i_kappa], alpha = 1, linewidth = 5, linestyle = ':', label = 'clone-size')
		ax_NC.plot(time, NC_avidity - normalization_avidity, color = colors_kappa[i_kappa], alpha = 1, linewidth = 5, linestyle = '-.', label = 'potency')

	print(np.max(np.array(NC_final_all) - normalization_all))
	print(np.max(np.array(NC_final_best) - normalization_best))
	print(np.max(np.array(NC_final_biggest) - normalization_biggest))
	print(np.max(np.array(NC_final_avidity) - normalization_avidity))

	#bins = np.linspace(-5, 7, 100)
	NC_data_all = np.histogram((NC_final_all) - normalization_all, bins = 'auto', density = False)
	NC_data_best = np.histogram((NC_final_best) - normalization_best, bins = np.linspace(-7, 7, 80), density = False)
	NC_data_biggest = np.histogram((NC_final_biggest) - normalization_biggest, bins = np.linspace(np.min((NC_final_biggest) - normalization_biggest), 7, 80), density = False)
	NC_data_avidity = np.histogram((NC_final_avidity) - normalization_avidity, bins = np.linspace(-7, 7, 80), density = False)
	
	ax_NC_scatter_affinity.scatter((NC_final_all) - normalization_all, (NC_final_best) - normalization_best, color = colors_kappa[i_kappa])
	ax_NC_scatter_clone_size.scatter((NC_final_all) - normalization_all, (NC_final_biggest) - normalization_biggest, color = colors_kappa[i_kappa])
	ax_NC_scatter_avidity.scatter((NC_final_all) - normalization_all, (NC_final_avidity) - normalization_avidity, color = colors_kappa[i_kappa])

	ax_NC_scatter_affinity.plot(np.linspace(-2, 6, 50), np.linspace(-2, 6, 50), color = 'black', alpha = .8)
	ax_NC_scatter_clone_size.plot(np.linspace(-2, 6, 50), np.linspace(-2, 6, 50), color = 'black', alpha = .8)
	ax_NC_scatter_avidity.plot(np.linspace(-2, 6, 50), np.linspace(-2, 6, 50), color = 'black', alpha = .8)

	ax_NC_distribution.plot(NC_data_all[1][:-1], NC_data_all[0]/Counter_all, color = colors_kappa[i_kappa], linestyle='-', marker = '', linewidth = 2, label = 'all')
	ax_NC_distribution.plot(NC_data_best[1][:-1], NC_data_best[0]/Counter_best, color = colors_kappa[i_kappa], linestyle='--', marker = '', linewidth = 2, label = 'affinity')
	ax_NC_distribution.plot(NC_data_biggest[1][:-1], NC_data_biggest[0]/Counter_biggest, color = colors_kappa[i_kappa], linestyle=':', marker = '', linewidth = 2, label = 'clone-size')
	ax_NC_distribution.plot(NC_data_avidity[1][:-1], NC_data_avidity[0]/Counter_avidity, color = colors_kappa[i_kappa], linestyle='-.', marker = '', linewidth = 2, label = 'potency')

	if(kappa!=1.0):
		ax_NC_distribution2.plot(NC_data_all[1][:-1], 1-np.cumsum(NC_data_all[0]/Counter_all), color = colors_kappa[i_kappa], linestyle='', marker = 'o', linewidth = 2, label = r'$%.1f$'%(kappa))
		#ax_NC_distribution2.plot(NC_data_best[1][:-1], 1-np.cumsum(NC_data_best[0]/Counter_best), color = colors_kappa[i_kappa], linestyle='', marker = '^', linewidth = 2, label = 'affinity')
		#ax_NC_distribution2.plot(NC_data_biggest[1][:-1], 1-np.cumsum(NC_data_biggest[0]/Counter_biggest), color = colors_kappa[i_kappa], linestyle='', marker = '*', linewidth = 2, label = 'clone-size')
		#ax_NC_distribution2.plot(NC_data_avidity[1][:-1], 1-np.cumsum(NC_data_avidity[0]/Counter_avidity), color = colors_kappa[i_kappa], linestyle='', marker = 's', linewidth = 2, label = 'potency')
	
	#Nb = np.exp(lambda_B*Tf)*((k_on*N_c)/(lambda_A*N_A))**(lambda_B/lambda_A)*(k_pr/k_on)**(kappa*lambda_B/lambda_A)*Kds**(-kappa*lambda_B/lambda_A)

Counter_best_renorm = counter_total
print(Counter_best_renorm)
#Counter_best_renorm = 1
NC_best_renorm_data = np.histogram((NC_final_best_renorm) - normalization_best, bins = 'auto', density = False)
ax_NC_distribution.plot(NC_best_renorm_data[1][:-1], NC_best_renorm_data[0]/Counter_best_renorm, color = 'gray', linestyle=':', marker = '', linewidth = 2)
#ax_NC_distribution2.plot(NC_best_renorm_data[1][:-1], 1-np.cumsum(NC_best_renorm_data[0]/Counter_best_renorm*np.diff(NC_best_renorm_data[1])), color = 'gray', linestyle='', marker = 'D', linewidth = 2)
ax_NC_distribution2.plot(NC_best_renorm_data[1][:-1], 1-np.cumsum(NC_best_renorm_data[0]/Counter_best_renorm), color = 'gray', linestyle='', marker = 'D', linewidth = 2)

# Printing K from Gumbel
Nb = C
#NC_array = np.log(1/(1+(Kds/((AA*(Nb))/N_A))))
NC_array = np.log((Nb/N_A)/Kds)
print(NC_array)
p_NC = P_min_e_Q0(N_r, Q0, dE)#/NC_array**2*(Nb/N_A)
print(np.sum(P_min_e_Q0(N_r, Q0, dE)*dE), np.sum(P_min_e_Q0(N_r, Q0, dE)[:-1]*abs(np.diff((NC_array)))))
p_NC = p_NC/np.sum(np.flip(p_NC[:-1])*abs(np.diff(np.flip(NC_array))))
ax_NC_distribution.plot((np.flip(NC_array[:-1]))-normalization_all, np.flip(p_NC[:-1]), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 2, alpha = .8, label = 'Gumbel')
ax_NC_distribution2.plot((np.flip(NC_array[:-1]))-normalization_all, 1-np.cumsum(np.flip(p_NC[:-1])*abs(np.diff(np.flip(NC_array)))), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 4, alpha = .8, label = 'Gumbel')

NC_array_tail = np.linspace(3, 7, 50)
exponent_tail  = beta_act+1
fit_tail = np.exp(-exponent_tail*(NC_array_tail))/np.exp(-exponent_tail*(5.5))
fit_tail *= (1-np.cumsum(np.flip(p_NC[:-1])*abs(np.diff(np.flip(NC_array)))))[(np.flip(NC_array)[:-1]-normalization_all)<5.5][-1]
print((1-np.cumsum(np.flip(p_NC[:-1])*abs(np.diff(np.flip(NC_array)))))[(np.flip(NC_array)[:-1]-normalization_all)<5.5][-1])
ax_NC_distribution2.plot(NC_array_tail, fit_tail, linewidth = 2, color = 'grey')

my_plot_layout(ax = ax_NC_distribution, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_distribution.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_NC_distribution.set_ylim(bottom = 2e-3, top = 1)
ax_NC_distribution.set_xlim(left = 0, right = 6.5)
#ax_NC_distribution.set_xticks([])
#ax_NC_distribution.set_yticks([])
#ax_NC_distribution.set_yticklabels([1, 0.1, 0.01])
fig_NC_distribution.savefig('../../Figures/1_Dynamics/Ensemble/NC_P_elite_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_distribution2, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_NC_distribution2.legend(fontsize = 28, title_fontsize = 30, loc = 3, title = r'$p$')
ax_NC_distribution2.set_ylim(bottom = 1e-4, top = 3)
ax_NC_distribution2.set_xlim(left = 1, right = 7)
#ax_NC_distribution2.set_xticks([])
#ax_NC_distribution2.set_yticks([])
#ax_NC_distribution2.set_yticklabels([1, 0.1, 0.01])
fig_NC_distribution2.savefig('../../Figures/1_Dynamics/Ensemble/NC_F_elite_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_scatter_affinity, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_scatter_affinity.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_NC_scatter_affinity.set_ylim(bottom = -4, top = 7)
ax_NC_scatter_affinity.set_xlim(left = -2, right = 7)
#ax_NC_scatter_affinity.set_xticks([])
#ax_NC_scatter_affinity.set_yticks([])
#ax_NC_scatter_affinity.set_yticklabels([1, 0.1, 0.01])
fig_NC_scatter_affinity.savefig('../../Figures/1_Dynamics/Ensemble/NC_scatter_affinity_elite_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_scatter_clone_size, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_scatter_clone_size.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_NC_scatter_clone_size.set_ylim(bottom = -4, top = 7)
ax_NC_scatter_clone_size.set_xlim(left = -2, right = 7)
#ax_NC_scatter_clone_size.set_xticks([])
#ax_NC_scatter_clone_size.set_yticks([])
#ax_NC_scatter_clone_size.set_yticklabels([1, 0.1, 0.01])
fig_NC_scatter_clone_size.savefig('../../Figures/1_Dynamics/Ensemble/NC_scatter_clone-size_elite_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_scatter_avidity, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_scatter_avidity.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_NC_scatter_avidity.set_ylim(bottom = -4, top = 7)
ax_NC_scatter_avidity.set_xlim(left = -2, right = 7)
#ax_NC_scatter_avidity.set_xticks([])
#ax_NC_scatter_avidity.set_yticks([])
#ax_NC_scatter_avidity.set_yticklabels([1, 0.1, 0.01])
fig_NC_scatter_avidity.savefig('../../Figures/1_Dynamics/Ensemble/NC_scatter_avidity_elite_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_avidity_order, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_avidity_order.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
#ax_NC_avidity_order.set_ylim(bottom = -4, top = 7)
#ax_NC_avidity_order.set_xlim(left = -2, right = 7)
#ax_NC_avidity_order.set_xticks([])
#ax_NC_avidity_order.set_yticks([])
#ax_NC_avidity_order.set_yticklabels([1, 0.1, 0.01])
fig_NC_avidity_order.savefig('../../Figures/1_Dynamics/Ensemble/NC_avidity-order_elite_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_avidity_cumulative, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_avidity_cumulative.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
#ax_NC_avidity_cumulative.set_ylim(bottom = -4, top = 7)
#ax_NC_avidity_cumulative.set_xlim(left = -2, right = 7)
#ax_NC_avidity_cumulative.set_xticks([])
#ax_NC_avidity_cumulative.set_yticks([])
#ax_NC_avidity_cumulative.set_yticklabels([1, 0.1, 0.01])
fig_NC_avidity_cumulative.savefig('../../Figures/1_Dynamics/Ensemble/NC_avidity-cum_elite_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_NC.legend(fontsize = 28, title_fontsize = 30, title = r'$p$')
ax_NC.set_xlim(left = 4.5, right = Tf)
ax_NC.set_ylim(bottom = -1, top = 3.5)
#ax_NC.set_yticks([1, 0.1, 0.01, 0.001])
#ax_NC.set_yticklabels([1, 0.1, 0.01])
fig_NC.savefig('../../Figures/1_Dynamics/Ensemble/NC_elite_'+energy_model+'.pdf')



