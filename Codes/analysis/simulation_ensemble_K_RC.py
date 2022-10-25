import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 200
N_r = 2e8
T0 = 3
Tf = 13
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


kappas = [1, 2, 3, 4, 5]
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
color_list = np.array([my_red, my_green, my_blue2, my_gold, my_purple])
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
fig_K, ax_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_distribution, ax_K_distribution = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_distribution2, ax_K_distribution2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_scatter_common, ax_K_scatter_common = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_scatter_rare, ax_K_scatter_rare = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_scatter_biggest_affinity, ax_K_scatter_biggest_affinity = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
#fig_K_avidity_order, ax_K_avidity_order = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
#fig_K_avidity_cumulative, ax_K_avidity_cumulative = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

# Printing K from Gumbel
Nb = C
#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
K_array = ((Nb/1)/Kds)
p_K = P_min_e_Q0(N_r, Q0, dE)#/K_array**2*(Nb/1)
print(np.sum(P_min_e_Q0(N_r, Q0, dE)*dE), np.sum(P_min_e_Q0(N_r, Q0, dE)[:-1]/K_array[:-1]*abs(np.diff((K_array)))))
p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
ax_K_distribution.plot(((np.flip(K_array[:-1]))/1), np.flip(p_K[:-1]), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 2, alpha = .8, label = 'Gumbel')
ax_K_distribution2.plot(((np.flip(K_array[:-1]))/1), 1-np.cumsum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array)))), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 4, alpha = .8, label = 'Gumbel')

K_array_tail = 10**np.linspace(9, 14.5, 50)
exponent_tail  = beta_r+1
fit_tail =(K_array_tail**(-exponent_tail))/((10**13.2)**(-exponent_tail))
fit_tail *= (1-np.cumsum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array)))))[np.log10(np.flip(K_array)[:-1]/1)<13.2][-1]
ax_K_distribution2.plot(K_array_tail, fit_tail, linewidth = 2, color = 'black', linestyle = 'dashed')

K_final_best_renorm = []
counter_total = 0

for i_kappa, kappa in enumerate(kappas):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0] 
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	beta_act = np.min([beta_r, beta_kappa])

	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
	#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
	data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
	
	K_all = np.zeros_like(time)
	K_final_all = []
	Counter_all = 0

	K_common = np.zeros_like(time)
	K_final_common = []
	Counter_common = 0

	K_rare = np.zeros_like(time)
	K_final_rare = []
	Counter_rare = 0

	final_biggest = []
	final_biggest_affinity = []

	index_rare_cases = []
	index_common_cases = []

	for i_ens in tqdm(np.arange(N_ens)):
		data_i = data.loc[data[4]==i_ens]
		data_active = data_i.loc[data_i[1]==1]
		t_act_data = np.min(data_active[3])
		data_active = data_active.loc[data_active[3]<(t_act_data+1.0+.1*(kappa-1))]
		
		data_active_all = data_active#.loc[data_active[3]<(t_act_theory)]
		activation_times_all = np.array(data_active_all[3])
		energies_all = np.array(data_active_all[0])
		data_active_common = data_active.loc[data_active[3]>=(t_act_theory)]

		#---------------------------- B cell linages ----------------------
		clone_sizes_all = get_clones_sizes_C(len(activation_times_all), time, activation_times_all, lambda_B, C, dT)
		#--------------------------t_C filter-------------------------
		lim_size = 2
		clone_sizes_C_all, activation_times_C_all, energies_C_all, filter_C_all, n_C_all = apply_filter_C(clone_sizes_all, activation_times_all, energies_all, lim_size)
				#-------Simulations-------
		if(len(energies_C_all)>0):
			Kds_C_all = np.exp(energies_C_all)
			Avidities = np.divide(((clone_sizes_C_all-1).T)/1, Kds_C_all).T

			filter_common = activation_times_C_all>=t_act_theory
			Kds_C_common = np.exp(energies_C_all[filter_common])
			Avidities_common = np.divide(((clone_sizes_C_all[filter_common,:]-1).T)/1, Kds_C_all[filter_common]).T

			filter_rare = activation_times_C_all<t_act_theory
			Kds_C_rare = np.exp(energies_C_all[filter_rare])
			Avidities_rare = np.divide(((clone_sizes_C_all[filter_rare,:]-1).T)/1, Kds_C_all[filter_rare]).T

			# final_potencies = (Avidities[:,-1])
			# final_potency = (np.sum(Avidities[:,-1]))
			# sort_inds_avidity = np.flip(final_potencies.argsort())
			# cum_avidity_freq = np.cumsum(final_potencies[sort_inds_avidity])/final_potency
			# ax_K_avidity_cumulative.plot(np.arange(len(final_potencies))+1, cum_avidity_freq, color = colors_kappa[i_kappa], alpha = .5 )
			# order = 0
			# K_i = final_potencies[sort_inds_avidity[order]]
			# while (((K_i)/(final_potency))<.9):
			# 	order+=1
			# 	K_i+=final_potencies[sort_inds_avidity[order]]
			# Avidity_orders.append(order+1)

			K_i_all = np.sum(Avidities, axis = 0)
			K_i_common = np.sum(Avidities_common, axis = 0)
			K_i_rare = np.sum(Avidities_rare, axis = 0)
			
			K_final_best_renorm.append((C/1)/np.min(Kds_C_all))
			final_biggest_affinity.append(energies_C_all[clone_sizes_C_all[:,-1]==np.max(clone_sizes_C_all[:,-1])][0])
			final_biggest.append(np.max(clone_sizes_C_all[:,-1]))

			if(np.sum(K_i_all)!=0):
				K_all += K_i_all
				K_final_all.append(K_i_all[-1])
				Counter_all+=1
			if(np.sum(K_i_common)!=0):
				K_common += K_i_common
				K_final_common.append(K_i_common[-1])
				Counter_common+=1
				index_common_cases.append(int(i_ens))
			if(np.sum(K_i_rare)!=0):
				K_rare += K_i_rare
				K_final_rare.append(K_i_rare[-1])
				Counter_rare+=1
				index_rare_cases.append(int(i_ens))

			#if(i_ens%1==0):
			#	ax_K.plot(time, K_i, color = colors_kappa[i_kappa], alpha = .1, linewidth = 1)

	#ax_K_avidity_order.hist(Avidity_orders, color = colors_kappa[i_kappa], bins = 'auto')
	counter_total+=Counter_all

	print(Counter_all, Counter_common, Counter_rare)
	print(len(K_final_all), len(K_final_common), len(K_final_rare))

	K_all = (K_all/Counter_all)
	K_common = (K_common/Counter_common)
	K_rare = (K_rare/Counter_rare)

	if(i_kappa==0):
		normalization_all = 1# K_all[-1]
		normalization_common = 1# K_all[-1]
		normalization_rare = 1# K_all[-1]
	
	ax_K.plot(time, (K_all/normalization_all), color = colors_kappa[i_kappa], alpha = 1, linewidth = 5, linestyle = '-', label = r'$%d$'%kappa)
	ax_K.plot(time, (K_common/normalization_common), color = colors_kappa[i_kappa], alpha = 1, linewidth = 5, linestyle = '--')
	ax_K.plot(time, (K_rare/normalization_rare), color = colors_kappa[i_kappa], alpha = 1, linewidth = 5, linestyle = ':')

	print(np.log(np.max(np.array(K_final_all)/normalization_all)))
	print(np.log(np.max(np.array(K_final_common)/normalization_common)))
	print(np.log(np.max(np.array(K_final_rare)/normalization_rare)))

	K_data_all = np.histogram(np.log10(np.array(K_final_all)/normalization_all), bins = np.linspace(6, 14.5, 120), density = False)
	K_data_common = np.histogram(np.log10(np.array(K_final_common)/normalization_common), bins = np.linspace(6, 14.5, 120), density = False)
	K_data_rare = np.histogram(np.log10(np.array(K_final_rare)/normalization_rare), bins = np.linspace(6, 14.5, 120), density = False)

	ax_K_scatter_common.scatter((np.array(K_final_all)[index_common_cases]/normalization_all), (np.array((K_final_common))/normalization_common), color = colors_kappa[i_kappa], marker = '*', alpha = .8)
	ax_K_scatter_rare.scatter((np.array(K_final_all)[index_rare_cases]/normalization_all), (np.array((K_final_rare))/normalization_rare), color = colors_kappa[i_kappa], marker = 'o', alpha = .8)
	#ax_K_scatter_biggest_affinity.scatter(np.exp(final_biggest_affinity), np.array(final_biggest)/C, color = colors_kappa[i_kappa], marker = 'o', alpha = .8, label = r'$%.1f$'%kappa)
	#ax_K_scatter_biggest_affinity.hist2d(np.exp(final_biggest_affinity), np.array(final_biggest)/C, bins = [np.logspace(-9, 3, 20), np.linspace(0, 1, 10)], density = True, color = colors_kappa[i_kappa], label = r'$%.1f$'%kappa)
	
	df = pd.DataFrame({'x':np.exp(final_biggest_affinity), 'y':np.array(final_biggest)/C})
	#df = pd.DataFrame({'x':np.array(final_biggest)/C, 'y':np.array(final_biggest)/C})

	sns.kdeplot(data = df, x='x', y='y', ax=ax_K_scatter_biggest_affinity, color = colors_kappa[i_kappa], log_scale = (True, False), cut = 0, label = r'$%d$'%kappa)#, clip = ((np.min(np.exp(final_biggest_affinity)), np.max(np.exp(final_biggest_affinity))),(0, 1)))
	
	ax_K_scatter_common.plot(np.logspace(6, 14.5, 50), np.logspace(6, 14.5, 50), color = 'black', alpha = .8)
	ax_K_scatter_rare.plot(np.logspace(6, 14.5, 50), np.logspace(6, 14.5, 50), color = 'black', alpha = .8)

	ax_K_distribution.plot(10**(K_data_all[1][:-1]), K_data_all[0]/Counter_all, color = colors_kappa[i_kappa], linestyle='-', marker = '', linewidth = 2, label = r'$%d$'%kappa)
	ax_K_distribution.plot(10**(K_data_common[1][:-1]), K_data_common[0]/Counter_common, color = colors_kappa[i_kappa], linestyle='--', marker = '', linewidth = 2)
	ax_K_distribution.plot(10**(K_data_rare[1][:-1]), K_data_rare[0]/Counter_rare, color = colors_kappa[i_kappa], linestyle=':', marker = '', linewidth = 2)

	ax_K_distribution2.plot(10**(K_data_all[1][:-1]), 1-np.cumsum(K_data_all[0]/Counter_all), color = colors_kappa[i_kappa], linestyle='', marker = 'o', linewidth = 2, label = r'$%d$'%kappa)
	ax_K_distribution2.plot(10**(K_data_common[1][:-1]), 1-np.cumsum(K_data_common[0]/Counter_common), color = colors_kappa[i_kappa], linestyle='', marker = '^', linewidth = 2)
	ax_K_distribution2.plot(10**(K_data_rare[1][:-1]), 1-np.cumsum(K_data_rare[0]/Counter_rare), color = colors_kappa[i_kappa], linestyle='', marker = 's', linewidth = 2)
	
	#Nb = np.exp(lambda_B*Tf)*((k_on*N_c)/(lambda_A*N_A))**(lambda_B/lambda_A)*(k_pr/k_on)**(kappa*lambda_B/lambda_A)*Kds**(-kappa*lambda_B/lambda_A)

Counter_best_renorm = counter_total

K_best_renorm_data = np.histogram(np.log10(np.array(K_final_best_renorm)/normalization_common), bins = np.linspace(10, 14.5, 60), density = True)
ax_K_distribution.plot(10**(K_best_renorm_data[1][:-1]), K_best_renorm_data[0]/Counter_best_renorm, color = 'gray', linestyle=':', marker = '', linewidth = 2)
ax_K_distribution2.plot(10**(K_best_renorm_data[1][:-1]), 1-np.cumsum(K_best_renorm_data[0]*np.diff(K_best_renorm_data[1])), color = 'gray', linestyle='', marker = 'D', linewidth = 2)


my_plot_layout(ax = ax_K_distribution, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K_distribution.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_K_distribution.set_ylim(bottom = 2e-3, top = 1)
ax_K_distribution.set_xlim(left = 10**(6), right = 10**(14.5))
#ax_K_distribution.set_xticks([])
#ax_K_distribution.set_yticks([])
#ax_K_distribution.set_yticklabels([1, 0.1, 0.01])
fig_K_distribution.savefig('../../Figures/1_Dynamics/Ensemble/K_RC_P_'+energy_model+'.pdf')

my_plot_layout(ax = ax_K_distribution2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K_distribution2.legend(fontsize = 28, title_fontsize = 30, title = r'$p$', loc = 3)
ax_K_distribution2.set_ylim(bottom = 1e-5, top = 5)
ax_K_distribution2.set_xlim(left = 10**(11), right = 1*10**(14.5))
#ax_K_distribution2.set_xticks([])
#ax_K_distribution2.set_yticks([])
#ax_K_distribution2.set_yticklabels([1, 0.1, 0.01])
fig_K_distribution2.savefig('../../Figures/1_Dynamics/Ensemble/K_RC_F_'+energy_model+'.pdf')

my_plot_layout(ax = ax_K_scatter_common, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K_scatter_common.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_K_scatter_common.set_ylim(bottom = 10**(6), top = 10**(14.5))
ax_K_scatter_common.set_xlim(left = 10**(6), right = 10**(14.5))
#ax_K_scatter_common.set_xticks([])
#ax_K_scatter_common.set_yticks([])
#ax_K_scatter_common.set_yticklabels([1, 0.1, 0.01])
fig_K_scatter_common.savefig('../../Figures/1_Dynamics/Ensemble/K_RC_scatter_common_'+energy_model+'.pdf')

my_plot_layout(ax = ax_K_scatter_rare, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K_scatter_rare.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_K_scatter_rare.set_ylim(bottom = 10**(6), top = 10**(14.5))
ax_K_scatter_rare.set_xlim(left = 10**(6), right = 10**(14.5))
#ax_K_scatter_rare.set_xticks([])
#ax_K_scatter_rare.set_yticks([])
#ax_K_scatter_rare.set_yticklabels([1, 0.1, 0.01])
fig_K_scatter_rare.savefig('../../Figures/1_Dynamics/Ensemble/K_RC_scatter_rare_'+energy_model+'.pdf')


my_plot_layout(ax = ax_K_scatter_biggest_affinity, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K_scatter_biggest_affinity.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 0)
ax_K_scatter_biggest_affinity.set_ylim(bottom = -0.05, top = 1.05)
ax_K_scatter_biggest_affinity.set_xlim(left = 1e-9, right = 2e-4)
#ax_K_scatter_biggest_affinity.set_xticks([])
#ax_K_scatter_biggest_affinity.set_yticks([])
#ax_K_scatter_biggest_affinity.set_yticklabels([1, 0.1, 0.01])
fig_K_scatter_biggest_affinity.savefig('../../Figures/1_Dynamics/Ensemble/K_RC_scatter_big_affi_'+energy_model+'.pdf')


# my_plot_layout(ax = ax_K_avidity_order, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
# #ax_K_avidity_order.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
# #ax_K_avidity_order.set_ylim(bottom = 1e-4)
# #ax_K_avidity_order.set_xlim(left = 1, right = 6.5)
# #ax_K_avidity_order.set_xticks([])
# #ax_K_avidity_order.set_yticks([])
# #ax_K_avidity_order.set_yticklabels([1, 0.1, 0.01])
# fig_K_avidity_order.savefig('../../Figures/1_Dynamics/Ensemble/K_avidity-order_RC_'+energy_model+'.pdf')

# my_plot_layout(ax = ax_K_avidity_cumulative, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
# #ax_K_avidity_cumulative.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
# #ax_K_avidity_cumulative.set_ylim(bottom = 1e-4)
# #ax_K_avidity_cumulative.set_xlim(left = 1, right = 6.5)
# #ax_K_avidity_cumulative.set_xticks([])
# #ax_K_avidity_cumulative.set_yticks([])
# #ax_K_avidity_cumulative.set_yticklabels([1, 0.1, 0.01])
# fig_K_avidity_cumulative.savefig('../../Figures/1_Dynamics/Ensemble/K_avidity-cum_RC_'+energy_model+'.pdf')

my_plot_layout(ax = ax_K, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K.legend(fontsize = 28, title_fontsize = 30, title = r'$p$')
ax_K.set_xlim(left = 4.5, right = Tf)
ax_K.set_ylim(bottom = 2e8, top = 2e12)
#ax_K.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K.set_yticklabels([1, 0.1, 0.01])
fig_K.savefig('../../Figures/1_Dynamics/Ensemble/K_RC_'+energy_model+'.pdf')




