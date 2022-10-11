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
Tf = 12
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2, 3, 4]
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
fig_NC_scatter, ax_NC_scatter = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
for i_kappa, kappa in enumerate(kappas):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0] 
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)

	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
	#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
	data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
	

	NC_rare = np.zeros_like(time)
	#NC_common = np.zeros_like(time)
	NC_final_rare = []
	#NC_final_common = []
	Counter_rare = 0
	#Counter_common = 0 

	NC_best = np.ones_like(time)
	NC_final_best = []
	Counter_best = 0
	for i_ens in tqdm(np.arange(N_ens)):
		data_i = data.loc[data[4]==i_ens]
		data_active = data_i.loc[data_i[1]==1]
		t_act_data = np.min(data_active[3])
		data_active = data_active.loc[data_active[3]<(t_act_data+0.9+0.1*(kappa-1))]
		
		data_active_rare = data_active#.loc[data_active[3]<(t_act_theory)]
		activation_times_rare = np.array(data_active_rare[3])
		energies_rare = np.array(data_active_rare[0])

		#data_active_common = data_active.loc[data_active[3]>(t_act_theory)]
		#activation_times_common = np.array(data_active_common[3])
		#energies_common = np.array(data_active_common[0])

		#---------------------------- B cell linages ----------------------
		clone_sizes_rare = get_clones_sizes_C(len(activation_times_rare), time, activation_times_rare, lambda_B, C, dT)
		#clone_sizes_common = get_clones_sizes_C(len(activation_times_common), time, activation_times_common, lambda_B, C, dT)

		#--------------------------t_C filter-------------------------
		lim_size = 2
		clone_sizes_C_rare, activation_times_C_rare, energies_C_rare, filter_C_rare, n_C_rare = apply_filter_C(clone_sizes_rare, activation_times_rare, energies_rare, lim_size)
		#clone_sizes_C_common, activation_times_C_common, energies_C_common, filter_C_common, n_C_common = apply_filter_C(clone_sizes_common, activation_times_common, energies_common, lim_size)
		#-------Simulations-------
		Kds_C_rare = np.exp(energies_C_rare)
		#Kds_C_common = np.exp(energies_C_common)

		#NC_i_rare = np.log(1-np.array([np.product(1-1/(1+(Kds_C_rare/((AA*(clone_sizes_C_rare[:,t]-1))/N_A)))) for t in np.arange(len(time))]))
		NC_i_rare = [np.log(np.sum(((clone_sizes_C_rare[:,t]-1)/N_A)/Kds_C_rare)) for t in np.arange(len(time))]
		#NC_i_best = [np.sum(np.log(((clone_sizes_C_rare[:,t]-1)/N_A)/Kds_C_rare[:])) for t in np.arange(len(time))]
		NC_i_best = [np.log(np.sum(((clone_sizes_C_rare[energies_C_rare==np.min(energies_C_rare),t]-1)/N_A)/Kds_C_rare[energies_C_rare==np.min(energies_C_rare)])) for t in np.arange(len(time))]
		#NC_i_common = np.log(1-np.array([np.product(1-1/(1+(Kds_C_common/((AA*(clone_sizes_C_common[:,t]-1))/N_A)))) for t in np.arange(len(time))]))

		if(np.sum(~np.isinf(NC_i_rare))!=0):
			NC_rare += NC_i_rare
			NC_final_rare.append(NC_i_rare[-1])
			Counter_rare+=1

			NC_best += NC_i_best
			NC_final_best.append(NC_i_best[-1])
			Counter_best+=1
		# if(np.sum(~np.isinf(NC_i_common))!=0):
		# 	NC_common += NC_i_common
		# 	NC_final_common.append(NC_i_common[-1])
		# 	Counter_common+=1

		#if(i_ens%1==0):
		#	ax_NC.plot(time, NC_i, color = colors_kappa[i_kappa], alpha = .1, linewidth = 1)
	print(Counter_rare, Counter_best)
	NC_rare = (NC_rare/Counter_rare)
	NC_best = (NC_best/Counter_best)
	#NC_common = NC_common/Counter_common

	if(i_kappa==0):
		normalization_rare = NC_rare[-1]
		normalization_best = NC_rare[-1]
	
	ax_NC.plot(time, NC_rare - normalization_rare, color = colors_kappa[i_kappa], alpha = 1, label = r'$%d$'%kappa, linewidth = 5)
	#ax_NC.plot(time, NC_common - normalization_common, color = colors_kappa[i_kappa], alpha = .8, linewidth = 5, linestyle = 'dashed')
	ax_NC.plot(time, NC_best - normalization_best, color = colors_kappa[i_kappa], alpha = 1, linewidth = 5, linestyle = ':')

	print(np.min(np.array(NC_final_rare) - normalization_rare))
	print(np.min(np.array(NC_final_best) - normalization_best))

	NC_rare_data = np.histogram((NC_final_rare) - normalization_rare, bins = np.linspace(-7, 7, 28), density = False)
	#NC_common_data = np.histogram(np.array(NC_final_common) - normalization_common, bins = np.linspace(-7, 7, 50), density = False)
	NC_best_data = np.histogram((NC_final_best) - normalization_best, bins = np.linspace(-7, 7, 28), density = False)
	
	ax_NC_scatter.scatter((NC_final_rare) - normalization_rare, (NC_final_best) - normalization_best, color = colors_kappa[i_kappa])
	ax_NC_scatter.plot(np.linspace(-2, 6, 50), np.linspace(-2, 6, 50), color = 'black', alpha = .8)

	ax_NC_distribution.plot(NC_rare_data[1][:-1], NC_rare_data[0]/Counter_rare, color = colors_kappa[i_kappa], linestyle='-', marker = '', label = r'$%d$'%kappa, linewidth = 2)
	#ax_NC_distribution.plot(NC_common_data[1][:-1], NC_common_data[0]/Counter_common, color = colors_kappa[i_kappa], linestyle='--', marker = '', linewidth = 2)
	ax_NC_distribution.plot(NC_best_data[1][:-1], NC_best_data[0]/Counter_best, color = colors_kappa[i_kappa], linestyle=':', marker = '', linewidth = 2)

	ax_NC_distribution2.plot(NC_rare_data[1][:-1], 1-np.cumsum(NC_rare_data[0]/Counter_rare), color = colors_kappa[i_kappa], linestyle='-', marker = '', label = r'$%d$'%kappa, linewidth = 2)
	#ax_NC_distribution2.plot(NC_common_data[1][:-1], np.cumsum(NC_common_data[0]/Counter_common), color = colors_kappa[i_kappa], linestyle='--', marker = '', linewidth = 2)
	ax_NC_distribution2.plot(NC_best_data[1][:-1], 1-np.cumsum(NC_best_data[0]/Counter_best), color = colors_kappa[i_kappa], linestyle=':', marker = '', linewidth = 2)
	
	#Nb = np.exp(lambda_B*Tf)*((k_on*N_c)/(lambda_A*N_A))**(lambda_B/lambda_A)*(k_pr/k_on)**(kappa*lambda_B/lambda_A)*Kds**(-kappa*lambda_B/lambda_A)

# Printing K from Gumbel
Nb = C
#NC_array = np.log(1/(1+(Kds/((AA*(Nb))/N_A))))
NC_array = np.log((Nb/N_A)/Kds)
p_NC = P_min_e_Q0(N_r, Q0, dE)#/NC_array**2*(Nb/N_A)
print(np.sum(P_min_e_Q0(N_r, Q0, dE)*dE), np.sum(P_min_e_Q0(N_r, Q0, dE)[:-1]*abs(np.diff((NC_array)))))
p_NC = p_NC/np.sum(np.flip(p_NC[:-1])*abs(np.diff(np.flip(NC_array))))
ax_NC_distribution.plot((np.flip(NC_array[:-1]))-normalization_rare, np.flip(p_NC[:-1]), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 2, alpha = .8, label = 'Gumbel')
ax_NC_distribution2.plot((np.flip(NC_array[:-1]))-normalization_rare, 1-np.cumsum(np.flip(p_NC[:-1])*abs(np.diff(np.flip(NC_array)))), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 4, alpha = .8, label = 'Gumbel')

my_plot_layout(ax = ax_NC_distribution, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_distribution.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_NC_distribution.set_ylim(bottom = 2e-3, top = 1)
ax_NC_distribution.set_xlim(left = 1, right = 6.5)
#ax_NC_distribution.set_xticks([])
#ax_NC_distribution.set_yticks([])
#ax_NC_distribution.set_yticklabels([1, 0.1, 0.01])
fig_NC_distribution.savefig('../../Figures/1_Dynamics/Ensemble/NC_P_RC_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_distribution2, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_NC_distribution2.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_NC_distribution2.set_ylim(bottom = 1e-4)
ax_NC_distribution2.set_xlim(left = 1, right = 6.5)
#ax_NC_distribution2.set_xticks([])
#ax_NC_distribution2.set_yticks([])
#ax_NC_distribution2.set_yticklabels([1, 0.1, 0.01])
fig_NC_distribution2.savefig('../../Figures/1_Dynamics/Ensemble/NC_F_RC_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_scatter, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_scatter.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
#ax_NC_scatter.set_ylim(bottom = 1e-4)
#ax_NC_scatter.set_xlim(left = 1, right = 6.5)
#ax_NC_scatter.set_xticks([])
#ax_NC_scatter.set_yticks([])
#ax_NC_scatter.set_yticklabels([1, 0.1, 0.01])
fig_NC_scatter.savefig('../../Figures/1_Dynamics/Ensemble/NC_scatter_RC_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_NC.legend(fontsize = 30, title_fontsize = 32, title = r'$p$')
ax_NC.set_xlim(left = 4.5, right = Tf)
ax_NC.set_ylim(bottom = -1, top = 3.5)
#ax_NC.set_yticks([1, 0.1, 0.01, 0.001])
#ax_NC.set_yticklabels([1, 0.1, 0.01])
fig_NC.savefig('../../Figures/1_Dynamics/Ensemble/NC_RC_'+energy_model+'.pdf')




