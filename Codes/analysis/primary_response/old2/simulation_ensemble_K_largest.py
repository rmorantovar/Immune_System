import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_enss = [500, 500, 500, 500]
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
color_list = np.array([my_red, my_blue2, my_green, my_gold, my_purple2])
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
fig_K_scatter_biggest_affinity, ax_K_scatter_biggest_affinity = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})


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
		
	ax_K_scatter_biggest_affinity.scatter(np.exp(final_biggest_affinity), np.array(final_biggest)/C, color = colors_kappa[i_kappa], marker = 'o', alpha = .6, label = r'$%.1f$'%kappa)
	#ax_K_scatter_biggest_affinity.hist2d(np.exp(final_biggest_affinity), np.array(final_biggest)/C, bins = [np.logspace(-9, 3, 20), np.linspace(0, 1, 10)], density = True, color = colors_kappa[i_kappa], label = r'$%.1f$'%kappa)
	
	df = pd.DataFrame({'x':np.exp(final_biggest_affinity), 'y':np.array(final_biggest)/C})
	#df = pd.DataFrame({'x':np.array(final_biggest)/C, 'y':np.array(final_biggest)/C})

	sns.kdeplot(data = df, x='x', y='y', ax=ax_K_scatter_biggest_affinity, color = colors_kappa[i_kappa], log_scale = (True, False), cut = 0)#, clip = ((np.min(np.exp(final_biggest_affinity)), np.max(np.exp(final_biggest_affinity))),(0, 1)))
	
my_plot_layout(ax = ax_K_scatter_biggest_affinity, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K_scatter_biggest_affinity.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 0)
ax_K_scatter_biggest_affinity.set_ylim(bottom = -0.05, top = 1.05)
ax_K_scatter_biggest_affinity.set_xlim(left = 1e-9, right = 2e-4)
#ax_K_scatter_biggest_affinity.set_xticks([])
#ax_K_scatter_biggest_affinity.set_yticks([])
#ax_K_scatter_biggest_affinity.set_yticklabels([1, 0.1, 0.01])
fig_K_scatter_biggest_affinity.savefig('../../Figures/1_Dynamics/Ensemble/K_largest_'+energy_model+'.pdf')


