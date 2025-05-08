import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 500
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

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2, 3, 4]
#kappas = [1, 3]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_red])
color_list = np.array([my_blue2, my_green, my_red, my_gold])

colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

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

#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
print('beta_r = %.2f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

print('--------')
print('Loops...')
t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
#--------------------------Loops--------------------------
fig_m_f, ax_m_f = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_m_bar, ax_m_bar = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
for i_kappa, kappa in enumerate(kappas):
	
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	delta_t_kappa = (E_kappa-E_pr)*kappa/lambda_A
	delta_t_r = (E_r - E_pr)*kappa/lambda_A
	t_kappa = t_prime + delta_t_kappa
	t_r = t_prime + delta_t_r

	if(kappa>1):
		t_act = t_r
	else:
		t_act = t_kappa - (1/lambda_A)*(np.log(N_r*Q0[Es[:-1]<E_kappa][-1]))
	print('Activation time:',t_act)
	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
	#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
	data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%(L)+parameters_path, data_type = 'm')
	
	if(return_data_type):
		m_bar = data[0]
		m_final = data[1]

	else:	
		m_bar = np.zeros_like(time_array)
		m_final = []
		for i_ens in tqdm(np.arange(N_ens)):
			data_i = data.loc[data['i_ens']==i_ens]
			data_active = data_i.loc[data_i['active']==1]
			t_act_data = np.min(data_active['act_time'])
			data_active = data_active.loc[data_active['act_time']<(t_act_data+1.3+0.3*(kappa-1))] # it was 1.0 + 0.1*...
			activation_times = np.array(data_active['act_time'])
			energies  = np.array(data_active['energy'])

			#---------------------------- B cell linages ----------------------
			clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)

			#--------------------------t_C filter-------------------------
			lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
			clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
			#print(clone_sizes_C)
			#print(clone_sizes_C-1)
			unique, counts = np.unique(np.nonzero((clone_sizes_C - 1))[1], return_counts=True)
			#-------Simulations-------
			m_final.append(len(energies_C))
			m_bar[unique]+=np.log(counts)

			ax_m_bar.plot(unique*dT + T0, counts, color = colors_kappa[i_kappa], alpha = .1)

		f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_m.pkl', 'wb')
		pickle.dump([m_bar, m_final], f, pickle.HIGHEST_PROTOCOL)	
	
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time_array]) 
	ax_m_bar.plot(time_array, np.exp(m_bar/(N_ens)), color = colors_kappa[i_kappa], alpha = 1, label = r'$%d$'%(kappa))
	ax_m_bar.plot(time_array, m_bar_theory, color = colors_kappa[i_kappa], alpha = 1, linestyle = 'dashed')
	ax_m_bar.scatter(t_act, 1, color = colors_kappa[i_kappa])
	#ax_m_f.hist(m_final, alpha = .8, color = colors_kappa[i_kappa], bins = np.logspace(0, 3+np.log10(5), 16), label = r'%d'%(kappa))
	ax_m_f.hist(m_final, alpha = .8, color = colors_kappa[i_kappa], bins = np.linspace(1, 2e3, 50), label = r'$%d$'%(kappa))

my_plot_layout(ax = ax_m_bar, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_m_bar.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
ax_m_bar.set_ylim(bottom = 1e-3, top = 1e5)
#ax_m_bar.set_xlim(left = -3, right = 8.5)
#ax_m_bar.set_xticks([])
#ax_m_bar.set_yticks([])
#ax_m_bar.set_yticklabels([1, 0.1, 0.01])
fig_m_bar.savefig('../../Figures/1_Dynamics/Ensemble/L%d/m_bar_'%(L)+energy_model+'.pdf')

my_plot_layout(ax = ax_m_f, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_m_f.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
#ax_m_f.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#ax_m_f.set_ylim(bottom = -2, top = 4.5)
#ax_m_f.set_yticks([1, 0.1, 0.01, 0.001])
#ax_m_f.set_yticklabels([1, 0.1, 0.01])
fig_m_f.savefig('../../Figures/1_Dynamics/Ensemble/L%d/m_f_'%(L)+energy_model+'.pdf')







