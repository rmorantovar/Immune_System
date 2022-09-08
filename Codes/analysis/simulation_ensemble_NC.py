import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 10
N_r = 1e8
T0 = 3
Tf = 8
Tf_sim = 7
#Tf = 10
dT = 0.01
lambda_A = 6
k_pr = 1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

ns = [2.2, 2.0, 1.8, 1.5]#, 1]
ns = [1.4, 1.8, 2.2]
ns = [3, 2, 1]

transparency_n = [1]

colors_theta = ['lightskyblue', 'tab:cyan','tab:green', 'tab:orange']
colors_theta = ['tab:cyan','tab:green', 'tab:orange']
colors_R = [['deepskyblue', 'lightskyblue', 'lightskyblue'], ['tab:purple', 'tab:cyan', 'tab:cyan'], ['tab:blue', 'tab:green', 'tab:green'], ['tab:red', 'tab:orange', 'tab:orange']]
colors_R = [['tab:purple', 'tab:cyan', 'tab:cyan'], ['tab:blue', 'tab:green', 'tab:green'], ['tab:red', 'tab:orange', 'tab:orange']]

lambda_B = lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
#N_c = 1e5
E_ms = -27.63
C = 3e4

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
PWM_data = get_motif(antigen, energy_model, Text_files_path)
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
for i_n, n in enumerate(ns):

	print('--------')
	print('n = %.2f...'%n)
	beta_n, E_n, Kd_n = get_n_properties(betas, Q0, Es, dE, n)

	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, n, linear, N_ens)+energy_model
	data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)

	NC = np.zeros_like(time)
	for i_ens in tqdm(np.arange(N_ens)):
		data_i = data.loc[data[4]==i_ens]
		data_active = data_i.loc[data_i[1]==1]
		t_act_data = np.min(data_active[3])
		data_active = data_active.loc[data_active[3]<(t_act_data+1.3)]
		activation_times = np.array(data_active[3])
		energies  = np.array(data_active[0])

		#---------------------------- B cell linages ----------------------
		clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

		#--------------------------t_C filter-------------------------
		lim_size = 2
		clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)

		#-------Simulations-------
		Kds_C = np.exp(energies_C)
		NC_i = np.log(1-np.array([np.product(1-1/(1+(Kds_C/((1e12*(clone_sizes_C[:,t]-1))/N_A)))) for t in np.arange(len(time))]))
		NC += NC_i
		if(i_ens%1==0):
			ax_NC.plot(time, NC_i, color = colors_theta[i_n], alpha = .1, linewidth = 1)

	NC = NC/N_ens		
	ax_NC.plot(time, NC, color = colors_theta[i_n], alpha = transparency_n[0], label = r'$%d$'%n, linewidth = 5)

my_plot_layout(ax = ax_NC, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_NC.legend(fontsize = 32, title_fontsize = 34, title = r'$n$')
#ax_NC.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
ax_NC.set_ylim(bottom = -10)
#ax_NC.set_yticks([1, 0.1, 0.01, 0.001])
#ax_NC.set_yticklabels([1, 0.1, 0.01])
fig_NC.savefig('../../Figures/1_Dynamics/Ensemble/NC_'+energy_model+'.pdf')







