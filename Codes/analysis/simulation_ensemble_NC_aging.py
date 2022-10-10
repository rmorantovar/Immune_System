import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 200
N_rs = [[2e8], [2e8, 2e8/2, 2e8/5, 2e8/20]]
linewidths_N_r = [[5], [5, 4, 3, 2]]
T0 = 3
Tf = 10
Tf_sim = 7.5
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2]
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
color_list = np.array([my_red, my_green])

colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

lambda_B = lambda_A
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

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

fig_NC, ax_NC = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_distribution, ax_NC_distribution = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_NC_distribution2, ax_NC_distribution2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------

for i_kappa, kappa in enumerate(kappas):

	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)

	for i_N_r, N_r in enumerate(N_rs[i_kappa]):
		#--------------------------Repertoire properties--------------------------
		beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
		print('N_r = %.e'%N_r)
		print('beta_r = %.1f'%beta_r)

		#-----------------Loading data----------------------------
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
		#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)

		NC = np.zeros_like(time)
		NC_final = []
		for i_ens in tqdm(np.arange(N_ens)):
			data_i = data.loc[data[4]==i_ens]
			data_active = data_i.loc[data_i[1]==1]
			t_act_data = np.min(data_active[3])
			data_active = data_active.loc[data_active[3]<(t_act_data+1.0)]
			activation_times = np.array(data_active[3])
			energies  = np.array(data_active[0])

			#---------------------------- B cell linages ----------------------
			clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

			#--------------------------t_C filter-------------------------
			lim_size = 2
			clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)

			#-------Simulations-------
			Kds_C = np.exp(energies_C)
			NC_i = np.log(1-np.array([np.product(1-1/(1+(Kds_C/((AA*(clone_sizes_C[:,t]-1))/N_A)))) for t in np.arange(len(time))]))
			NC += NC_i
			NC_final.append(NC_i[-1])
			#if(i_ens%1==0):
			#	ax_NC.plot(time, NC_i, color = colors_kappa[i_kappa], alpha = .1, linewidth = 1)

		NC = NC/N_ens
		if(i_kappa==0):
			normalization = NC[-1]
			print(normalization)
		ax_NC.plot(time, NC-normalization, color = colors_kappa[i_kappa], alpha = 1, label = r'$%.e$'%N_r, linewidth = linewidths_N_r[i_kappa][i_N_r])

		NC_data = np.histogram(np.array(NC_final)-normalization, bins = np.linspace(-5, 7, 40), density = False)
		ax_NC_distribution.plot(NC_data[1][:-1], NC_data[0]/N_ens, color = colors_kappa[i_kappa], linestyle='-', marker = 's', label = r'$%.e$'%N_r, linewidth = linewidths_N_r[i_kappa][i_N_r])

		ax_NC_distribution2.plot(NC_data[1][:-1], np.cumsum(NC_data[0]/N_ens), color = colors_kappa[i_kappa], linestyle='-', marker = '', label = r'$%.e$'%N_r, linewidth = linewidths_N_r[i_kappa][i_N_r])
		#Nb = np.exp(lambda_B*Tf)*((k_on*N_c)/(lambda_A*N_A))**(lambda_B/lambda_A)*(k_pr/k_on)**(kappa*lambda_B/lambda_A)*Kds**(-kappa*lambda_B/lambda_A)

		if(kappa==2):
			Nb = C
			NC_array = np.log(1/(1+(Kds/((AA*(Nb))/N_A))))
			p_NC = P_min_e_Q0(N_r, Q0, dE)*(Nb*AA/N_A)/NC_array**1
			p_NC = p_NC/np.sum(p_NC[:-1]*abs(np.diff(NC_array)))
			ax_NC_distribution.plot(np.flip(NC_array[:-1]-normalization), np.flip(p_NC[:-1]), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = linewidths_N_r[i_kappa][i_N_r], alpha = .8)
			ax_NC_distribution2.plot(np.flip(NC_array[:-1]-normalization), np.cumsum(np.flip(p_NC[:-1])*abs(np.diff(np.flip(NC_array)))), linestyle = '--', marker = '',  color = 'black', ms = 2, linewidth = linewidths_N_r[i_kappa][i_N_r], alpha = .8)

my_plot_layout(ax = ax_NC_distribution, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_NC_distribution.legend(fontsize = 32, title_fontsize = 34, title = r'$p$', loc = 4)
#ax_NC_distribution.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
ax_NC_distribution.set_ylim(bottom = 2e-3, top = 1)
ax_NC_distribution.set_xlim(left = -3, right = 7.5)
#ax_NC_distribution.set_xticks([])
#ax_NC_distribution.set_yticks([])
#ax_NC_distribution.set_yticklabels([1, 0.1, 0.01])
fig_NC_distribution.savefig('../../Figures/1_Dynamics/Ensemble/NC_P_aging_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC_distribution2, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_NC_distribution2.legend(fontsize = 32, title_fontsize = 34, title = r'$N_r$', loc = 4)
#ax_NC_distribution2.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
ax_NC_distribution2.set_ylim(bottom = 1e-20)
ax_NC_distribution2.set_xlim(left = -3, right = 8.5)
#ax_NC_distribution2.set_xticks([])
#ax_NC_distribution2.set_yticks([])
#ax_NC_distribution2.set_yticklabels([1, 0.1, 0.01])
fig_NC_distribution2.savefig('../../Figures/1_Dynamics/Ensemble/NC_F_aging_'+energy_model+'.pdf')

my_plot_layout(ax = ax_NC, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_NC.legend(fontsize = 30, title_fontsize = 34, title = r'$N_r$', loc = 4)
#ax_NC.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
ax_NC.set_ylim(bottom = -2, top = 2.5)
#ax_NC.set_yticks([1, 0.1, 0.01, 0.001])
#ax_NC.set_yticklabels([1, 0.1, 0.01])
fig_NC.savefig('../../Figures/1_Dynamics/Ensemble/NC_aging_'+energy_model+'.pdf')







