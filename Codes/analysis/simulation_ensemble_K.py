import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 200
N_r = 2e8
transparencies_p = [.8, 1, .8, .8, .8]

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
kappas = [1, 2, 3, 4]#, 5]
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
color_list = np.array([my_red, my_blue2, my_green, my_gold, my_brown])
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
print('N_r = %.e'%N_r)
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_K, ax_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_max, ax_K_max = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

max_potency_simulations = dict()
max_potency_simulations_std = dict()
max_potency_theory = dict()

for i_kappa, kappa in enumerate(kappas):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0] 
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)

	if kappa==1:
		t_act_1 = t_act_theory
	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
	#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
	#data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
	data, return_data_type = get_data_ensemble_K(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)

	if(return_data_type):
		K = data[0]
		Counter = data[1]
	else:

		K = np.zeros_like(time)
		K_final = []
		Counter = 0
		for i_ens in tqdm(np.arange(N_ens)):
			data_i = data.loc[data[4]==i_ens]
			data_active = data_i.loc[data_i[1]==1]
			t_act_data = np.min(data_active[3])
			data_active = data_active.loc[data_active[3]<(t_act_data+1.0+0.1*(kappa-1))]
			activation_times = np.array(data_active[3])
			energies  = np.array(data_active[0])

			#---------------------------- B cell linages ----------------------
			clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

			#--------------------------t_C filter-------------------------
			lim_size = 2
			clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
			#-------Simulations-------
			if(len(energies_C)>0):
				Kds_C = np.exp(energies_C)

				K_i = [np.sum(((clone_sizes_C[:,t]-1)/1)/Kds_C) for t in np.arange(len(time))]

				if(np.sum(K_i)!=0):
					#K += K_i
					K += K_i
					K_final.append(K_i[-1])
					Counter+=1

				#if(i_ens%1==0):
				#	ax_K.plot(time, K_i, color = colors_kappa[i_kappa], alpha = .1, linewidth = 1)

		f = open(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/processed_data_K.pkl', 'wb')
		pickle.dump([K, Counter], f, pickle.HIGHEST_PROTOCOL)	

	K = (K/Counter)

	normalization = 1
	if(kappa==1):
		normalization = 1# K[-1]
	
	ax_K.plot(time, (K/normalization), color = colors_kappa[i_kappa], alpha = .9, linewidth = 5, linestyle = '-', label = r'$%d$'%(kappa))
	
	max_potency_simulations[kappa] = (K[-1]/normalization)
	#max_potency_simulations_std[kappa] = np.sqrt(np.var(np.log(np.array(K_final)/normalization)))

	#Nb = np.exp(lambda_B*Tf)*((k_on*N_c)/(lambda_A*N_A))**(lambda_B/lambda_A)*(k_pr/k_on)**(kappa*lambda_B/lambda_A)*Kds**(-kappa*lambda_B/lambda_A)

	if(kappa==1):
		# Printing K from Gumbel
		Nb = C
		#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
		K_array = ((Nb/1)/Kds)
		p_K = P_min_e_Q0(N_r, Q0, dE)#/K_array**2*(Nb/1)
		p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
		QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_1))/N_A, Es, kappa, lambda_A, N_c, dE)[3]
		Kd_1_renorm = Kds[(QR/1)==np.max(QR/1)]
		ax_K.hlines((C*1)/Kd_1_renorm, 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')
		Kd_r_renorm = Kds[(P_min_e_Q0(N_r, Q0, dE)/Kds)==np.max(P_min_e_Q0(N_r, Q0, dE)/Kds)]
		ax_K.hlines(C/Kd_r_renorm, 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')

print('--------')
# kappas_theory = np.linspace(1, 5.5, 30)
# for kappa in tqdm(kappas_theory):
# 	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
# 	t_act_theory = time[m_bar_theory>=1][0] 
# 	QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_theory))/N_A, Es, kappa, lambda_A, N_c, dE)[3]
# 	numerator = np.sum(dE[:]*QR[:]*np.exp(-lambda_B*kappa/lambda_A*Es[:-1][:] - Es[:-1][:]))
# 	denominator = np.sum(dE[:]*QR[:]*np.exp(-lambda_B*kappa/lambda_A*Es[:-1][:]))
# 	#Es[:-1]>E_r
# 	K = (C*(numerator/denominator))
# 	#print(K)
# 	if(kappa == 1):
# 		normalization_theory = 1
# 	max_potency_theory[kappa] = K - normalization_theory

# ax_K_max.plot(kappas_theory, np.array(list(max_potency_theory.values())), color = my_purple, linestyle = '-', marker = '', linewidth = 3, label = 'theory')

ax_K.plot(time, ((1.08*C*np.exp(1.02*lambda_B*(time-t_act_1+.45)**(0.695)))/(1.08*C+(np.exp(1.02*lambda_B*(time-t_act_1+.45)**(0.695))-1)))/Kd_r_renorm, linewidth = 5, color = 'black', linestyle = 'dotted', zorder = -20)

ax_K_max.plot(kappas, np.array(list(max_potency_simulations.values())), color = my_purple2, linestyle = '', marker = 'D', linewidth = 3, label = 'simulations')
#ax_K_max.errorbar(x=kappas, y=np.array(list(max_potency_simulations.values())), yerr = np.array(list(max_potency_simulations_std.values())), ls = 'none')

print(Kds[betas[:-1]>1][-1])
t_growth = (1/lambda_B)*np.log(C/50)
#ax_K.plot(t_growth+t_prime+kappas_theory/lambda_A*(E_r - E_pr), max_potency_theory.values(), color = 'grey', linewidth = 2)
#ax_K.hlines([C/Kds[betas[:-1]>1][-1], C/Kd_r], 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')


my_plot_layout(ax = ax_K, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K.legend(fontsize = 28, title_fontsize = 30, title = r'$p$', loc = 4)
ax_K.set_xlim(left = 4.5, right = Tf)
ax_K.set_ylim(bottom = 1e9, top = 2e12)
#ax_K.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K.set_yticklabels([1, 0.1, 0.01])
fig_K.savefig('../../Figures/1_Dynamics/Ensemble/K_'+energy_model+'.pdf')

my_plot_layout(ax = ax_K_max, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K_max.legend(fontsize = 28, title_fontsize = 30)
ax_K_max.set_xlim(left = 0.8, right = 4.5)
#ax_K_max.set_ylim(bottom = 20, top = 32)
#ax_K_max.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K_max.set_yticklabels([1, 0.1, 0.01])
fig_K_max.savefig('../../Figures/1_Dynamics/Ensemble/K_max_'+energy_model+'.pdf')



