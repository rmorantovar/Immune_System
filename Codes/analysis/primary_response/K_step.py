import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------
N_ens = 400
L_0 = 1e9
transparencies_p = [.8, 1, .8, .8, .8]

T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
p = 4.0
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*10
#N_c = 1e5
#E_m = -27.63
E_m = -24
C = 1e4

AA = 1
time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0


k_steps = np.array([3, 30, 300,  1200])#]) # hour^-1
E_lim = -12
t_lims = [8.5, 7, 5.5,  4.5]#]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_green, my_red, my_purple2])
color_list = np.array([my_purple2, my_red, my_green, my_blue2, my_gold])

colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(k_steps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18

L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path + "primary_immune_response/in/")
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_m, L)
Kds = np.exp(Es[:-1])

#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
print('L_0 = %.e'%L_0)
print('beta_r = %.1f'%beta_r)

print('p = %.2f...'%p)
beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
#fig_K, ax_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K, ax_K = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})


max_potency_simulations = dict()
max_potency_simulations_std = dict()
max_potency_theory = dict()

t50s = dict()
Z50s = dict()
Kd_r_renorm_bool = False


# Printing K from Gumbel
Nb = C
#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
K_array = ((Nb/1)/Kds)
p_K = P_min_e_Q0(L_0, Q0, dE)#/K_array**2*(Nb/1)
p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
Kd_r_renorm = Kds[(P_min_e_Q0(L_0, Q0, dE)/Kds)==np.max(P_min_e_Q0(L_0, Q0, dE)/Kds)]

for i_k_step, k_step in enumerate(k_steps):

	k_step = k_step*24 #days^-1
	t_lim = t_lims[i_k_step]
	m_bar_theory = np.array([np.sum(L_0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
	t_act_theory = time_array[m_bar_theory>1][0] 
	print('--------')
	#--------------------------Proofreading properties--------------------------
	beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
	print('beta_pr = %.2f'%beta_pr)

	if p==1:
		t_act_1 = t_act_theory
		print(t_act_1)
	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
	return_data_type = 0
	data, return_data_type = get_data_b(folder_path = '../../out/primary_immune_response', data_type = 'K_2_k_step-%d'%k_step)

	if(return_data_type):
		K = data[0]
		Counter = data[1]
	else:
		data = pd.read_csv(Text_files_path + 'primary_immune_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
		K = np.zeros_like(time_array)
		K_final = []
		Counter = 0
		for i_ens in tqdm(np.arange(N_ens)):
			data_active = data.loc[data['ens_id']==i_ens]
			#data_active = data_i.loc[data_i['active']==1]
			t_act_data = np.min(data_active['time'])
			data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
			activation_times = np.array(data_active['time'])
			energies  = np.array(data_active['E'])

			#---------------------------- B cell linages ----------------------
			clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)

			#--------------------------t_C filter-------------------------
			lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
			clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
			#-------Simulations-------
			if(len(energies_C)>0):
				Kds_C = np.exp(energies_C)

				K_i = [np.sum(((clone_sizes_C[:,t]-1)/1)/Kds_C) for t in np.arange(len(time_array))]

				if(np.sum(K_i)!=0):
					#K += K_i
					K += K_i
					K_final.append(K_i[-1])
					Counter+=1

				#if(i_ens%1==0):
				#	ax_K.plot(time_array, K_i, color = colors_p[i_k_step], alpha = .1, linewidth = 1)

		f = open('../../out/primary_immune_response/processed_data_K_2_k_step-%d.pkl'%k_step, 'wb')
		pickle.dump([K, Counter], f, pickle.HIGHEST_PROTOCOL)	

	K = (K/Counter)

	normalization = 1

	ax_K.plot(time_array, (K/normalization)/(C/Kd_r_renorm), color = colors_p[i_k_step], alpha = .9, linewidth = 5, linestyle = '-', label = r'$%.2f$'%(k_step/(24*60)))#label = r'$%.2f\cdot 10^{%d}$'%(10**(np.log10(k_step/(24*60))%1), int(np.log10(k_step/(24*60)))))

	t50 = time_array[K>(K[-1]/2)][0]
	Z50 = K[K>(K[-1]/2)][0]
	ax_K.plot(t50, Z50/(C/Kd_r_renorm), color = colors_p[i_k_step], ls = '', lw = '3', marker = 'o', ms = 12)
	t50s[k_step] = t50
	Z50s[k_step] = Z50
	print(Z50s[k_step]/Z50s[k_steps[0]*24])

#ax_K.plot(t50s.values(), np.array(list(Z50s.values()))/(C/Kd_r_renorm), color = 'black', ls = '', lw = '2', marker = 'o')

def my_logistic_func(x, a, b, c):

	return ((1.08*C*np.exp(c*lambda_B*(x-t_act_1+a)**(b)))/(1.08*C+(np.exp(c*lambda_B*(x-t_act_1+a)**(b))-1)))/Kd_r_renorm/(C/Kd_r_renorm)

#popt, pcov = curve_fit(my_logistic_func, list(t50s.values()), np.array(list(Z50s.values()))/(C/Kd_r_renorm), p0 = [.3, .7, 1])

#print(popt)
#ax_K.plot(time_array, my_logistic_func(time_array, *popt), color = 'black', ls = '-', lw = '2', marker = '')

print('--------')

a = .39
b = .74
c = 1.16
#ax_K.plot(time_array, ((1.08*C*np.exp(c*lambda_B*(time_array-t_act_1+a)**(b)))/(1.08*C+(np.exp(c*lambda_B*(time_array-t_act_1+a)**(b))-1)))/Kd_r_renorm/(C/Kd_r_renorm), linewidth = 5, color = 'black', linestyle = 'dotted', zorder = -20)

print(Kds[betas[:-1]>1][-1])
t_growth = (1/lambda_B)*np.log(C/50)
#ax_K.plot(t_growth+t_prime+ps_theory/lambda_A*(E_r - E_pr), max_potency_theory.values(), color = 'grey', linewidth = 2)
#ax_K.hlines([C/Kds[betas[:-1]>1][-1], C/Kd_r], 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')


my_plot_layout(ax = ax_K, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K.legend(fontsize = 22, title_fontsize = 24, title = r'$k_{\rm{step }} [\rm{min}^{-1}]$', loc = 0)
ax_K.set_xlim(left = 2, right = Tf-1)
ax_K.set_ylim(bottom = 1e-3, top = 1e0)
#ax_K.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K.set_yticklabels([1, 0.1, 0.01])
fig_K.savefig('../../../Figures/primary_immune_response/1_Dynamics/CSV/K_step_'+energy_model+'.pdf')




