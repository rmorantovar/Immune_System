import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 200
N_enss = [[500], [500], [500, N_ens, N_ens], [500]]#, 502, 503, 504, 505, 506, 507, 508, 509, 400, 300, 200, 100, 50], [301, 302, 303, 304, 305, 306, 307, 308, 309]]
N_rs = [[1e8], [1e8], [1e8, 1e8/10, 1e8/100], [1e8]]
#N_rs = [[1e8, 1e8/2, 1e8/5, 1e8/10, 1e8/20, 1e8/50, 1e8/100]]
linewidths_N_r = [[5], [5], [5, 3, 2], [5], [5]]
ms_N_r = [[12], [12], [12, 10, 8], [12]]
linestyles_N_r = [['-'], ['-'], ['-', '--', '--'], ['-'], ['-']]
transparencies_N_r = [[1], [1], [1, 1, 1], [1], [1]]

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

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0


kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2, 3.0, 4]
#kappas = [3]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_green, my_red, my_gold])
#color_list = np.array([my_green, my_blue2, my_gold])

colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])


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

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
#fig_K, ax_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K, ax_K = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

m_bar_theory = np.array([np.sum(1e8*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, 1, lambda_A, N_c, dE)[3]*dE) for t in time_array])
t_act_theory = time_array[m_bar_theory>1][0] 
t_act_1 = 1.5075376884422111

for i_kappa, kappa in enumerate(kappas):
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_p_properties(betas, Q0, Es, dE, kappa)
	custom_lines = []
	custom_labels = []
	t50s = dict()
	Z50s = dict()
	for i_N_r, N_r in enumerate(N_rs[i_kappa]):
		N_ens = N_enss[i_kappa][i_N_r]
		#--------------------------Repertoire properties--------------------------
		beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
		print('N_r = %.e'%N_r)
		print('beta_r = %.1f'%beta_r)

		m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time_array])
		t_act_theory = time_array[m_bar_theory>1][0] 

		#-----------------Loading data----------------------------
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
		#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		#data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
		data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'K_aging')

		if(N_r == 1e8):
			normalization = 1#K[-1]
			# Printing K from Gumbel
			Nb = C
			#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
			K_array = ((Nb/1)/Kds)
			p_K = P_min_e_Q0(N_r, Q0, dE)#/K_array**2*(Nb/1)
			p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
			Kd_r_renorm = Kds[(P_min_e_Q0(N_r, Q0, dE)/Kds)==np.max(P_min_e_Q0(N_r, Q0, dE)/Kds)]
			#ax_K.hlines(1, 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')
			
			QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_1))/N_A, Es, 1, lambda_A, N_c, dE)[3]
			Kd_1_renorm = Kds[(QR/1)==np.max(QR/1)]
			#ax_K.hlines((C*1)/Kd_1_renorm/(C/Kd_r_renorm), 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')

		if(return_data_type):
			K = data[0]
			Counter = data[1]
			K_is = data[2]
		else:

			K = np.zeros_like(time_array)
			K_final = []
			Counter = 0
			K_is = np.array([], dtype = object)

			for i_ens in tqdm(np.arange(N_ens)):
				data_active = data.loc[data['i_ens']==i_ens]
				#data_active = data_i.loc[data_i['active']==1]
				t_act_data = np.min(data_active['act_time'])
				data_active = data_active.loc[data_active['act_time']<(t_act_data+1.2+0.3*(kappa-1))] # it was 1.0 + 0.1*...
				activation_times = np.array(data_active['act_time'])
				energies  = np.array(data_active['energy'])

				#---------------------------- B cell linages ----------------------
				clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)

				#--------------------------t_C filter-------------------------
				lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
				clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
				#-------Simulations-------
				Kds_C = np.exp(energies_C)

				K_i = [(np.sum(((clone_sizes_C[:,t]-1)/1)/Kds_C)) for t in np.arange(len(time_array))]

				if(np.sum(K_i)!=0):
					#K += K_i
					K += K_i
					K_final.append(K_i[-1])
					Counter+=1
				
				if(((K_i[-1]/(C/Kd_r_renorm))>5e0) and (N_r == 1e8)) :
					ax_K.plot(time_array, K_i/(C/Kd_r_renorm), color = 'orange', linewidth = 3, alpha = .9)
					K_is = np.append(K_is, K_i)

			f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_K_aging.pkl', 'wb')
			pickle.dump([K, Counter, K_is], f, pickle.HIGHEST_PROTOCOL)	

		K = (K/Counter)
		
		normalization = 1
		
		t50 = time_array[K>(K[-1]/2)][0]
		Z50 = K[K>(K[-1]/2)][0]
		t50s[N_r] = t50
		Z50s[N_r] = Z50
		print('changes:', t50s[N_r] - t50s[1e8], Z50s[N_r]/Z50s[1e8])
		if(N_r==1e8):
			ax_K.plot(t50, Z50/(C/Kd_r_renorm), color = colors_kappa[i_kappa], ls = '', lw = 3, marker = 'o', ms = ms_N_r[i_kappa][i_N_r], label = r'$%.1f$'%kappa)
		else:
			ax_K.plot(t50, Z50/(C/Kd_r_renorm), color = colors_kappa[i_kappa], ls = '', lw = 3, marker = 'o', ms = ms_N_r[i_kappa][i_N_r])

		if(kappa==3.0):

			ax_K.plot(time_array, K/(C/Kd_r_renorm), color = colors_kappa[i_kappa], alpha = .9, linewidth = linewidths_N_r[i_kappa][i_N_r], linestyle = linestyles_N_r[i_kappa][i_N_r])
			#custom_lines.append(Line2D([0], [0], color=colors_kappa[i_kappa], lw=linewidths_N_r[i_kappa][i_N_r], ls = linestyles_N_r[i_kappa][i_N_r]))
			#custom_labels.append(r'$%.0f \cdot 10^{%d}$'%(10**(np.log10(N_r)%1), int(np.log10(N_r))))

			for j in range(int(len(K_is)/len(time_array))):
				ax_K.plot(time_array, K_is[j*len(time_array):(j+1)*len(time_array)]/(C/Kd_r_renorm), color = colors_kappa[i_kappa], linewidth = 3, alpha = .9, ls = 'dashdot')

		#Nb = np.exp(lambda_B*Tf)*((k_on*N_c)/(lambda_A*N_A))**(lambda_B/lambda_A)*(k_pr/k_on)**(kappa*lambda_B/lambda_A)*Kds**(-kappa*lambda_B/lambda_A)



K_0 = 4.5e-1
delta_E_r = -(1/3.8)*np.log(0.1)# - np.log(np.max(N_rs[1]))*(1/1.8 - 1/2)
#ax_K.vlines(8.90, (K_0 - K_0*(1-1/np.exp(delta_E_r))), K_0, color = 'darkgrey', linewidth = 4)
ax_K.hlines(2.e5/(C/Kd_r_renorm), 4.2, 4.2 + 0.5*np.log(10)/2, color = 'darkgrey', linewidth = 4)

#custom_lines.append(Line2D([0], [0], color = 'orange', lw=3))
#custom_labels.append('Elite')

#Kd_r_renorm = Kds[(P_min_e_Q0(np.max(N_rs[0]), Q0, dE)/Kds)==np.max(P_min_e_Q0(np.max(N_rs[0]), Q0, dE)/Kds)]
#Kd_r_renorm = Kds[(P_min_e_Q0(np.min(N_rs[1]), Q0, dE)/Kds)==np.max(P_min_e_Q0(np.min(N_rs[1]), Q0, dE)/Kds)] #for the aged repertoire
def my_logistic_func(x, a, b, c):

	return ((1.08*C*np.exp(c*lambda_B*(x-t_act_1+a)**(b)))/(1.08*C+(np.exp(c*lambda_B*(x-t_act_1+a)**(b))-1)))/Kd_r_renorm/(C/Kd_r_renorm)

popt = [-1.17961706, 0.42941023, 2.00880388]

ax_K.plot(time_array, my_logistic_func(time_array, *popt), color = 'black', ls = '-', lw = 2, marker = '', zorder = -20)
#ax_K.plot(time_array, ((C*np.exp(lambda_B*(time_array-t_act_1)))/(C+(np.exp(lambda_B*(time_array-t_act_1))-1)))/Kd_r_renorm, linewidth = 4, color = 'black', linestyle = 'dotted')


my_plot_layout(ax = ax_K, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K.legend(handles = custom_lines, labels = custom_labels, fontsize = 26, title_fontsize = 28, title = r'$N_r$')
ax_K.legend(fontsize = 22, title_fontsize = 24, title = r'$p$')
ax_K.set_xlim(left = 2, right = Tf-1)
ax_K.set_ylim(bottom = 1e-3, top = 10)
#ax_K.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K.set_yticklabels([1, 0.1, 0.01])
fig_K.savefig('../../Figures/1_Dynamics/Ensemble/L%d/K_aging_'%L+energy_model+'.pdf')




