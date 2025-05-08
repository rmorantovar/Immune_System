import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------
N_ens = 400
N_enss = [[N_ens], [N_ens], [N_ens], [4000, N_ens, N_ens], [N_ens]]#, 502, 503, 504, 505, 506, 507, 508, 509, 400, 300, 200, 100, 50], [301, 302, 303, 304, 305, 306, 307, 308, 309]]
L_0s = [[1e9], [1e9], [1e9], [1e9, 1e9/10, 1e9/100], [1e9]]
#L_0s = [[1e8, 1e8/2, 1e8/5, 1e8/10, 1e8/20, 1e8/50, 1e8/100]]
linewidths_L_0 = [[5], [5], [5], [5, 2, 2], [5]]
ms_L_0 = [[12], [12], [12], [12, 8, 8], [12]]
linestyles_L_0 = [['-'], ['-'], ['-'], ['-', '-', '-'], ['-']]
transparencies_L_0 = [[1], [1], [1], [1, 1, 1], [1]]

T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*10
#N_c = 1e5
#E_m = -27.63
E_m = -24
E_lims = [[-8.5], [-12], [-12],[-12, -12, -12],[-12]]
t_lims = [[4.0], [5.2], [6.2], [7., 9, 9], [8.2]]

#E_lims = [-13]
C = 1e4
AA = 1


time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0


ps = [2.2, 2.0, 1.8, 1.5]#, 1]
ps = [1.4, 1.8, 2.2]
ps = [1, 2, 3, 4, 5]
#ps = [3]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_green, my_purple2, my_red, my_gold])
#color_list = np.array([my_green, my_blue2, my_gold])

colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])


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
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path + "primary_immune_response/in/")
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.001, 80, 1000000, PWM_data, E_m, L)
Kds = np.exp(Es[:-1])

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_step = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
#fig_K, ax_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K, ax_K = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

m_bar_theory = np.array([np.sum(1e9*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, 1, lambda_A, N_c, dE)[3]*dE) for t in time_array])
t_act_theory = time_array[m_bar_theory>1][0] 
t_act_1 = 1.9095477386934674

for i_p, p in enumerate(ps):
	print('--------')
	print('p = %.2f...'%p)
	beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
	custom_lines = []
	custom_labels = []
	t50s = dict()
	Z50s = dict()
	for i_L_0, L_0 in enumerate(L_0s[i_p]):
		N_ens = N_enss[i_p][i_L_0]
		E_lim = E_lims[i_p][i_L_0]
		t_lim = t_lims[i_p][i_L_0]
		#--------------------------Repertoire properties--------------------------
		beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
		print('L_0 = %.e'%L_0)
		print('beta_r = %.1f'%beta_r, 'beta_r = %.1f'%E_r)

		m_bar_theory = np.array([np.sum(L_0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
		t_act_theory = time_array[m_bar_theory>1][0] 

		#-----------------Loading data----------------------------
		#parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
		#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		#data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
		return_data_type = 0
		data, return_data_type = get_data_b(folder_path = '../../out/primary_immune_response', data_type = 'K_aging_p-%.1f_L0-%d'%(p, L_0))

		if(L_0 == 1e9):
			normalization = 1#K[-1]
			# Printing K from Gumbel
			Nb = C
			#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
			K_array = ((Nb/1)/Kds)
			p_K = P_min_e_Q0(L_0, Q0, dE)#/K_array**2*(Nb/1)
			p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
			Kd_r_renorm = Kds[(P_min_e_Q0(L_0, Q0, dE)/1)==np.max(P_min_e_Q0(L_0, Q0, dE)/1)]
			print(C/Kd_r_renorm)
			#ax_K.hlines(1, 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')
			
			QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t_act_1))/N_A, Es, 1, lambda_A, N_c, dE)[3]
			Kd_1_renorm = Kds[(QR/1)==np.max(QR/1)]
			#ax_K.hlines((C*1)/Kd_1_renorm/(C/Kd_r_renorm), 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')

		if(return_data_type):
			K = data[0]
			counter = data[1]
			K_is = data[2]
		else:
			data = pd.read_csv(Text_files_path + 'primary_immune_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
			K = np.zeros_like(time_array)
			K_final = []
			counter = 0
			K_is = np.array([], dtype = object)
			counter3 = False
			counter4 = False

			for i_ens in tqdm(np.arange(N_ens)):
				data_active = data.loc[data['ens_id']==i_ens]
				#data_active = data_i.loc[data_i['active']==1]
				t_act_data = np.min(data_active['time'])
				data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.2*(p-1))] # it was 1.0 + 0.1*...
				activation_times = np.array(data_active['time'])
				energies  = np.array(data_active['E'])
				if len(energies>0):
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
						counter+=1
					
					if(((K_i[-1]/(C/Kd_r_renorm))>4e-0) and ((K_i[-1]/(C/Kd_r_renorm))<5e-0) and (L_0 == 1e9) and (p==4) and (counter3 == False)) :
						print('3', K_i[-1]/(C/Kd_r_renorm))
						K_is = np.append(K_is, K_i)
						K_i = np.array(K_i)
						ax_K.plot(time_array, K_i/(C/Kd_r_renorm), color = colors_p[i_p], linewidth = 2, alpha = .9, ls = 'dashdot')
						ax_K.plot(time_array[K_i<K_i[-1]/2][-1], K_i[K_i<K_i[-1]/2][-1]/(C/Kd_r_renorm), lw = 3, marker = 'o', ms = 8, color = 'lightcoral')
						counter3 = True
					if(((K_i[-1]/(C/Kd_r_renorm))>10e-0) and (L_0 == 1e9) and (p==4) and (counter4 == False)) :
						print('4', K_i[-1]/(C/Kd_r_renorm))
						K_is = np.append(K_is, K_i)
						K_i = np.array(K_i)
						ax_K.plot(time_array, K_i/(C/Kd_r_renorm), color = colors_p[i_p], linewidth = 2, alpha = .9, ls = 'dashdot')
						ax_K.plot(time_array[K_i<K_i[-1]/2][-1], K_i[K_i<K_i[-1]/2][-1]/(C/Kd_r_renorm), lw = 3, marker = 'o', ms = 8, color = 'lightcoral')
						counter4 = True

			f = open('../../out/primary_immune_response' + '/processed_data_K_aging_p-%.1f_L0-%d.pkl'%(p, L_0), 'wb')
			pickle.dump([K, counter, K_is], f, pickle.HIGHEST_PROTOCOL)	

		K = (K/counter)
		
		normalization = 1
		
		t50 = time_array[K>(K[-1]/2)][0]
		Z50 = K[K>(K[-1]/2)][0]
		t50s[L_0] = t50
		Z50s[L_0] = Z50
		print('changes:', t50s[L_0] - t50s[1e9], Z50s[L_0]/Z50s[1e9])
		if(L_0==1e9):
			ax_K.plot(t50, Z50/(C/Kd_r_renorm), color = colors_p[i_p], ls = '', lw = 3, marker = 'o', ms = ms_L_0[i_p][i_L_0], label = r'$%.d$'%p)
		else:
			ax_K.plot(t50, Z50/(C/Kd_r_renorm), color = 'darkred', ls = '', lw = 3, marker = 'o', ms = ms_L_0[i_p][i_L_0])

		if(p==4.0):
			if L_0==1e9:
				ax_K.plot(time_array, K/(C/Kd_r_renorm), color = colors_p[i_p], alpha = .9, linewidth = linewidths_L_0[i_p][i_L_0], linestyle = linestyles_L_0[i_p][i_L_0])
				for j in range(int(len(K_is)/len(time_array))):
					K_i = K_is[j*len(time_array):(j+1)*len(time_array)]
					ax_K.plot(time_array, K_i/(C/Kd_r_renorm), color = 'lightcoral', linewidth = 2, alpha = 1, ls = '-')
					ax_K.plot(time_array[K_i<np.array(K_i)[-1]/2][-1], K_i[K_i<np.array(K_i)[-1]/2][-1]/(C/Kd_r_renorm), lw = 2, marker = 'o', ms = 8, color = 'lightcoral')
			else:
				ax_K.plot(time_array, K/(C/Kd_r_renorm), color = 'darkred', alpha = .9, linewidth = linewidths_L_0[i_p][i_L_0], linestyle = linestyles_L_0[i_p][i_L_0])
			#custom_lines.append(Line2D([0], [0], color=colors_p[i_p], lw=linewidths_L_0[i_p][i_L_0], ls = linestyles_L_0[i_p][i_L_0]))
			#custom_labels.append(r'$%.0f \cdot 10^{%d}$'%(10**(np.log10(L_0)%1), int(np.log10(L_0))))

			

		#Nb = np.exp(lambda_B*Tf)*((k_on*N_c)/(lambda_A*N_A))**(lambda_B/lambda_A)*(k_step/k_on)**(p*lambda_B/lambda_A)*Kds**(-p*lambda_B/lambda_A)

K_0 = 4.5e-1
delta_E_r = -(1/3.8)*np.log(0.1)# - np.log(np.max(L_0s[1]))*(1/1.8 - 1/2)
#ax_K.vlines(8.90, (K_0 - K_0*(1-1/np.exp(delta_E_r))), K_0, color = 'darkgrey', linewidth = 4)
ax_K.hlines(2.e5/(C/Kd_r_renorm), 4.2, 4.2 + 0.5*np.log(10)/2, color = 'darkgrey', linewidth = 4)

#custom_lines.append(Line2D([0], [0], color = 'orange', lw=3))
#custom_labels.append('Elite')

#Kd_r_renorm = Kds[(P_min_e_Q0(np.max(L_0s[0]), Q0, dE)/Kds)==np.max(P_min_e_Q0(np.max(L_0s[0]), Q0, dE)/Kds)]
#Kd_r_renorm = Kds[(P_min_e_Q0(np.min(L_0s[1]), Q0, dE)/Kds)==np.max(P_min_e_Q0(np.min(L_0s[1]), Q0, dE)/Kds)] #for the aged repertoire
def my_logistic_func(x, a, b, c):

	return ((1.*C*np.exp(c*lambda_B*(x+a)**(b)))/(1.*C+(np.exp(c*lambda_B*(x+a)**(b))-1)))/Kd_r_renorm/(C/Kd_r_renorm)

popt = [-3.76090378, 0.31828948, 2.57462231]


ax_K.plot(time_array, my_logistic_func(time_array, *popt), color = 'black', ls = '-', lw = 2, marker = '', zorder = -20)
#ax_K.plot(time_array, ((C*np.exp(lambda_B*(time_array-t_act_1)))/(C+(np.exp(lambda_B*(time_array-t_act_1))-1)))/Kd_r_renorm, linewidth = 4, color = 'black', linestyle = 'dotted')


my_plot_layout(ax = ax_K, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_K.legend(handles = custom_lines, labels = custom_labels, fontsize = 26, title_fontsize = 28, title = r'$L_0$')
ax_K.legend(fontsize = 22, title_fontsize = 24, title = r'$p$')
ax_K.set_xlim(left = 2, right = Tf-1)
ax_K.set_ylim(bottom = 1e-3, top = 14)
#ax_K.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K.set_yticklabels([1, 0.1, 0.01])
fig_K.savefig('../../../Figures/primary_immune_response/1_Dynamics/CSV/K_aging_'+energy_model+'.pdf')




