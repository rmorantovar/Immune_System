import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 500
N_enss = [[500, 500, 500, 500]]#, 502, 503, 504, 505, 506, 507, 508, 509, 400, 300, 200, 100, 50], [301, 302, 303, 304, 305, 306, 307, 308, 309]]
#N_rs = [[2e8], [2e8, 2e8/2, 2e8/5, 2e8/10], [2e8], [2e8]]
N_rs = [[1e8, 1e8/2, 1e8/5, 1e8/10]]
N_rs = [[1e8, 1e8/20]]
linewidths_N_r = [[5, 4, 3, 2], [5], [5]]
linestyles_N_r = [['-', '--', '--', '--'], ['-'], ['-']]
transparencies_N_r = [[1, 1, 1, 1], [.4], [.4]]

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

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0


kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [3.0]
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
color_list = np.array([my_red, my_blue])
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
fig_K, ax_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
m_bar_theory = np.array([np.sum(1e8*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, 1, lambda_A, N_c, dE)[3]*dE) for t in time])
t_act_theory = time[m_bar_theory>1][0] 
t_act_1 = t_act_theory
for i_kappa, kappa in enumerate(kappas):
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	custom_lines = []
	custom_labels = []
	for i_N_r, N_r in enumerate(N_rs[i_kappa]):
		N_ens = N_enss[i_kappa][i_N_r]
		#--------------------------Repertoire properties--------------------------
		beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
		print('N_r = %.e'%N_r)
		print('beta_r = %.1f'%beta_r)

		m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
		t_act_theory = time[m_bar_theory>1][0] 

		#-----------------Loading data----------------------------
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
		#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
		#data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
		data, return_data_type = get_data_ensemble_K_aging(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path)

		if(N_r == 1e8):
			normalization = 1#K[-1]
			# Printing K from Gumbel
			Nb = C
			#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
			K_array = ((Nb/1)/Kds)
			p_K = P_min_e_Q0(N_r, Q0, dE)#/K_array**2*(Nb/1)
			p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
			Kd_r_renorm = Kds[(P_min_e_Q0(N_r, Q0, dE)/Kds)==np.max(P_min_e_Q0(N_r, Q0, dE)/Kds)]
			ax_K.hlines(1, 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')
			
			QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_1))/N_A, Es, 1, lambda_A, N_c, dE)[3]
			Kd_1_renorm = Kds[(QR/1)==np.max(QR/1)]
			ax_K.hlines((C*1)/Kd_1_renorm/(C/Kd_r_renorm), 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')

		if(return_data_type):
			K = data[0]
			Counter = data[1]
			K_is = data[2]
		else:

			K = np.zeros_like(time)
			K_final = []
			Counter = 0
			K_is = np.array([], dtype = object)

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
				Kds_C = np.exp(energies_C)

				K_i = [(np.sum(((clone_sizes_C[:,t]-1)/1)/Kds_C)) for t in np.arange(len(time))]

				if(np.sum(K_i)!=0):
					#K += K_i
					K += K_i
					K_final.append(K_i[-1])
					Counter+=1
				
				if(((K_i[-1]/(C/Kd_r_renorm))>5e0) and (N_r == 1e8)) :
					ax_K.plot(time, K_i/(C/Kd_r_renorm), color = 'orange', linewidth = 3, alpha = .9)
					K_is = np.append(K_is, K_i)

			f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_K_aging.pkl', 'wb')
			pickle.dump([K, Counter, K_is], f, pickle.HIGHEST_PROTOCOL)	

		K = (K/Counter)
		
		normalization = 1
		
		if(kappa!=1.0):

			ax_K.plot(time, K/(C/Kd_r_renorm), color = colors_kappa[i_kappa], alpha = .9, linewidth = linewidths_N_r[i_kappa][i_N_r], linestyle = linestyles_N_r[i_kappa][i_N_r])
			custom_lines.append(Line2D([0], [0], color=colors_kappa[i_kappa], lw=linewidths_N_r[i_kappa][i_N_r], ls = linestyles_N_r[i_kappa][i_N_r]))
			custom_labels.append(r'$%.0f \cdot 10^{%d}$'%(10**(np.log10(N_r)%1), int(np.log10(N_r))))

			for j in range(int(len(K_is)/len(time))):
				ax_K.plot(time, K_is[j*len(time):(j+1)*len(time)]/(C/Kd_r_renorm), color = 'orange', linewidth = 3, alpha = .9)
	
		#Nb = np.exp(lambda_B*Tf)*((k_on*N_c)/(lambda_A*N_A))**(lambda_B/lambda_A)*(k_pr/k_on)**(kappa*lambda_B/lambda_A)*Kds**(-kappa*lambda_B/lambda_A)

K_0 = 4.5e-1
delta_E_r = -(1/3.8)*np.log(0.1)# - np.log(np.max(N_rs[1]))*(1/1.8 - 1/2)
ax_K.vlines(8.90, (K_0 - K_0*(1-1/np.exp(delta_E_r))), K_0, color = 'darkgrey', linewidth = 4)
ax_K.hlines(2.e5/(C/Kd_r_renorm), 4.2, 4.2 + 0.5*np.log(10)/2, color = 'darkgrey', linewidth = 4)

custom_lines.append(Line2D([0], [0], color = 'orange', lw=3))
custom_labels.append('Elite')

a = .39
b = .74
c = 1.16

#ax_K.plot(time, ((C*np.exp(lambda_B*(time-t_act_1)))/(C+(np.exp(lambda_B*(time-t_act_1))-1)))/Kd_r_renorm, linewidth = 4, color = 'black', linestyle = 'dotted')
Kd_r_renorm = Kds[(P_min_e_Q0(np.max(N_rs[0]), Q0, dE)/Kds)==np.max(P_min_e_Q0(np.max(N_rs[0]), Q0, dE)/Kds)]
ax_K.plot(time, ((1.08*C*np.exp(c*lambda_B*(time-t_act_1+a)**(b)))/(1.08*C+(np.exp(c*lambda_B*(time-t_act_1+a)**(b))-1)))/Kd_r_renorm/(C/Kd_r_renorm), linewidth = 5, color = 'black', linestyle = 'dotted', zorder = -20)
#Kd_r_renorm = Kds[(P_min_e_Q0(np.min(N_rs[1]), Q0, dE)/Kds)==np.max(P_min_e_Q0(np.min(N_rs[1]), Q0, dE)/Kds)] #for the aged repertoire
#ax_K.plot(time, ((1.08*C*np.exp(1.02*lambda_B*(time-t_act_1+.45-0.5*np.log(10)/2)**(0.695)))/(1.08*C+(np.exp(1.02*lambda_B*(time-t_act_1+.45-0.5*np.log(10)/2)**(0.695))-1)))/Kd_r_renorm, linewidth = 5, color = 'black', linestyle = 'dotted', zorder = -20)

my_plot_layout(ax = ax_K, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K.legend(handles = custom_lines, labels = custom_labels, fontsize = 26, title_fontsize = 28, title = r'$N_r$')
ax_K.set_xlim(left = 3, right = Tf-1)
ax_K.set_ylim(bottom = 1e-2, top = 10)
#ax_K.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K.set_yticklabels([1, 0.1, 0.01])
fig_K.savefig('../../Figures/1_Dynamics/Ensemble/L%d/K_aging_'%L+energy_model+'.pdf')




