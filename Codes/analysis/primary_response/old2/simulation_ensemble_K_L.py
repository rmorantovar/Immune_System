import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 500
N_enss = [[500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500]]#, 502, 503, 504, 505, 506, 507, 508, 509, 400, 300, 200, 100, 50], [301, 302, 303, 304, 305, 306, 307, 308, 309]]
#N_rs = [[2e8], [2e8, 2e8/2, 2e8/5, 2e8/10], [2e8], [2e8]]
N_rs = [[1e4, 1e5, 5e5, 1e6, 5e6, 1e7, 2e7, 5e7, 1e8, 5e8, 1e9, 1e10]]
linewidths_N_r = [[5, 4, 3, 2, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]]
linestyles_N_r = [['--', '--', '-', '--', '--', '--', '--', '--', '--', '--', '--', '--']]
transparencies_N_r = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

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
Es, dE, Q0, betas = calculate_Q0(0.1, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_K, ax_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_K_bar, ax_K_bar = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

m_bar_theory = np.array([np.sum(1e8*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, 1, lambda_A, N_c, dE)[3]*dE) for t in time])
t_act_theory = time[m_bar_theory>1][0] 
t_act_1 = t_act_theory

# Printing K from Gumbel
Nb = C
#K_array = np.log(1/(1+(Kds/((AA*(Nb))/1))))
K_array = ((Nb/1)/Kds)
p_K = P_min_e_Q0(1e10, Q0, dE)#/K_array**2*(Nb/1)
p_K = p_K/np.sum(np.flip(p_K[:-1])/np.flip(K_array[:-1])*abs(np.diff(np.flip(K_array))))
Kd_r_renorm = Kds[(P_min_e_Q0(1e10, Q0, dE)/Kds)==np.max(P_min_e_Q0(1e10, Q0, dE)/Kds)]
ax_K.hlines(1, 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')

QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_1))/N_A, Es, 1, lambda_A, N_c, dE)[3]
Kd_1_renorm = Kds[(QR/1)==np.max(QR/1)]
ax_K.hlines((C*1)/Kd_1_renorm/(C/Kd_r_renorm), 0, Tf, linewidth = 2, color = 'black', linestyle = 'dashed')

for i_kappa, kappa in enumerate(kappas):
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	custom_lines = []
	custom_labels = []
	max_potency_simulations = dict()
	max_potency_simulations_std = dict()
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
		data, return_data_type = get_data_ensemble_K_L(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path)	

		if(return_data_type):
			K= data[0]
			Counter = data[1]
			K_final = data[2]
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
				Kds_C = np.exp(energies_C)

				K_i = [(np.sum(((clone_sizes_C[:,t]-1)/1)/Kds_C)) for t in np.arange(len(time))]

				if(np.sum(K_i)!=0):
					K += K_i
					K_final.append(K_i[-1])
					Counter+=1
				
			f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_K_L.pkl', 'wb')
			pickle.dump([K, Counter, K_final], f, pickle.HIGHEST_PROTOCOL)	

		K = (K/Counter)
		
		if(kappa!=1.0):

			ax_K.plot(time, K/(C/Kd_r_renorm), color = colors_kappa[i_kappa], alpha = .9, linewidth = linewidths_N_r[i_kappa][i_N_r], linestyle = linestyles_N_r[i_kappa][i_N_r])
			custom_lines.append(Line2D([0], [0], color=colors_kappa[i_kappa], lw=linewidths_N_r[i_kappa][i_N_r], ls = linestyles_N_r[i_kappa][i_N_r]))
			custom_labels.append(r'$%.0f \cdot 10^{%d}$'%(10**(np.log10(N_r)%1), int(np.log10(N_r))))
	
		max_potency_simulations[N_r] = np.mean(np.array(K_final)/(C/Kd_r_renorm))
		print(N_r, max_potency_simulations[N_r])
		max_potency_simulations_std[N_r] = np.sqrt(np.var((np.array(K_final)/(C/Kd_r_renorm))))


K_0 = 6e-1
delta_E_r = -(1/2)*np.log(0.1)# - np.log(np.max(N_rs[1]))*(1/1.8 - 1/2)
ax_K.vlines(8.90, (K_0 - K_0*(1-1/np.exp(delta_E_r))), K_0, color = 'darkgrey', linewidth = 4)
ax_K.hlines(2.e5/(C/Kd_r_renorm), 4.2, 4.2 + 0.5*np.log(10)/2, color = 'darkgrey', linewidth = 4)

# custom_lines.append(Line2D([0], [0], color = 'orange', lw=3))
# custom_labels.append('Elite')

a = .39
b = .74
c = 1.16
ax_K.plot(time, ((1.08*C*np.exp(c*lambda_B*(time-t_act_1+a)**(b)))/(1.08*C+(np.exp(c*lambda_B*(time-t_act_1+a)**(b))-1)))/Kd_r_renorm/(C/Kd_r_renorm), linewidth = 5, color = 'black', linestyle = 'dotted', zorder = -20) 


N_r_theory = np.logspace(np.log10(1e3), 10, 50)

beta_tail = np.ones_like(N_r_theory)
for n, N_r in enumerate(N_r_theory):
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
	#beta_tail[n] = np.exp(np.mean(np.log(betas[betas>beta_r])))
	#beta_tail[n] = np.mean(betas[betas>beta_r])
	beta_tail[n] = (np.max(betas)-beta_r)/2
print(1/beta_tail)
fit = N_r_theory**(1/beta_tail)
fit/=(N_r_theory[N_r_theory<1e4][-1]**(1/beta_tail[N_r_theory<1e4][-1]))
fit*=np.array(list(max_potency_simulations.values()))[0]

ax_K_bar.plot(N_rs[0], (np.array(list(max_potency_simulations.values()))), color = 'indigo', linestyle = '', marker = 'D', linewidth = 3, ms = 10, alpha = 1)
#ax_K_bar.errorbar(x = N_rs[0], y = (np.array(list(max_potency_simulations.values()))), yerr = (np.array(list(max_potency_simulations_std.values()))), ls = 'none', color = 'indigo', alpha = .6)

popt, pcov = curve_fit(my_linear_func, np.log(N_rs[0][0:-4]), np.log(np.array(list(max_potency_simulations.values()))[0:-4]))

ax_K_bar.plot(N_r_theory, fit, color = 'indigo', linestyle = '-', marker = '', linewidth = 3, label = r'$1/\bar \beta^{\rm tail}(L)$', ms = 10, alpha = 1)
ax_K_bar.plot(N_r_theory, np.exp(my_linear_func(np.log(N_r_theory), *popt)), color = 'indigo', linestyle = '--', marker = '', linewidth = 3, label = '$%.2f\pm%.2f$ (fit)'%(popt[1], np.sqrt(pcov[1,1])), ms = 10, alpha = .6)

print('%.2e'%np.min(Kds))
print('%.2e'%Kds[betas[:-1]<(1/popt[1])][0])
print('%.2e'%Kds[betas[:-1]<(beta_r)][0])

print(np.exp(np.mean(np.log(betas[betas>beta_r]))), np.mean((betas[betas>beta_r])), (np.max(betas)-beta_r)/2)

print(betas[betas>beta_r])

my_plot_layout(ax = ax_K, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K.legend(handles = custom_lines, labels = custom_labels, fontsize = 26, title_fontsize = 28, title = r'$N_r$')
ax_K.set_xlim(left = 2, right = Tf)
ax_K.set_ylim(bottom = 1e-3, top = 2e0)
#ax_K.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K.set_yticklabels([1, 0.1, 0.01])
fig_K.savefig('../../Figures/1_Dynamics/Ensemble/L%d/K_L_'%L+energy_model+'.pdf')

my_plot_layout(ax = ax_K_bar, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_K_bar.legend(fontsize = 26, title_fontsize = 28)
#ax_K_bar.set_xlim(left = 1, right = Tf-1)
#ax_K_bar.set_ylim(bottom = 1e-3, top = 2e0)
#ax_K_bar.set_yticks([1, 0.1, 0.01, 0.001])
#ax_K_bar.set_yticklabels([1, 0.1, 0.01])
fig_K_bar.savefig('../../Figures/1_Dynamics/Ensemble/L%d/K_L_bar_'%L+energy_model+'.pdf')




