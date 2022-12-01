import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 500
N_r = 1e8
transparencies_p = [.8, 1, .8, .8, .8]

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
kappas = [1, 2, 2.5, 3, 4]#, 5]
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
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])
#color_list = np.array([my_green, my_blue2, my_gold])

colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

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
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path)
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.1, 50, 400000, PWM_data, E_ms, L)
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
fig_times, ax_times = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

max_potency_simulations = dict()
max_potency_simulations_std = dict()
max_potency_theory = dict()

t_act_p = []
t_f_p = []

t_act_p2 = []
t_f_p2 = []

for i_kappa, kappa in enumerate(kappas):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0] 
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)

	if kappa==1:
		t_act_1 = t_act_theory
	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
	#data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
	#data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
	data, return_data_type = get_data_ensemble_times(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path)

	if(return_data_type):
		t_acts = data[0]
		t_fs = data[1]
	else:
		t_acts = []
		t_fs = []
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
			if(len(energies_C)>1):
				t_acts.append(np.min(activation_times_C))
				t_fs.append(np.max(activation_times_C))

		f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_times.pkl', 'wb')
		pickle.dump([t_acts, t_fs], f, pickle.HIGHEST_PROTOCOL)	


	t_act_p.append(np.mean(t_acts))
	t_f_p.append(np.mean(t_fs))

	t_act_p2.append(np.std(t_acts))
	t_f_p2.append(np.std(t_fs))

	parts = ax_times.violinplot(t_acts, positions = [kappa])

	for pc in parts['bodies']:
		pc.set_facecolor(colors_kappa[3])
		pc.set_edgecolor('black')
		pc.set_alpha(.6)
	for partname in ('cbars','cmins','cmaxes'):
		vp = parts[partname]
		vp.set_edgecolor('darkred')
		vp.set_linewidth(1)

kappas_theory = np.linspace(1, 4, 50)
t_theory = np.ones_like(kappas_theory)
t_theory2 = np.ones_like(kappas_theory)

for i, kappa in enumerate(kappas_theory):
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_numeric = time[m_bar_theory>1][0]
	t_theory2[i] = t_act_numeric
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	if(kappa<beta_r):
		t_theory[i] = t_prime + kappa/lambda_A*np.log(Kd_kappa/Kd_pr) - 1/lambda_A*np.log(N_r*Q0[Es[:-1]<E_kappa][-1])
	else:
		t_theory[i] = t_prime + kappa/lambda_A*np.log(Kd_r/Kd_pr)

#t_theory = t_theory/t_theory[0]*t_act_p[0]
#t_theory = t_theory/t_theory[-1]*t_act_p[-1]
print(t_prime)

ax_times.plot(kappas, t_act_p, color = colors_kappa[3], alpha = 1, linewidth = 3, linestyle = '', marker = 'D', ms = 15, label = r'$t_{\rm act}$')
ax_times.errorbar(x = kappas, y = t_act_p, yerr = t_act_p2, ls = 'none', color = colors_kappa[3], alpha = .6)
ax_times.plot(kappas_theory, t_theory, color = colors_kappa[3], alpha = 1, linewidth = 4, linestyle = '-', marker = '', ms = 15)
ax_times.plot(kappas_theory, t_theory2, color = colors_kappa[3], alpha = 1, linewidth = 4, linestyle = ':', marker = '', ms = 15)

ax_times.plot(kappas, t_f_p, color = colors_kappa[2], alpha = 1, linewidth = 3, linestyle = '', marker = 'o', ms = 15, label = r'$t_f$')
ax_times.errorbar(x = kappas, y = t_f_p, yerr = t_f_p2, ls = 'none', color = colors_kappa[3], alpha = .6)

ax_times.hlines(t_prime, 1, 4, color = 'gray', ls = '--')

print('t_acts =', t_act_p )
print('t_fs = ', t_f_p)
print('--------')

my_plot_layout(ax = ax_times, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_times.legend(fontsize = 28, title_fontsize = 30, title = r'$p$', loc = 4)
#ax_times.set_xlim(left = 2, right = Tf-1)
#ax_times.set_ylim(bottom = 1e-3, top = 2e0)
#ax_times.set_yticks([1, 0.1, 0.01, 0.001])
#ax_times.set_yticklabels([1, 0.1, 0.01])
fig_times.savefig('../../Figures/1_Dynamics/Ensemble/L%d/times_'%L+energy_model+'.pdf')




