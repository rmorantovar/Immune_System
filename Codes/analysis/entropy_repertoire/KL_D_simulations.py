import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Immune_system/repertoire_entropy/H/'

#--------------- PARAMETERS ---------------------
N_ens = 100
L0 = 1e9
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
E_m = -24
C = 1e4
AA = 1

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

E_lims = [-5.0, -5.5, -6.0, -6.5, -7, -7.5, -8, -8.5, -9.0]
t_lims = [4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8.0]
ps = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5.0]

# E_lims = [-5, -8]
# t_lims = [4, 7]
# ps = [1, 4]

color_list = np.array([my_blue2, my_green, my_purple2, my_red, my_gold])

color_vals = np.linspace(0, 1, 9)
cmap = plt.get_cmap('Reds')
my_colors = [cmap(val) for val in color_vals] 


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
PWM_data, M, Alphabet = get_motif(antigen, energy_model, '../../in/')
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_m, L)
Kds = np.exp(Es[:-1])
#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
print('L0 = %.e'%L0, 'beta_r = %.1f'%beta_r, 'K_r = %.1e'%Kd_r, 'E_r = %.1f'%E_r)
#--------------------------Proofreading properties--------------------------
beta_step, E_step, Kd_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_step = %.2f'%beta_step)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('Loops...')
#--------------------------Loops--------------------------


D0_lineages_p = {}
DR_lineages_p = {}

D0_cells_p = {}
DR_cells_p = {}

D_theory = {}

D_approx = {}

fig_total, ax_total = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
fig0, ax0 = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
for i_p, p in enumerate(ps):
	fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.98})
	fig2, ax2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.98})

	E_lim = E_lims[i_p]
	t_lim = t_lims[i_p]
	m_bar_theory = np.array([np.sum(L0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
	t_act_theory = time_array[m_bar_theory>1][0] 
	
	print('--------')
	print('p = %.2f...'%p)
	beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
	print('K_p = %.1e...'%Kd_p, 'E_p = %.1e...'%E_p)
	print('t_act_theoty:', t_act_theory)
	if p==1:
		t_act_1 = t_act_theory
		print(t_act_1)
	#-----------------Loading data----------------------------
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
	return_data_type = 0
	data, return_data_type = get_data(folder_path = '../../out/repertoire_entropy', data_type = 'cells_p-%.1f'%p)

	if(return_data_type):
		energies_lineages = data[0]
		energies_cells = data[1]
		L_act_p = data[2]
	else:
		data = pd.read_csv(Text_files_path + 'output_N_ens_%d_N_a_1_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f_seqs_0/activated_population_'%(N_ens, L0, p, k_step, E_lim, t_lim, E_m)+antigen+'.csv')
		energies_lineages =[]
		energies_cells =[]
		L_act_p = 0
		counter = 0
		for i_ens in tqdm(np.arange(N_ens)):
			data_active = data.loc[data['ens_id']==i_ens]
			#data_active = data_i.loc[data_i['active']==1]
			t_act_data = np.min(data_active['time'])
			data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.2*(p-1))] # it was 1.0 + 0.1*...
			activation_times = np.array(data_active['time'])
			energies  = np.array(data_active['E'])

			#---------------------------- B cell linages ----------------------
			clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)

			#--------------------------t_C filter-------------------------
			lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
			clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
			#-------Simulations-------
			if(len(energies_C)>0):
				L_act_i = len(energies_C)
				counter+=1

				energies_lineages+=list(energies_C)

				for j in range(len(energies_C)):
					energies_cells+=([energies_C[j]]*int(clone_sizes_C[j, -1]))

			L_act_p+=L_act_i

		f = open('../../out/repertoire_entropy'+'/processed_data_cells_p-%.1f.pkl'%p, 'wb')
		pickle.dump([energies_lineages, energies_cells, L_act_p/counter], f, pickle.HIGHEST_PROTOCOL)	

	print('L_act:', L_act_p)
	D_approx[p] = np.log10(L0/L_act_p)

	bins = np.linspace(-22, E_lim+1, 30)
	n_coarse = 1
	# bins = Es[::n_coarse]
	data_lineages = np.histogram(energies_lineages, density = False, bins = bins)#, cumulative = True, alpha = 0)
	data_cells = np.histogram(energies_cells, density = False, bins = bins)#, cumulative = True, alpha = 0)
	
	P_l = np.diff(np.cumsum(data_lineages[0]))/np.diff(data_lineages[1][:-1])#/(len(energies_lineages))
	print(np.sum(P_l*np.diff(data_lineages[1][:-1])))
	P_l /= np.sum(P_l*np.diff(data_lineages[1][:-1]))
	print(np.sum(P_l*np.diff(data_lineages[1][:-1])))
	P_c = np.diff(np.cumsum(data_cells[0]))/np.diff(data_cells[1][:-1])#/(len(energies_cells))
	P_c /= np.sum(P_c*np.diff(data_lineages[1][:-1]))
	Q_0 = Q0[::n_coarse][:-1]#*L0#/len(energies_lineages)

	P_l = np.interp(Es[:-2], data_lineages[1][:-2], P_l)
	P_c = np.interp(Es[:-2], data_cells[1][:-2], P_c)

	D_temp = []
	ts = np.linspace(t_lim-3.0, t_lim, 40)
	for i_t, t in enumerate(ts):
		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)
		Q_R = QR[::n_coarse][:-1]#*L0#/len(energies_lineages)
		Q_R = Q_R/np.sum(Q_R*dE[::n_coarse][:-1])
		#D = abs(np.sum((np.cumsum(P_l[P_l!=0]))*np.log((np.cumsum(P_l[P_l!=0]))/(np.cumsum(Q_R[P_l!=0]))) * dE[::n_coarse][:-1][P_l!=0]))
		D = abs(np.sum(((P_l[P_l!=0]))*np.log(((P_l[P_l!=0]))/((Q_R[P_l!=0]))) * dE[::n_coarse][:-1][P_l!=0]))
		D_temp.append(D)

	t_optimal = ts[np.argmin(D_temp)]
	print('t_opt', t_optimal)

	u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t_optimal))/N_A, Es, p, lambda_A, N_c, dE)
	Q_R = QR[::n_coarse][:-1]#*L0#/len(energies_lineages)
	Q_R = Q_R/np.sum(Q_R*dE[::n_coarse][:-1])

	ax.plot(np.exp(Es[::n_coarse][:-2][::1]), Q_0[::1], color = 'grey', ls = '--', label = r'$\Omega_{0}$', lw = 4)

	ax.plot(np.exp(Es[:-2][P_l!=0][::80]), P_l[P_l!=0][::80], label = r'$\textrm{Simulation}$', color = my_blue, alpha = .8, ls = '', marker = 'D')
	# ax.plot(np.exp(Es[:-2][P_l!=0]), P_c[P_c!=0], label = r'$\textrm{cells}$', color = my_red, alpha = .8, ls = '-', marker = '')
	# ax_total.plot(np.exp(Es[:-2][P_l!=0][::80]), P_l[P_l!=0][::80], label = r'$%.1f$'%p, color = my_colors[i_p], alpha = .8, ls = '', marker = 'D')

	ax.plot(np.exp(Es[::n_coarse][:-2][::1]), Q_R[::1], color = 'black', ls = '--', label = r'$\Omega_{\textrm{\cal{B}}}$', lw = 4)
	ax_total.plot(np.exp(Es[::n_coarse][:-2][::1]), Q_R[::1], color = my_colors[i_p], ls = '--', lw = 4, label = r'$%.1f$'%p)
	# ax.plot(np.exp(Es[:-1]), np.exp(beta_r * Es[:-1]) * (Q_R[Es[:-2]<(E_r)][-1]) / (np.exp(beta_r * (E_r))), color = my_blue, ls = '--')
	# ax.plot(np.exp(Es[:-1]), np.exp((beta_r - (lambda_B*p)/(lambda_A*1)) * Es[:-1]) * (Q_R[Es[:-2]<(E_r+2)][-1]) / (np.exp((beta_r - (lambda_B*p)/(lambda_A*1)) * (E_r+2))), color = my_red, ls = ':')
	# ax.plot(Es[:-1], np.exp((-4+beta_r) * Es[:-1]) * (Q0[Es[:-1]<(E_r+1.4)][-1]*1e8/len(energies_lineages)) / (np.exp((-4+beta_r) * (E_r+1.4))), color = my_blue, ls = '--')

	ax2.plot(Es[::n_coarse][:-2], np.log10((P_l/1)/(Q_0)), ls = '--', marker = 'o', color = my_blue, alpha = .8, label = r'$\hat{\Omega}_{\textrm{act}}/\Omega_0$')
	ax2.plot(Es[::n_coarse][:-2], np.log10((P_l/1)/(Q_R)), ls = '--', marker = '^', color = my_blue, alpha = .8, label = r'$\hat{\Omega}_{\textrm{act}}/\Omega_{\textrm{act}}$')
	# ax2.plot(Es[::n_coarse][:-2], np.log10((P_c/1)/(Q_0)), ls = '--', marker = 'o', color = my_red, alpha = .8, label = r'$\hat{\Pi}_{\textrm{act}}/\Omega_0$')
	# ax2.plot(Es[::n_coarse][:-2], np.log10((P_c/1)/(Q_R)), ls = '--', marker = '^', color = my_red, alpha = .8, label = r'$\hat{\Pi}_{\textrm{act}}/\Omega_{\textrm{act}}$')

	D0_lineages_p[p] = np.sum((P_l[P_l!=0])*np.log((P_l[P_l!=0]/1)/(Q_0[P_l!=0])) * dE[::n_coarse][:-1][P_l!=0])# - np.sum((P_l[P_l!=0])*np.log((Q_0[P_l!=0])/(1)) * dE[::n_coarse][:-1][P_l!=0])
	DR_lineages_p[p] = np.sum((P_l[P_l!=0])*np.log10((P_l[P_l!=0]/1)/(Q_R[P_l!=0])) * dE[::n_coarse][:-1][P_l!=0])# - np.sum((P_l[P_l!=0])*np.log((Q_R[P_l!=0])/(1)) * dE[::n_coarse][:-1][P_l!=0])

	D0_cells_p[p] = np.sum((P_c[P_c!=0])*np.log10((P_c[P_c!=0]/1)/(Q_0[P_c!=0])) * dE[::n_coarse][:-1][P_c!=0])# - np.sum((P_c[P_c!=0])*np.log((Q_0[P_c!=0])/(1)) * dE[::n_coarse][:-1][P_c!=0])
	DR_cells_p[p] = np.sum((P_c[P_c!=0])*np.log10((P_c[P_c!=0]/1)/(Q_R[P_c!=0])) * dE[::n_coarse][:-1][P_c!=0])# - np.sum((P_c[P_c!=0])*np.log((Q_R[P_c!=0])/(1)) * dE[::n_coarse][:-1][P_c!=0])

	D_theory[p] = np.sum((Q_R[Q_R!=0])*np.log((Q_R[Q_R!=0])/(Q_0[Q_R!=0])) * dE[::n_coarse][:-1][Q_R!=0])# - np.sum((P_c[P_c!=0])*np.log((Q_R[P_c!=0])/(1)) * dE[::n_coarse][:-1][P_c!=0])
	
	my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 28, x_fontsize=30, y_fontsize=30 )
	ax.set_xlim(left = -22, right = -5+1+8+2)
	ax.set_ylim(bottom = 2e-7, top = 1)
	#ax.set_yticks([1, 0.1, 0.01, 0.001])
	#ax.set_yticklabels([1, 0.1, 0.01])
	fig.savefig('../../../Figures/entropy_repertoire/simulations/Pi_vs_Omega_act_'+energy_model+'_p-%.1f_log.pdf'%(p))

	my_plot_layout(ax = ax, xscale='log', yscale= 'linear', ticks_labelsize= 28, x_fontsize=30, y_fontsize=30 )
	ax.legend(fontsize = 18, title_fontsize = 24, loc = 2)
	ax.set_xlim(left = -22, right = -5+1+8+2)
	ax.set_ylim(bottom = 2e-7, top = np.max(Q_R)+.2)
	#ax.set_yticks([1, 0.1, 0.01, 0.001])
	#ax.set_yticklabels([1, 0.1, 0.01])
	fig.savefig('../../../Figures/entropy_repertoire/simulations/Pi_vs_Omega_act_'+energy_model+'_p-%.1f.pdf'%(p))

	my_plot_layout(ax = ax2, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	ax2.legend(fontsize = 16, title_fontsize = 20, loc = 0)
	ax2.set_xlim(left = -21, right = -5+1)
	#ax2.set_ylim(bottom = 2e-4, top = 1e0)
	#ax2.set_yticks([1, 0.1, 0.01, 0.001])
	#ax2.set_yticklabels([1, 0.1, 0.01])
	fig2.savefig('../../../Figures/entropy_repertoire/simulations/Pi_vs_Omega_act_distance_'+energy_model+'_p-%.1f.pdf'%(p))

ax_total.plot(np.exp(Es[::n_coarse][:-2][::1]), Q_0[::1], color = 'grey', ls = '--', label = r'$\Omega_{0}$', lw = 4, zorder = 0)

# print(D_approx)
print(D0_lineages_p.values())

ax0.plot(ps, [np.exp(i) for i in D0_lineages_p.values()], marker = 'D', color = my_blue, label = r'$\textrm{Simulation}$', lw = 3, ms = 15, alpha = .8)
# ax0.plot(ps, DR_lineages_p.values(), marker = 'D', color = my_red, label = r'$D_{\hat{\Omega}_{\textrm{act}}||\Omega_{\textrm{act}}}$')
# ax0.plot(ps, D0_cells_p.values(), marker = 'D', color = my_green, label = r'$D_{\hat{\Pi}_{\textrm{act}}||\Omega_{0}}$')
# ax0.plot(ps, DR_cells_p.values(), marker = 'D', color = my_purple, label = r'$D_{\hat{\Pi}_{\textrm{act}}||\Omega_{\textrm{act}}}$')
ax0.plot(ps, [np.exp(i) for i in D_theory.values()], marker = 'D', color = 'black', label = r'$\textrm{Theory}$', ls = '--', lw = 3, ms = 15, alpha = .8)
# ax0.plot(ps, D_approx.values(), marker = 'o', color = 'grey', label = r'$\textrm{Approx}$', ls = '--')

my_plot_layout(ax = ax0, xscale='linear', yscale= 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax0.legend(fontsize = 34, title_fontsize = 32, loc = 4)
#ax0.set_xlim(left = -21, right = -7)
#ax0.set_ylim(bottom = 2e-7, top = 3e-2)
#ax0.set_yticks([1, 0.1, 0.01, 0.001])
#ax0.set_yticklabels([1, 0.1, 0.01])
fig0.savefig('../../../Figures/entropy_repertoire/simulations/D_'+energy_model+'.pdf')

my_plot_layout(ax = ax_total, xscale='log', yscale= 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_total.set_xlim(left = 2e-10, right = 8e-1)
ax_total.set_ylim(bottom = 2e-5, top = 2)
#ax_total.set_yticks([1, 0.1, 0.01, 0.001])
#ax_total.set_yticklabels([1, 0.1, 0.01])
# ax_total.legend(fontsize = 24, title = r'$p$', title_fontsize = 26, loc = 0)
fig_total.savefig('../../../Figures/entropy_repertoire/simulations/Omega_p_'+energy_model+'_log.pdf')

