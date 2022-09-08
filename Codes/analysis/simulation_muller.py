import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'
#For simulation in C++
# ./EF_response_v2.x -a 6 -b 0.5 -k 1 -t 0 -T 6.5 -E MJ -C 10000 -B 100000000 -s TACNSEYPNTTKCGRWYC -q 2.0
# ./EF_response_v3.x -a 6 -b 0.5 -k 1 -t 0 -T 7.5 -E TCRen -C 10000 -B 100000000 -s EYTACNSEYPNTTKCGRWYCGRYPN -q 2.0 --ensemble -N 10
#--------------- PARAMETERS ---------------------
N_ens = 1
N_r = 1e8
T0 = 3
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.005
lambda_A = 6
k_pr = 1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

ns = [2.2, 2.0, 1.8, 1.5]#, 1]
ns = [1.4, 1.8, 2.2]
ns = [3, 2, 1]

colors_n = ['tab:cyan','tab:green', 'tab:orange', 'orange', 'darkred']
#colors_n = ['tab:orange']
colors_R = [['tab:purple', 'tab:cyan', 'tab:cyan'], ['tab:blue', 'tab:green', 'tab:green'], ['tab:red', 'tab:orange', 'tab:orange'], ['tab:red', 'tab:orange', 'tab:orange']]
#colors_R = [['tab:red', 'tab:orange', 'tab:orange']]

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
for i_n, n in enumerate(ns):
	print('--------')
	print('n = %.2f...'%n)
	beta_n, E_n, Kd_n = get_n_properties(betas, Q0, Es, dE, n)
	for rep in range(1):
		fig_muller, ax_muller = plt.subplots(figsize=(18,6), gridspec_kw={'left':0.06, 'right':.98, 'bottom':.1, 'top': 0.96}, dpi = 600)

		#-----------------Loading data----------------------------
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, n, linear, N_ens)+energy_model
		data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies%d.txt'%rep, sep = '\t', header=None)

		#-----------------Filtering data----------------------------
		min_e_data = np.min(data[0])
		max_e_data = np.max(data[0])

		data_active = data.loc[data[1]==1]
		print('Activated clones before filter:%d'%len(data_active[0]))
		t_act_data = np.min(data_active[3])
		print('t_act_data: %.2f'%t_act_data)
		data_active = data_active.loc[data_active[3]<(t_act_data+1.3)]
		activation_times = np.array(data_active[3])
		print('Activated clones after filter:%d'%len(activation_times))
		energies  = np.array(data_active[0])
		ar1, ar2 = np.histogram(activation_times, bins = time)
		m_data = np.cumsum(ar1)

		#---------------------------- B cell linages ----------------------
		clone_sizes = get_clones_sizes_C(int(m_data[-1]), time, activation_times, lambda_B, C, dT)
		#-----------------------------Activation time------------------------
		t_act = get_t_act(time, N_r, Q0, k_on, k_pr, lambda_A, Es, dE, n, N_c)
		#-----------------------------m(t)-----------------------------------
		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf_sim)/N_A, Es, n, lambda_A, N_c, dE)
		m_f_expected = np.sum(N_r*QR*dE)
		print('Activated clones expected at Tf_sim:%.d'%m_f_expected)

		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_data+1.3))/N_A, Es, n, lambda_A, N_c, dE)
		m_f_expected = np.sum(N_r*QR*dE)
		print('Activated clones expected at t_act_data + 1.3:%.d'%m_f_expected)

		#---------------------- Total Pop size ----------------------
		total_pop = np.sum(clone_sizes, axis = 0)
		total_pop_active = total_pop - total_pop[0] + 1
		t_C = time[total_pop_active<(C-1)][-1] # Calculate time for reaching carrying capacity

		#--------------------------t_C filter-------------------------
		print('Applying filter...')
		lim_size = 10
		clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
		print('Activated clones:', np.shape(clone_sizes_C))

		total_pop_active = np.sum(clone_sizes_C, axis = 0) #C
		bcell_freqs = clone_sizes_C/total_pop_active
		bcell_freqs = clone_sizes_C/(np.sum(clone_sizes_C[:,-1]))
		entropy = -np.sum(bcell_freqs*np.log(bcell_freqs), axis = 0)

		#------------------------- Stackplots -------------------------
		greys = plt.cm.get_cmap('Greys', 50)
		min_bell_freq = np.min(bcell_freqs[:,-1])
		
		for c in np.invert(range(len(clone_sizes_C[:,0]))):
			if bcell_freqs[c, -1]>(0.05):
				ax_muller.stackplot(time, [(bcell_freqs[c, -1] - bcell_freqs[c, :])/2 + np.ones_like(bcell_freqs[0, :])*np.sum(bcell_freqs[:c, -1]), bcell_freqs[c, :], (bcell_freqs[c, -1] - bcell_freqs[c, :])/2], colors = ['white', colors_n[i_n], 'white']);
				ax_muller.scatter(activation_times_C[c], (bcell_freqs[c, -1] - bcell_freqs[c, 0])/2 + np.sum(bcell_freqs[:c, -1]), marker = 'D', edgecolor='black', linewidth=1, facecolor = colors_n[i_n], s = 40)

			else:
				col = greys(np.random.randint(10, 40))
				ax_muller.stackplot(time, [(bcell_freqs[c, -1] - bcell_freqs[c, :])/2 + np.ones_like(bcell_freqs[0, :])*np.sum(bcell_freqs[:c, -1]), bcell_freqs[c, :], (bcell_freqs[c, -1] - bcell_freqs[c, :])/2], colors = ['white', col, 'white']);


		cumsum_freqs = np.cumsum(bcell_freqs, axis = 0)

		# if(i_n!=4):
		# 	for c in range(int(len(clone_sizes_C[:,0]))):
		# 		ax_muller.plot(time, cumsum_freqs[c, :], linewidth = .00001*n, color = 'black')

			
		my_plot_layout(ax = ax_muller, ticks_labelsize=38, yscale = 'linear')
		ax_muller.set_yticks([])
		#ax_muller.set_xticks(np.arange(Tf))
		ax_muller.set_xticks([])
		ax_muller.set_xlim(T0, Tf-2)
		ax_muller.set_ylim(0, 1)
		fig_muller.savefig('../../Figures/1_Dynamics/Trajectories/Muller/B_cell_clones_n-%.2f_%d_'%(n, rep)+energy_model+'.pdf')
		fig_muller.savefig('../../Figures/1_Dynamics/Trajectories/Muller/B_cell_clones_n-%.2f_%d_'%(n, rep)+energy_model+'.png')
		plt.close(fig_muller)
