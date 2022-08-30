import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 1
N_r = 1e8
T0 = 3
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.01
lambda_A = 6
k_pr = 0.1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

ns = [2.2, 2.0, 1.8, 1.5]#, 1]
ns = [1.4, 1.8, 2.2]
ns = [2, 1]

colors_theta = ['tab:blue','darkblue', 'olive', 'orange', 'darkred']
colors_theta = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
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
colors = ['tab:blue', 'tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = [0]
linear = 0

# antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
# antigen = 'FMLFMAVFVMTSWYC'
# antigen = 'FTSENAYCGR'
# antigen = 'TACNSEYPNTTK'
# antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'
antigen = 'TACNSEYPNTTKCGRWYC'
L=len(antigen)
#----------------------------------------------------------------

#--------------------------Energy Motif--------------------------
PWM_data = get_motif(antigen, 'MJ2', Text_files_path)
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

#--------------------------Loops--------------------------
for i_n, n in enumerate(ns):
	print('n = %.2f...'%n)
	for rep in range(1):
		fig_clones2, ax_clones2 = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		beta_theta, E_theta, Kd_theta = get_theta_properties(betas, Q0, Es, dE, n)

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

		#---- B cell linages ----
		clone_sizes = get_clones_sizes_C(int(m_data[-1]), time, activation_times, lambda_B, C, dT, )
		#-----------------------------Activation time------------------------
		t_act = get_t_act(time, N_r, Q0, k_on, k_pr, lambda_A, Es, dE, n, N_c)
		#-----------------------------m(t)-----------------------------------
		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf_sim)/N_A, Es, n, lambda_A, N_c, dE)
		m_f_expected = np.sum(N_r*QR*dE)
		print('Activated clones expected:%.d'%m_f_expected)

		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_data+1.3))/N_A, Es, n, lambda_A, N_c, dE)
		m_f_expected = np.sum(N_r*QR*dE)
		print('Activated clones expected:%.d'%m_f_expected)

#----------------------------------------------------------------

#--------------------------Entropy function--------------------------

#----------------------------------------------------------------

#--------------------------Entropy function--------------------------

#----------------------------------------------------------------

