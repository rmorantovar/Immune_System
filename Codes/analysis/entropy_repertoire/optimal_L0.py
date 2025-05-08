import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")
from collections import defaultdict


Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Immune_system/repertoire_entropy/H/'

#--------------- PARAMETERS ---------------------
L0s = np.logspace(4, 18, 15)
T0 = 0
Tf = 20
Tf_sim = 7
#Tf = 10
dT = 0.02
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
ls = range(15, 25, 1)

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0
replicates = range(1, 10)

# print('--------')
# print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
dumb_antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
PWM_data, M, Alphabet = get_motif(dumb_antigen, energy_model, '../../in/')

fig_D, ax_D = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_p, ax_p = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_t, ax_t = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_step, ax_step = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_K_star, ax_K_star = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})


data = defaultdict(list)
for l in ls:
	print('l=', l)
	for L0 in tqdm(L0s):
		for rep in replicates:
			antigen = ''.join(np.random.choice(Alphabet[:-1], l))
			PWM_data, M, Alphabet = get_motif(antigen, energy_model, '../../in/')
			#Change values by the minimum
			for i in np.arange(l):
			    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)
			#--------------------------Entropy function--------------------------
			Es, dE, Q0, betas = calculate_Q0(0.01, 50, 5000, PWM_data, E_m, l)
			Kds = np.exp(Es[:-1])
			#--------------------------Repertoire properties--------------------------
			beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
			# select p and k_step as the optimal parameters.			
			p = beta_r
			k_step = Kd_r*k_on

			#--------------------------Proofreading properties--------------------------
			beta_step, E_step, Kd_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)

			# Calculate the activation time
			m_bar_theory = np.array([np.sum(L0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
			t_act_theory = time_array[m_bar_theory>1][0] 
			
			beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)

			# n_coarse = 1
			# u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t_act_theory))/N_A, Es, p, lambda_A, N_c, dE)
			# Q_R = QR[::n_coarse][:-1]#*L0#/len(energies_lineages)
			# Q_R = Q_R/np.sum(Q_R*dE[::n_coarse][:-1])
			# D = np.sum((Q_R[Q_R!=0])*np.log10((Q_R[Q_R!=0])/(Q0[:-1][Q_R!=0])) * dE[::n_coarse][:-1][Q_R!=0])# - np.sum((P_c[P_c!=0])*np.log((Q_R[P_c!=0])/(1)) * dE[::n_coarse][:-1][P_c!=0])

			data['L0'].append(L0)
			data['l'].append(l)
			data['p'].append(p)
			data['step'].append(k_step/(24*60))
			data['t'].append(t_act_theory)
			data['K_star'].append(Kd_r)
			#data['antigen'].append(antigen)
			data['rep'].append(rep)

data_df = pd.DataFrame(data)
grouped = data_df.groupby(['L0', 'l'])
agg_funcs = {'p': 'mean', 'step': 'mean', 't': 'mean', 'K_star': 'mean'}#, 'antigen': 'first'}
averages = grouped.agg(agg_funcs)
data_df.to_csv('../../out/repertoire_entropy/data.csv')
averages.to_csv('../../out/repertoire_entropy/data_processed.csv')

