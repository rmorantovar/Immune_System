import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_r = 1e5
transparencies_p = [.8, 1, .8, .8, .8]
N_ens = 200
T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
k_pr = 120
k_pr = k_pr*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
b0 = 1e5
#N_c = 1e5
#E_ms = -27.63
E_ms = -25
C = 1e4
p=2
times = np.linspace(0, Tf, 1000)

AA = 1
time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18

L=len(antigen)
#print('--------')
#print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path)
#print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.08, 'right':.96, 'bottom':.08, 'top': 0.96})

for N0 in [1e2, 1e3, 1e4]:

	N_r = 1e11/N0
	N_c = N0*b0

	#--------------------------Repertoire properties--------------------------
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
	#-----------------Loading data----------------------------
	lineages_id_sorted = []
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, p, 1e8, linear, N_ens)+energy_model
	data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'N0')
	for i_ens in tqdm(np.arange(N_ens)):
		data_active = data.loc[data['i_ens']==i_ens]
		energies  = np.sort(data_active['energy'])[:100]
		ks = np.exp(energies)*k_on
		lineages = np.arange(len(ks))
		sampled_times = []
		lineages_id = []
		for l in lineages:
			F1 = 1-np.exp(-(b0*1*k_on)/(lambda_A*N_A*(1+ks[l]/k_pr)**p)*(np.exp(lambda_A*times)-1))
			for i in range(int(N0)):
				r1 = np.random.random()
				t1 = times[F1<r1][-1]
				sampled_times.append(t1)
				lineages_id.append(l)

		lineages_id_sorted = np.concatenate((lineages_id_sorted, np.array(lineages_id)[np.argsort(sampled_times)][:100]))

	data_times = np.histogram(lineages_id_sorted, bins = np.arange(101))
	ax.plot(data_times[1][:-1], data_times[0]/N_ens)

fig.savefig('../../Figures/1_Dynamics/times/times_sim_p-%d.pdf'%(p))
