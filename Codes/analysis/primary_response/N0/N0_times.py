import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
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
b0 = 1e5
#N_c = 1e5
#E_ms = -27.63
E_ms = -25
C = 1e4
p=3

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

fig, ax = plt.subplots(3, 1, figsize=(5*1.62,5*3/2), gridspec_kw={'left':0.08, 'right':.96, 'bottom':.08, 'top': 0.96})
fig2, ax2 = plt.subplots(3, 1, figsize=(5*1.62,5*3/2), gridspec_kw={'left':0.08, 'right':.96, 'bottom':.08, 'top': 0.96})
N_ens = 1000
N0 = 1e2

times = np.linspace(0, Tf, 10000)

#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, 1e11/N0)

Tc = 5
ki = Kd_r*k_on

for i_k, k in enumerate(np.array([ki, ki*np.exp(1), ki*np.exp(2)])):

	F1 = 1-np.exp(-(b0*1*k_on)/(lambda_A*N_A*(1+k/k_pr)**p)*(np.exp(lambda_A*times)-1))
	FN0 = 1-np.exp(-(b0*N0*k_on)/(lambda_A*N_A*(1+k/k_pr)**p)*(np.exp(lambda_A*times)-1))
	P1 = ((b0*1*k_on)/(lambda_A*N_A*(1+k/k_pr)**p)*(np.exp(lambda_A*times)-1))*np.exp(-(b0*1*k_on)/(lambda_A*N_A*(1+k/k_pr)**p)*(np.exp(lambda_A*times)-1))
	P1/=np.sum(P1[:-1]*np.diff(times))
	PN0 = ((b0*N0*k_on)/(lambda_A*N_A*(1+k/k_pr)**p)*(np.exp(lambda_A*times)-1))*np.exp(-(b0*N0*k_on)/(lambda_A*N_A*(1+k/k_pr)**p)*(np.exp(lambda_A*times)-1))
	PN0/=np.sum(PN0[:-1]*np.diff(times))

	sampled_times_1 = []
	sampled_times_i_tot = []

	sampled_sizes_1 = []
	sampled_sizes_i_tot = []

	for j in tqdm(range(N_ens)):
		sampled_times_i = []
		sampled_sizes_i = []
		for i in range(int(N0)):
			r1 = np.random.random()
			t1 = times[F1<r1][-1]
			sampled_times_i.append(t1)
			sampled_sizes_i.append(np.exp(lambda_B*(Tf-t1)))
			sampled_times_i_tot.append(t1)

		sampled_times_1.append(np.min(sampled_times_i))
		
		sampled_sizes_1.append(np.exp(lambda_B*(Tf-np.min(sampled_times_i))))
		sampled_sizes_i_tot.append(np.sum(np.array(sampled_sizes_i)[np.array(sampled_times_i)<Tc]))

	ax[i_k].hist(sampled_times_1, bins = 20, color = my_blue, alpha = .5, density = True)
	ax[i_k].hist(sampled_times_i_tot, bins = 20, color = my_blue, alpha = .8, density = True)
	ax2[i_k].hist(np.array(sampled_sizes_1)/np.array(sampled_sizes_i_tot), bins = 20, color = my_blue, alpha = .8, density = True)

	sampled_times_N0 = []
	for j in tqdm(range(N_ens)):
		r1 = np.random.random()
		t1 = times[FN0<r1][-1]
		sampled_times_N0.append(t1/(1))

	ax[i_k].hist(sampled_times_N0, bins = 20, color = my_red, alpha = .8, density = True)


	my_plot_layout(ax = ax[i_k])
	ax[i_k].set_xlim(0, Tf)

	ax[i_k].plot(times, P1, color = my_blue, alpha = .4)
	ax[i_k].plot(times, PN0, color = my_red, alpha = .8)

	print(np.sum(P1[:-1]*np.diff(times)))
fig.savefig('../../Figures/1_Dynamics/times/times_N0-%d.pdf'%N0)
fig2.savefig('../../Figures/1_Dynamics/times/sizes_N0-%d.pdf'%N0)

# plt.plot(times/(24*3600), F)
# plt.show()















