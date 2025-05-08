import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
transparencies_p = [.8, 1, .8, .8, .8]
N_ens = 500
T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
#k_pr = 120
k_pr = k_pr*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
b0 = 1e5
#N_c = 1e5
#E_ms = -27.63
E_ms = -25
C = 1e4
p=3.0
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

beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)

for N0 in [1e0, 1e1, 1e2, 1e3]:
	fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.08, 'right':.96, 'bottom':.08, 'top': 0.96})
	fig2, ax2 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.08, 'right':.96, 'bottom':.08, 'top': 0.96})
	fig3, ax3 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.08, 'right':.96, 'bottom':.08, 'top': 0.96})

	N_r = 1e8
	N_c = N0*b0

	#--------------------------Repertoire properties--------------------------
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
	beta_act = np.min([beta_r, beta_p])
	#-----------------Loading data----------------------------
	lineages_id_sorted = []
	parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, p, 1e8, linear, N_ens)+energy_model
	data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'N0')
	
	#data_csv = pd.read_csv('data_N0-%d_p-%d.csv'%(N0, p))
	n_first_clones = 50
	final_Nb = np.zeros(n_first_clones)
	counts_final_Nb = np.zeros(n_first_clones)

	final_Nb_l = np.zeros(n_first_clones)
	counts_final_Nb_l = np.zeros(n_first_clones)
	energies_total = dict()
	time_difference = []
	size_difference = []
	max_rank = 50
	for i_ens in tqdm(np.arange(int(N_ens/1))):
		data_active = data.loc[data['i_ens']==i_ens]
		energies  = np.sort(data_active['energy'])[:300]
		#energies = data_csv[str(i_ens)]

		ks = np.exp(energies)*k_on
		lineages = np.arange(len(ks))
		sampled_times = []
		lineages_id = []
		
		energies_total[i_ens] = energies

		for l in lineages:
			F1 = 1-np.exp(-(b0*1*k_on)/(lambda_A*N_A*(1+ks[l]/k_pr)**p)*(np.exp(lambda_A*times)-1))
			for i in range(int(N0)):
				r1 = np.random.random()
				t1 = times[F1<r1][-1]
				sampled_times.append(t1)
				lineages_id.append(l)
		first_activations_times = np.array(sampled_times)[np.argsort(sampled_times)][:300]
		first_activations_ids = np.array(lineages_id)[np.argsort(sampled_times)][:300]
		first_activations_energies = energies[first_activations_ids]

		#---------------------------- B cell linages ----------------------
		clone_sizes = get_clones_sizes_C(len(first_activations_times), time_array, first_activations_times, lambda_B, C, dT)

		#--------------------------t_C filter-------------------------
		lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
		clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, first_activations_times, first_activations_energies, lim_size)
		activations_ids_C = first_activations_ids[filter_C]


		sort_inds = clone_sizes_C[:, -1].argsort()
		clone_sizes_C_sorted = clone_sizes_C[sort_inds, :][-int(n_first_clones):, :]
	
		biggest_clone_i = clone_sizes_C_sorted[-1, -1]
		
		sorted_clones = np.flip(clone_sizes_C_sorted[:, -1])/biggest_clone_i
		max_rank_i = len(sorted_clones)
		if max_rank_i>1:
			for i in range(max_rank_i):
				#final_Nb[i]+= np.log(sorted_clones[i])
				final_Nb[i]+= (sorted_clones[i])
				counts_final_Nb[i] += 1
			if(max_rank_i<max_rank):
				max_rank = max_rank_i


		clone_sizes_lineages = np.zeros(300)
		#sum clone sizes of the same lineage:
		for l in range(300):
			clone_sizes_l_array = clone_sizes_C[:, -1][activations_ids_C==l]
			if(len(clone_sizes_l_array)>1):
				time_difference.append(abs(np.log(np.sort(clone_sizes_l_array)[-1]/np.sort(clone_sizes_l_array)[-2])/lambda_B))
				size_difference.append(np.sort(clone_sizes_l_array)[-1]/np.sum(clone_sizes_l_array))
			clone_sizes_lineages[l] = np.sum(np.concatenate(([0], clone_sizes_l_array)))

		sort_inds_l = clone_sizes_lineages[:].argsort()
		clone_sizes_C_sorted_l = clone_sizes_lineages[sort_inds_l][-int(n_first_clones):]
		
		biggest_clone_i_l = clone_sizes_C_sorted_l[-1]

		sorted_clones_l = np.flip(clone_sizes_C_sorted_l[:])/biggest_clone_i_l
		max_rank_i_l = len(sorted_clones_l)
		if max_rank_i_l>1:
			for i in range(max_rank_i_l):
				#final_Nb[i]+= np.log(sorted_clones[i])
				final_Nb_l[i]+= (sorted_clones_l[i])
				counts_final_Nb_l[i] += 1
			if(max_rank_i_l<max_rank):
				max_rank = max_rank_i_l
				
	df = pd.DataFrame.from_dict(energies_total)
	df.to_csv('data_N0-%d_p-%d.csv'%(N0, p))

	final_Nb = final_Nb/counts_final_Nb
	final_Nb_l = final_Nb_l/counts_final_Nb_l
	ranking = np.arange(1, n_first_clones+1)
	fit0 = ranking**(-lambda_B/(lambda_A))
	fit = ranking**(-p*lambda_B/(lambda_A*beta_act))
	#fit2 = ranking**(-lambda_B/(lambda_A)-1)
	ax.plot(ranking, final_Nb, marker = '*', color = my_blue, lw = 3, ls = '')
	ax.plot(ranking, final_Nb_l, marker = '^', ms = 5, color = my_red, lw = 3, ls = '')
	ax.plot(ranking, fit0, color = my_blue, lw = 2, ls = '--')
	ax.plot(ranking, fit, color = my_blue, lw = 3, ls = '-')
	#ax.plot(ranking, fit2, color = my_red, lw = 3, ls = '-')

	ax2.hist(time_difference, alpha = .8)
	ax2.vlines([1/lambda_A, np.mean(time_difference)], 0, ax2.get_ylim()[1], color = 'black', ls = ['-', '--'])

	ax3.hist(size_difference, alpha = .8, bins = np.logspace(-2, 0, 20))

	ax.set_yscale('log')
	ax.set_xscale('log')
	fig.savefig('../../Figures/1_Dynamics/times/ranking_N0-%d_p-%d.pdf'%(N0, p))

	fig2.savefig('../../Figures/1_Dynamics/times/delta_t12_N0-%d_p-%d.pdf'%(N0, p))

	ax3.set_xscale('log')
	fig3.savefig('../../Figures/1_Dynamics/times/delta_N1all_N0-%d_p-%d.pdf'%(N0, p))






