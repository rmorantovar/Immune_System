import sys
sys.path.append('../../lib/')
from funcs import*
from matplotlib.colors import LogNorm
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.signal import medfilt

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=80, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = -1)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**7, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--new', type=int, default=0, help="run Z values again.")
	parser.add_argument('--exp', type=int, default=1, help="experiment.")
	args = parser.parse_args()

	# Parameters -----------------------------------------------------

	N_ant = args.N_ant
	N_ens = args.N_ens
	N_inf = args.N_inf
	N_evo = args.N_evo
	N_epi = args.N_epi
	L0 = args.L0
	l = args.l
	new = args.new

	if N_evo == -1:
		N_evo = 'R'

	E_lim = -11.  # Threshold for the sum of entries
	t_lim = 8.  # Threshold for the sum of entries
	chunk_size = 1e6  # Size of each chunk
	p = 4
	k_step = 1/(60*2) #s^-1
	k_step = k_step*3600 # hour^-1
	k_step = k_step*24 #days^-1
	k_on = 1e6*24*3600; #(M*days)^-1
	n_jobs = -1

	T = 12.
	lamA = 6.0
	lamB = 3 * np.log(2) #(days)^-1
	lamB = 2.
	dT = 0.002
	C = 2e4
	v = lamA/p
	b0 = 1e5
	t0 = np.log(lamA/(b0*k_on/N_A))/lamA
	Kstep = k_step/k_on
	ab = 1e-9
	time_array = np.linspace(0, T, int((T-0)/dT))#[::100]

	cmap = plt.get_cmap('cividis_r')
	# Create array of n evenly spaced values between 0 and 1
	color_vals = np.linspace(0, 1, 19)
	my_affinity_colors = [cmap(val) for val in color_vals] 
	color_vals2 = np.linspace(0, 1, 7)
	cmap2 = plt.get_cmap('managua')
	my_pmem_colors = [cmap2(val) for val in color_vals2] 

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	DDEs = [0]
	pmems = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
	exps = [1, 2, 3, 4, 5, 6, 7]
	# pmems = [1, 1.5, 2.0, 2.5, 4.0]
	# exps = [1, 2, 3, 4, 7]
	# figZ, axZ = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
	figNaive, axNaive = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
	figMemory, axMemory = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
	figtheta, axtheta = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})

	max_rank = 5

	theta_naive = []
	theta_naive2 = []

	for i_p, p in enumerate(tqdm([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])):
		# print(i_exp)
		N_ens = args.N_ens
		project = 'memory_response'
		subproject = 'multi-epitope'
		subproject = 'Z_dynamics'
		experiment = 0
		root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
		pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{int(k_step)}_lamA-{lamA}_lamB-{lamB}"
		pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}_N_evo-{N_evo}"
		antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv", converters={"antigen": literal_eval})
		
		dataNaive = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(1) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})

		E_range = (-17, -13.5)

		U_eq = []
		SlogZeq = []

		U_p = []
		SlogZp = []

		U_p2 = []
		SlogZp2 = []

		for i in range(N_ens):
			# print(i)
			#----------NAIVE----------
			dataNaive_i = dataNaive[dataNaive['ens_id']==i].reset_index()

			nsmallest_E_i = dataNaive_i.nsmallest(max_rank, 'E')[['E','N']]
			nsmallest_E_i = nsmallest_E_i.reset_index()
			
			nsmallest_E_i['boltz'] = 1/np.exp(nsmallest_E_i['E'])
			Z_eq_i = np.sum(nsmallest_E_i['boltz'])
			nsmallest_E_i['eq'] = nsmallest_E_i['boltz']/Z_eq_i
			U_eq_i = np.sum(nsmallest_E_i['eq']*nsmallest_E_i['E'])
			S_eq_i = -np.sum(nsmallest_E_i['eq']*np.log(nsmallest_E_i['eq']))

			nsmallest_E_i['z'] = nsmallest_E_i['N']*nsmallest_E_i['boltz']
			Z_p_i = np.sum(nsmallest_E_i['z'])
			nsmallest_E_i['p'] = nsmallest_E_i['z']/Z_p_i
			U_p_i = np.sum(nsmallest_E_i['p']*nsmallest_E_i['E'])
			S_p_i = -np.sum(nsmallest_E_i['p']*np.log(nsmallest_E_i['p']))
			
			U_eq.append(U_eq_i)
			SlogZeq.append(S_eq_i - np.log(Z_eq_i))

			U_p.append(U_p_i)
			SlogZp.append(S_p_i - np.log(Z_p_i))

			#-------------
			nlargest_N_i = dataNaive_i.nlargest(max_rank, 'N')[['E','N']]
			nlargest_N_i = nlargest_N_i.reset_index()
			
			nlargest_N_i['boltz'] = 1/np.exp(nlargest_N_i['E'])
			Z_eq_i = np.sum(nlargest_N_i['boltz'])
			nlargest_N_i['eq'] = nlargest_N_i['boltz']/Z_eq_i
			U_eq_i = np.sum(nlargest_N_i['eq']*nlargest_N_i['E'])
			S_eq_i = -np.sum(nlargest_N_i['eq']*np.log(nlargest_N_i['eq']))

			nlargest_N_i['z'] = nlargest_N_i['N']*nlargest_N_i['boltz']
			Z_p_i = np.sum(nlargest_N_i['z'])
			nlargest_N_i['p'] = nlargest_N_i['z']/Z_p_i
			U_p_i = np.sum(nlargest_N_i['p']*nlargest_N_i['E'])
			S_p_i = -np.sum(nlargest_N_i['p']*np.log(nlargest_N_i['p']))

			U_p2.append(U_p_i)
			SlogZp2.append(S_p_i - np.log(Z_p_i))

		thetap, interceptp = np.polyfit(U_p, SlogZp, deg=1)
		thetap2, interceptp2 = np.polyfit(U_p2, SlogZp2, deg=1)
		thetaeq, intercepteq = np.polyfit(U_eq, SlogZeq, deg=1)

		
		axNaive.scatter(U_p, SlogZp, color = my_pmem_colors[i_p], alpha = .8, marker = 'o', label = r'$%.1f \,;\, %.2f$'%(p, thetap))
		# axNaive.scatter(U_p2, SlogZp2, color = my_pmem_colors[i_p], alpha = .8, marker = '*')
		
		axNaive.plot(np.linspace(*E_range, 100), thetap*np.linspace(*E_range, 100) + interceptp, ls = '--', color = my_pmem_colors[i_p])

		theta_naive.append(thetap)
		theta_naive2.append(thetap2)

	axNaive.scatter(U_eq, SlogZeq, color = 'grey', alpha = .8, label = r'$\textrm{Eq} \,;\, %.2f$'%(thetaeq))
	axNaive.plot(np.linspace(*E_range, 100), thetaeq * np.linspace(*E_range, 100), ls = '--', color = 'grey')
	for n in range(len(exps)):
		axtheta.plot(pmems[n], 1/np.array(theta_naive[n]), color = my_pmem_colors[n], marker = 'D', ms = 12, ls = '', alpha = 0.8)
		axtheta.plot(pmems[n], 1/np.array(theta_naive2[n]), color = my_pmem_colors[n], marker = '*', ms = 15, ls = '', alpha = 0.8)

	for p in [4.0, 3.0, 2.0]:
		theta_memory = []
		for i_exp, exp in enumerate(tqdm(exps)):
			# print(i_exp)
			N_ens = args.N_ens
			pmem = pmems[i_exp]
			project = 'memory_response'
			subproject = 'multi-epitope'
			subproject = 'Z_dynamics'
			experiment = exp
			root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
			pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{int(k_step)}_lamA-{lamA}_lamB-{lamB}"
			pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}_N_evo-{N_evo}"
			antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv", converters={"antigen": literal_eval})
			antigens = antigens_data['antigen']
			
			# subproject = 'Z_dynamics_C5'
			# output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
			# os.makedirs(output_plot, exist_ok=True)
			
			WT = antigens.iloc[0]
			dataNaive = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(1) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})
			#--------------------------Energy Motif--------------------------
			motif = get_motif(WT, energy_model, '../../')*1.2

			#Change values by the minimum
			E_m = -3
			for i in np.arange(l):
				E_m+=np.min(motif[:,i], axis=0)
				motif[:,i]-=np.min(motif[:,i], axis=0)
			# print('Em:', E_m)
			#--------------------------Entropy function--------------------------
			Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
			Kds = np.exp(Es[:-1])
			#--------------------------Repertoire properties--------------------------
			beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
			# if i_exp == 0:
			# 	axtheta.vlines(beta_r, 0.45, 1, color = 'k', ls = '--', lw = 0.5)
			# 	print(beta_r)
			dataNaive = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(1) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})
			dataMemory = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(1) + '/%d'%(0) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})
			E_ms = []
			j=0
			E_range = (-17, -13.5)

			U_q = []
			SlogZq = []

			U_q2 = []
			SlogZq2 = []

			# N_ens = 10
			for i in range(N_ens):
				# print(i)
				#----------NAIVE----------
				dataNaive_i = dataNaive[dataNaive['ens_id']==i].reset_index()

				nsmallest_E_i = dataNaive_i.nsmallest(max_rank, 'E')[['E','N', 'id']]
				nsmallest_E_i = nsmallest_E_i.reset_index()
				
				nsmallest_E_i['boltz'] = 1/np.exp(nsmallest_E_i['E'])
				Z_eq_i = np.sum(nsmallest_E_i['boltz'])
				nsmallest_E_i['eq'] = nsmallest_E_i['boltz']/Z_eq_i
				U_eq_i = np.sum(nsmallest_E_i['eq']*nsmallest_E_i['E'])
				S_eq_i = -np.sum(nsmallest_E_i['eq']*np.log(nsmallest_E_i['eq']))

				nsmallest_E_i['z'] = nsmallest_E_i['N']*nsmallest_E_i['boltz']
				Z_p_i = np.sum(nsmallest_E_i['z'])
				nsmallest_E_i['p'] = nsmallest_E_i['z']/Z_p_i
				U_p_i = np.sum(nsmallest_E_i['p']*nsmallest_E_i['E'])
				S_p_i = -np.sum(nsmallest_E_i['p']*np.log(nsmallest_E_i['p']))
				
				U_eq.append(U_eq_i)
				SlogZeq.append(S_eq_i - np.log(Z_eq_i))

				U_p.append(U_p_i)
				SlogZp.append(S_p_i - np.log(Z_p_i))

				#----------MEMORY----------
				dataMemory_i = dataMemory[dataMemory['ens_id']==i].reset_index()
				# dataMemory_i_t = ensemble_of_expansions_time(dataMemory_i, N_ens, p, time_array, lamB, 1e6, dT)
				# dataMemory_i['N'] = dataMemory_i_t['N_t'].apply(lambda x: x[-1])
				dataMemory_i = dataMemory_i[dataMemory_i['m']==1].reset_index()
				filtered_dataMemory_i = dataMemory_i[dataMemory_i['id'].isin(nsmallest_E_i['id'].to_numpy())]
				# tMemory = dataMemory_i['t'].min()
				nsmallest_E2_i = filtered_dataMemory_i.groupby('id', as_index=False).agg({'N': 'sum', 'E': 'mean'})
				nsmallest_E2_i = nsmallest_E2_i.reset_index()
				# grouped['N_sh'] = np.random.permutation(grouped['N'].values)
				
				nsmallest_E2_i['boltz'] = 1/np.exp(nsmallest_E2_i['E'])
				Z_eq2_i = np.sum(nsmallest_E2_i['boltz'])
				nsmallest_E2_i['eq'] = nsmallest_E2_i['boltz']/Z_eq2_i
				U_eq_i = np.sum(nsmallest_E2_i['eq']*nsmallest_E2_i['E'])
				S_eq_i = -np.sum(nsmallest_E2_i['eq']*np.log(nsmallest_E2_i['eq']))

				nsmallest_E2_i['z'] = nsmallest_E2_i['N']*nsmallest_E2_i['boltz']
				Z_q_i = np.sum(nsmallest_E2_i['z'])
				nsmallest_E2_i['q'] = nsmallest_E2_i['z']/Z_q_i
				U_q_i = np.sum(nsmallest_E2_i['q']*nsmallest_E2_i['E'])
				S_q_i = -np.sum(nsmallest_E2_i['q']*np.log(nsmallest_E2_i['q']))

				U_q.append(U_q_i)
				SlogZq.append(S_q_i - np.log(Z_q_i))

			
			thetaq, interceptq = np.polyfit(U_q, SlogZq, deg=1)
			
			axMemory.scatter(U_q, SlogZq, color = my_pmem_colors[i_exp], alpha = .8, label = r'$%.1f \,;\, %.2f$'%(pmem, thetaq))
			axMemory.plot(np.linspace(*E_range, 100), thetaq*np.linspace(*E_range, 100) + interceptq, ls = '--', color = my_pmem_colors[i_exp])

			theta_memory.append(thetaq)

		axMemory.scatter(U_eq, SlogZeq, color = 'grey', alpha = .8, label = r'$\textrm{Eq} \,;\, %.2f$'%(thetaeq))
		axMemory.plot(np.linspace(*E_range, 100), thetaeq * np.linspace(*E_range, 100), ls = '--', color = 'grey')
		# axNaive.plot(np.linspace(*E_range, 100), thetap*np.linspace(*E_range, 100) + interceptp, ls = '--', color = 'k')
		
		
		if p == 4.0:
			axtheta.plot(pmems, 1/np.array(theta_memory), color = my_pmem_colors[pmems.index(p)], marker = 'o', ms = 12, ls = '', alpha = 0.8)
		else:
			axtheta.plot(pmems, 1/np.array(theta_memory), color = my_pmem_colors[pmems.index(p)], marker = 'o', ms = 12, ls = '', alpha = 0.8)

	
	axtheta.plot(np.linspace(1, 4.0, 10), 1/(1 + lamB * np.linspace(1, 4.0, 10) / lamA), color = 'k', marker = '', ls = '--', label = r'$\theta = 1 + \lambda_Bp/\lambda_A$')

	output_plot = '../../../Figures/'+project+'/'+subproject
	os.makedirs(output_plot, exist_ok=True)

	# my_plot_layout(ax = axZ, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# # axZ.set_xticks([])
	# # axZ.set_yticks([])
	# axZ.set_ylim(bottom = -0.3, top = 7.3)
	# axZ.set_xlim(left = -0.3, right = 7.3)
	# axZ.legend(fontsize = 16, loc = 0, title = r'$p_m$', title_fontsize = 18)
	# figZ.savefig(output_plot + '/dTNbar.pdf')

	my_plot_layout(ax = axNaive, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# axRankCS.set_xticks([])
	# axRankCS.set_yticks([])
	# axRankCS.set_ylim(bottom = 4e-2)
	# axRankCS.set_xlim(left = -0.3, right = 7.3)
	axNaive.legend(fontsize = 16, loc = 0, title = r'$p\, ; \, \theta$', title_fontsize = 18)
	figNaive.savefig(output_plot + '/Naive.pdf')

	my_plot_layout(ax = axMemory, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# axMemory.set_xticks([])
	# axMemory.set_yticks([])
	# axMemory.set_ylim(bottom = 4e-2)
	# axMemory.set_xlim(left = -0.3, right = 7.3)
	axMemory.legend(fontsize = 18, loc = 0, title = r'$p_\textrm{mem}\, ; \, \theta$', title_fontsize = 20)
	figMemory.savefig(output_plot + '/Memory.pdf')

	my_plot_layout(ax = axtheta, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# axtheta.set_xticks([])
	# axtheta.set_yticks([])
	axtheta.set_ylim(top = 1.05)
	# axtheta.set_xlim(left = -0.3, right = 7.3)
	axtheta.legend(fontsize = 20, loc = 0, title_fontsize = 20)
	figtheta.savefig(output_plot + '/theta.pdf')

	
if __name__ == "__main__":
    main()

