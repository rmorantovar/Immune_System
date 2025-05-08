import sys
sys.path.append('../../lib/')
from funcs import*
from matplotlib.colors import LogNorm

def main():
	# Setting up command-line argument parser
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=40, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 0)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**8, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default=8., help="Threshold for activation time.") # Use 8 for L0>1e6
	parser.add_argument('--E_lim', type=float, default=-6., help="Threshold for the sum of entries.") # Use -6 for L0>1e6
	parser.add_argument('--E_m', type=float, default=-24, help="Threshold for the sum of entries.")
	parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=4, help="# steps.")
	parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
	parser.add_argument('--lamA', type=float, default=6., help="Antigen growth rate.")
	parser.add_argument('--lamB', type=float, default=2., help="Antigen growth rate.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
	parser.add_argument('--random_antigen', type=int, default=0)
	# parser.add_argument('--antigen', type=str, default='TACNSEYPNTTRAKCGRWYR')
	parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR')
	parser.add_argument('--energy_model', type=str, default = 'TCRen')
	parser.add_argument('--seqs', type=int, default = 1)
	parser.add_argument('--one_WT', type=int, default = 1)
	parser.add_argument('--secondary', type=int, default = 0)
	parser.add_argument('--secondary_all', type=int, default = 0)
	parser.add_argument('--pro', type=str, default='epitope_complexity', help="project.")
	parser.add_argument('--subpro', type=str, default='minimal_ID', help="subproject.")
	parser.add_argument('--exp', type=int, default=1, help="experiment.")
	args = parser.parse_args()

	# ------------ PARAMETERS AND INPUTS ------------
	N_ant = args.N_ant
	N_ens = args.N_ens
	N_inf = args.N_inf
	N_evo = args.N_evo
	N_epi = args.N_epi
	L0 = args.L0
	l = args.l
	E_lim = args.E_lim
	t_lim = args.t_lim
	E_m = args.E_m
	if L0>=1e6:
		chunk_size = args.chunk_size
	else:
		chunk_size = args.L0
	p = args.p
	k_step = args.k_step
	lamA = args.lamA
	lamB = args.lamB
	n_jobs = args.n_jobs
	random_antigen = args.random_antigen
	antigen = args.antigen
	energy_model = args.energy_model
	
	seqs = args.seqs
	one_WT = args.one_WT
	secondary = args.secondary
	secondary_all = args.secondary_all

	if N_evo == -1:
		N_evo = 'R'

	dT = 0.002
	C = 1e4
	T = 12
	time_array = np.linspace(0, T, int((T-0)/dT))
	Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)

	project = args.pro
	subproject = args.subpro
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
		
	start_time = time.time()
	print('Starting simulation ...')

	fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	figZ, axZ = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	fig_betas, ax_betas = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	fig_L0, ax_L0 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	results = defaultdict(list)
	for k, L0 in enumerate([10**6, 2*10**6, 10**7, 2*10**7, 10**8]):
		pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
		colors = [my_blue, my_red, my_green, my_cyan, my_gold]
		alpha = 1e-14

		output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
		os.makedirs(output_plot, exist_ok=True)
		N_epi_array = np.arange(3, 6)
		for N_epi in N_epi_array:
			pars_dir_2 = f"/N_ens-{N_ens}_N_epi-{N_epi}"
			antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv", converters={"antigen": literal_eval})
			antigens = antigens_data['antigen']
			# ------------ ------------ ------------
			if one_WT:
				WTs = [antigens.iloc[0]]
			else:
				WTs = list(antigens)
			
			lss = ['-', '--']
			lws = [1, 2]

			for a1, antigen1_ in enumerate(WTs):
				print('primary infection')
				output_dir1 = root_dir + pars_dir_1 + pars_dir_2 + "/%d"%(a1+1)
				output_file1 = os.path.join(output_dir1, 'potency_' + str(a1+1) + '.csv')
				# Calculate motif
				motif = get_motif(antigen1_, energy_model, '../../')*1.2
				E_ms = np.zeros(N_epi)
				E_rs = np.zeros(N_epi)
				beta_rs = np.zeros(N_epi)
				for epi in range(N_epi):
					E_m = -3
					# Normalize motif
					for i in range(l):
						E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
						motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
					# E_m = -23
					E_ms[epi] = E_m
					# Calculate Q0, Es, dE, betas
					Es, dE, Q0, betas = calculate_Q0(0.01, 50, 100000, motif[:, epi*l:(epi+1)*l], 0, l)
					#--------------------------Repertoire properties--------------------------
					beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
					E_rs[epi] = E_r
					beta_rs[epi] = beta_r
					Es = Es - E_r - 18
					# beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
					# print(E_m, E_r, beta_r)

				Kstar_dom = np.exp(-np.mean(E_rs))

				if os.path.isfile(output_file1):
					# read file and adapt format
					data_activation = pd.read_csv(output_file1, converters={"Z_t": literal_eval})
					data_activation['Z_t'] = data_activation['Z_t'].apply(lambda x: np.array(x, dtype=np.float32))
					
					Z_epi = []
					ID_epi = 0

					# fig_1, ax_1 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
					# ax_1.plot(time_array, np.ones_like(time_array)*1, color = 'k', ls = '--', lw = 1, alpha = 1)
					for ens_id in range(N_ens):
						data_ens = data_activation.loc[data_activation['ens_id']==ens_id]
						# for index, row in data_ens.iterrows():
							# ax_1.plot(time_array[::100][row['Z_t']>0], 1/(1+(alpha*row['Z_t'][row['Z_t']>0])**-1), color = colors[row['epi']-1], ls = '-', lw = (row['m']*1)+1, alpha = 0.3)
							# ax_1.plot(time_array[::100][row['Z_t']>0], row['Z_t'][row['Z_t']>0], color = colors[row['epi']-1], ls = '-', lw = (row['m']*1)+1, alpha = 0.3)
						data_ens_mod = data_ens.groupby(['m']).agg({'E':'mean', 
																	't':'mean', 
																	# 'Z_t':lambda x: 1-(1-1/(1+(alpha*x)**-1)).product()}).reset_index() 
																	'Z_t':lambda x: x.sum()}).reset_index() 
						# ax_1.plot(time_array[::100], data_ens_mod['Z_t'][0], color = 'k', ls = '-', lw = (row['m']*1)+1, alpha = .8)
						# ax.scatter(N_epi, data_ens_mod['Z_t'][0][-1], color = colors[N_epi-1], ls = '-', lw = (row['m']*1)+1, alpha = .6)
						Z_epi.append(data_ens_mod['Z_t'][0][-1])
						ID_i = np.array([data_ens[data_ens['epi']==i]['Z_t'].iloc[0][-1] for i in range(1, N_epi+1)])/data_ens_mod['Z_t'][0][-1]
						ID_epi += np.exp(-np.sum(ID_i*np.log(ID_i)))

					ID_epi/=N_ens
					results['L0'].append(L0)
					results['epi_max'].append(N_epi)
					results['betas'].append(np.mean(beta_rs))
					results['g'].append(ID_epi)
					results['Z'].append(np.exp(np.mean(np.log(Z_epi))))
					# my_plot_layout(ax = ax_1, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30, bottom = None, top = None)
					# fig_1.savefig(output_plot+'/Z1_'+str(a1)+f'_L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}'+'_N_epi-'+str(N_epi)+'_all.pdf')

				# for inf in tqdm(range(N_inf-1)): # How many WT recurrent infections?
				# 	antigen2_seq = from_aa_to_i(antigen1, energy_model, '../../')
					
				# 	input_file2 = os.path.join(output_dir1, 'potency_' + str(a1+1) + '.csv')
				# 	output_dir1 = output_dir1 + "/%d"%(a1+1)
				# 	output_file2 = os.path.join(output_dir1, 'potency_' + str(a1+1) + '.csv')

				# 	if os.path.isfile(input_file2):
				# 		if os.path.isfile(output_file2):
				# 			data_activation = pd.read_csv(output_file2, converters={"Z_t": literal_eval})
				# 			data_activation['Z_t'] = data_activation['Z_t'].apply(lambda x: np.array(x, dtype=np.float32))
							
				# 			data_activation.to_csv(output_dir1 + '/potency_' + str(a1+1) + '.csv', index = False)
				# 			fig_inf, ax_inf = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
				# 			for m in [0, 1]:
				# 				ax_inf.plot(time_array, np.ones_like(time_array)*1, color = 'k', ls = '--', lw = 1, alpha = 1)
				# 				for index, row in data_activation.loc[data_activation['m']==m].iterrows():	
				# 					ax_inf.plot(time_array[::100],row['Z_t']/Kstar_dom, color = colors[row['epi']-1], ls = '-', lw = (row['m']*4)+1, alpha = 0.8)
				# 			my_plot_layout(ax = ax_inf, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30, bottom = 1e-1, top = 1e5)
				# 			fig_inf.savefig('../../../Figures/'+project+'/'+subproject+'/'+str(experiment)+'/Z'+str(inf+2)+'_'+str(a1)+'_'+str(0)+'_L0-'+str(L0)+'_N_epi-'+str(N_epi)+'_all.pdf')
		
		# ax.plot(N_epi_array, ID, color = colors[k], label = r'$10^{%d}$'%np.log10(L0))
		# ax.hlines(np.sqrt(np.log(L0)/20), 1, N_epi, color = colors[k], ls = ':')

	results_df = pd.DataFrame(results)

	for k, N_epi in enumerate(N_epi_array):
		results_epi = results_df.loc[results_df['epi_max'] == N_epi]
		ax.plot(results_epi['L0'], results_epi['g']/N_epi, ls = '-', color = colors[k], label = r'$%d$'%N_epi)
		ax_betas.scatter(results_epi['g']/N_epi, np.sqrt(results_epi['betas']), color = colors[k], label = r'$%d$'%N_epi)
		ax_L0.scatter(results_epi['g']/N_epi, np.sqrt(np.log(results_epi['L0'])), marker = '*', color = colors[k], label = r'$%d$'%N_epi)
		axZ.plot(results_epi['L0'], results_epi['Z'], color = colors[k], label = r'$%d$'%N_epi)
		# ax.hlines(np.log(N_epi), 1e6, 1e8, color = colors[k], ls = ':')

	my_plot_layout(ax = axZ, xscale='log', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)#, bottom = 1e-1, top = 1.1)
	axZ.legend(fontsize = 18, title = r'$g_{\textrm{max}}$', title_fontsize = 20)
	figZ.savefig(output_plot+'/Z1_all.pdf')	

	my_plot_layout(ax = ax, xscale='log', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)#, bottom = 1e-1, top = 1.1)
	ax.legend(fontsize = 18, title = r'$g_{\textrm{max}}$', title_fontsize = 20)
	fig.savefig(output_plot+'/g1_all.pdf')

	my_plot_layout(ax = ax_betas, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)#, bottom = 1e-1, top = 1.1)
	ax_betas.legend(fontsize = 18, title = r'$g_{\textrm{max}}$', title_fontsize = 20)
	fig_betas.savefig(output_plot+'/g1_beta_all.pdf')	

	my_plot_layout(ax = ax_L0, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)#, bottom = 1e-1, top = 1.1)
	ax_L0.legend(fontsize = 18, title = r'$g_{\textrm{max}}$', title_fontsize = 20)
	fig_L0.savefig(output_plot+'/g1_L0_all.pdf')		


	# Print Final execution time
	end_time = time.time()
	print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()