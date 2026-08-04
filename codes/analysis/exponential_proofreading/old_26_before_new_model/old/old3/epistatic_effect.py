import sys
sys.path.append('../../lib/')
from funcs import*
from matplotlib.colors import LogNorm

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=10, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=5, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 0)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
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

	E_lim = -11.  # Threshold for the sum of entries
	t_lim = 8.  # Threshold for the sum of entries
	chunk_size = 1e6  # Size of each chunk
	p = 3
	k_step = 720
	n_jobs = -1

	lamA = 6.0
	lamB = 3 * np.log(2) #(days)^-1
	lamB = 2.
	dT = 0.05
	C = 1e4
	time_array = np.linspace(0, 15, int((10-0)/dT))
	colors_inf = plt.cm.jet(np.linspace(0,1,N_inf))
	colors_mut = [my_blue, my_red]

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	#energy_model = 'MJ2'
	# antigen = args.antigen
	# epitopes = antigen.split('-')
	# l=len(epitopes[0])
	
	project = 'memory_response'
	subproject = 'multi-epitope'
	subproject = 'Panel'
	subproject = 'epistatic_effect'
	experiment = 1
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"

	output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
	os.makedirs(output_plot, exist_ok=True)
	
	cmap = mpl.colormaps['viridis']
	colors = cmap(np.linspace(0, 1, 15))
	n_muts = [1, 2, 3]
	for N_evo in [1]:
		if N_evo == -1:
			N_evo = 'R'
		pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
		antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")#, converters={"antigen": literal_eval})
		antigens = antigens_data['antigen']

		for i in range(len(antigens)):
			antigens[i] = from_aa_to_i(antigens[i], energy_model, '../../')

		N_ant = len(antigens)
		
		if new:
			for kappa, ag_kappa in enumerate(tqdm(antigens)):
				data1 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(kappa+1) + '/activated_repertoire.csv')
				DDG_df = data1[['seq']].copy()
				for alpha, ag_alpha in enumerate(antigens):
					if hamming_distance(ag_kappa, ag_alpha) in n_muts:
						DDG_df = DDG_distributions(data1['E'], DDG_df, 'xxxx', 0, ag_kappa, kappa+1, ag_alpha, alpha+1, energy_model, N_epi, l)
						# print(DDG_df[DDG_df['seq'] == 'KQCWYPPWVMLYECET'])
				DDG_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(kappa+1) + '/DDGs_epistasis.csv', index = False)

		fig0, ax0 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		DEs_memory_mean = []
		DEs_ms_mean = []

		for n_mut in n_muts:
			fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
			DEs = []
			DEs_ms = []

			for kappa, ag_kappa in enumerate(antigens):

				mseq = find_complementary_seq_min(from_i_to_aa(ag_kappa, energy_model, '../../'), energy_model, '../../')

				# antigen_seq_test = from_aa_to_i(ag_kappa, energy_model, '../../')

				motif = get_motif(ag_kappa, energy_model, '../../')*1.2
				E_ms = np.ones(N_epi)
				for epi in range(N_epi):
					E_m = -3
					# Normalize motif
					for i in range(l):
						E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
						motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
					E_ms[epi] = E_m
					E_kappa = calculate_energy(motif, from_aa_to_i(mseq, energy_model, '../../')) + E_m


				DDG_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(kappa+1) + '/DDGs_epistasis.csv')
				data1 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(kappa+1) + '/activated_repertoire.csv')

				for alpha, ag_alpha in enumerate(antigens):
					if hamming_distance(ag_kappa, ag_alpha) == n_mut:
						# print(kappa, alpha)
						ax.scatter(np.exp(data1['E']), DDG_df[str(alpha+1)], color = my_blue, alpha = 0.8, s = 20)
						DEs = np.concatenate((DEs, DDG_df[str(alpha+1)]))
						# antigen_seq_test = from_aa_to_i(ag_alpha, energy_model, '../../')
						motif = get_motif(ag_alpha, energy_model, '../../')*1.2
						E_ms = np.ones(N_epi)
						for epi in range(N_epi):
							E_m = -3
						    # Normalize motif
							for i in range(l):
								E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
								motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
							E_ms[epi] = E_m
							E_alpha = calculate_energy(motif, from_aa_to_i(mseq, energy_model, '../../')) + E_m

						ax.scatter(np.exp(E_kappa), E_alpha - E_kappa, color = my_red, marker = '*', alpha = 0.8, s = 50)
						DEs_ms = np.concatenate((DEs_ms, [E_alpha - E_kappa]))
				            

			my_plot_layout(ax = ax, xscale='log', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
			# ax.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
			# ax.set_xlim(left = -17 + 11, right = 2.5 + 11)
			# ax.set_ylim(bottom = -17, top = 2.5)
			# ax.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
			# ax.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
			fig.savefig(output_plot + '/Epistatic_effect_'+energy_model+f'_L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}'+'_N_ens-'+str(N_ens)+'_n_mut-%d.pdf'%(n_mut))

			fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
			
			ax.hist(DEs, alpha = 0.8, bins = np.linspace(-6, 15, 16), color = my_blue, density = True)
			ax.hist(DEs_ms, alpha = 0.8, bins = np.linspace(-6, 15, 16), color = my_red, density = True)
			ax.vlines(np.mean(DEs), 0, ax.get_ylim()[1], color = my_blue, lw = 4, ls = '--', label = r'$%.2f$'%np.mean(DEs))
			ax.vlines(np.mean(DEs_ms), 0, ax.get_ylim()[1], color = my_red, lw = 4, ls = '--', label = r'$%.2f$'%np.mean(DEs_ms))
			DEs_memory_mean.append(np.mean(DEs))
			DEs_ms_mean.append(np.mean(DEs_ms))

			my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
			ax.legend(fontsize = 24)
			fig.savefig(output_plot + '/Epistatic_effect_distrib_'+energy_model+f'_L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}'+'_N_ens-'+str(N_ens)+'_n_mut-%d.pdf'%(n_mut))

		ax0.plot(n_muts, DEs_memory_mean, color = my_blue, marker = 'D', lw = 4, label = 'Memory')
		ax0.plot(n_muts, DEs_ms_mean, color = my_red, marker = 'D', lw = 4, label = 'MS')
		my_plot_layout(ax = ax0, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax0.set_xticks(n_muts)
		ax0.legend(fontsize = 24)
		fig0.savefig(output_plot + '/Epistatic_effect_comparison_'+energy_model+f'_L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}'+'_N_ens-'+str(N_ens)+'.pdf')

if __name__ == "__main__":
    main()





