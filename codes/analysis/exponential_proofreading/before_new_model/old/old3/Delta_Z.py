import sys
sys.path.append('../../lib/')
from funcs import*
from matplotlib.colors import LogNorm

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=10, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=40, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 0)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**5, help="Number of random sequences.")
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
	subproject = 'DeltaZ'
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"

	output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
	os.makedirs(output_plot, exist_ok=True)

	fig_all, ax_all = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	
	cmap = mpl.colormaps['viridis']
	colors = cmap(np.linspace(0, 1, 15))

	for N_evo in [1]:
		if N_evo == -1:
			N_evo = 'R'
		pars_dir_2 = f"/N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
		antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
		antigens = antigens_data['antigen']

		N_ant = len(antigens)
		
		if new:
			sera = defaultdict(list)
			for kappa, antigen_kappa in enumerate(antigens):
				if kappa == 0:
					data1 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(kappa+1) + '/activated_repertoire.csv')
					for alpha, antigen_alpha in enumerate(tqdm(antigens)):
						sera = run_essay(data1, 'xxxx', 0, antigen_kappa, kappa+1, antigen_alpha, alpha+1, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, time_array, sera)
						data2 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(kappa+1) + '/%d'%(alpha+1) + '/activated_repertoire.csv')
						sera = run_essay(data2, antigen_kappa, kappa+1, antigen_alpha, alpha+1, antigen_alpha, alpha+1, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, time_array, sera)
						sera = run_essay(data2, antigen_kappa, kappa+1, antigen_alpha, alpha+1, antigen_kappa, kappa+1, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, time_array, sera)
						
			sera_df = pd.DataFrame(sera)
			sera_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Z.csv', index = False)


		fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

		sera_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Z.csv')
		for kappa, antigen_kappa in enumerate(antigens):
			if kappa == 0:
				sera_a1 = sera_df.loc[sera_df['past id']==0]
				Z1_homo = sera_a1.loc[(sera_a1['current id']==kappa+1) & (sera_a1['test id']==kappa+1),'Z'].iloc[0]
				for alpha, antigen_alpha in enumerate(antigens):

					d = []
					for i in range(len(antigen_kappa)):
						if antigen_kappa[i]==antigen_alpha[i]:
							d.append(0)
						else:   
							d.append(1)

					Z1_hete = sera_a1.loc[(sera_a1['current id']==kappa+1) & (sera_a1['test id']==alpha+1),'Z'].iloc[0]
					
					sera_a2 = sera_df.loc[sera_df['past id']==kappa+1]
					Z2_homo_1 = sera_a2.loc[(sera_a2['current id']==kappa+1) & (sera_a2['test id']==kappa+1),'Z'].iloc[0]
					Z2_homo_2 = sera_a2.loc[(sera_a2['current id']==alpha+1) & (sera_a2['test id']==alpha+1),'Z'].iloc[0]
					Z2_homo_1_memory = sera_a2.loc[(sera_a2['current id']==kappa+1) & (sera_a2['test id']==kappa+1),'Z memory'].iloc[0]
					Z2_homo_2_memory = sera_a2.loc[(sera_a2['current id']==alpha+1) & (sera_a2['test id']==alpha+1),'Z memory'].iloc[0]

					# if sum(d)<=14:
					# ax_all.scatter(np.log2(Z1_hete/Z1_homo), np.log2(Z2_homo_2/Z2_homo_1), edgecolor = colors[sum(d)], facecolor="None", alpha = .4, marker = '*')
					ax_all.scatter(np.log2(Z1_hete/Z1_homo), np.log2(Z2_homo_2_memory/Z2_homo_1_memory), edgecolor = colors[sum(d)], facecolor="None", alpha = .8, marker = 'D')

		df_pivot = sera_a1.pivot_table(index='test id', columns='current id', values='Z').sort_index(ascending=True)
		diagonal = df_pivot.to_numpy().diagonal()

		# sns.heatmap(df_pivot.divide(diagonal, axis=1), ax = ax, norm=LogNorm(),
		# center = 1e0, vmin = 1e-6, vmax = 1e2, cmap = 'seismic')
		ax_all.plot(np.linspace(-16, 2.5, 10), np.linspace(-16, 2.5, 10), color = 'grey', alpha = 1, ls = '--')

		ax_all.plot(np.linspace(-16, 1, 10), (1 + lamB*p/lamA)*np.linspace(-16, 1, 10) + 3, color = 'k', alpha = 0.8, ls = '-')

	print('modified temperature =', (1 + lamB*p/lamA))
	my_plot_layout(ax = ax_all, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# ax_all.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	ax_all.set_xlim(left = -12, right = 2.5)
	ax_all.set_ylim(bottom = -12, top = 2.5)
	fig_all.savefig(output_plot + '/Z2Z1_'+energy_model+f'_L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}'+'_N_ens-'+str(N_ens)+'_all.pdf')


if __name__ == "__main__":
    main()

