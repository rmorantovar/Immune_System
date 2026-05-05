import sys
sys.path.append('../../lib/')
from functions_memory import*
# from classes import*
#from functions_2 import*
from matplotlib.colors import LogNorm

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=10, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=40, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 1)
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

	n_lines = 4 + 3*N_evo

	

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
	subproject = 'PS'
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	

	for N_evo in [1]:
		if N_evo == -1:
			N_evo = 'R'
			n_lines = 17
		pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
		antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
		antigens = antigens_data['antigen']

		N_ant = len(antigens)
		
		if new:
			sera = defaultdict(list)
			for kappa, antigen_kappa in enumerate(antigens):
				if kappa == 0:
					data1 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(kappa+1) + '/activated_repertoire.csv')
					for alpha, antigen_alpha in enumerate(tqdm(antigens)):
						sera = DDG(data1, 'xxxx', 0, antigen_kappa, kappa+1, antigen_alpha, alpha+1, energy_model, N_epi, l, sera)
						data2 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(kappa+1) + '/%d'%(alpha+1) + '/activated_repertoire.csv')
						sera = DDG(data2, antigen_kappa, kappa+1, antigen_alpha, alpha+1, antigen_alpha, alpha+1, energy_model, N_epi, l, sera)
						sera = DDG(data2, antigen_kappa, kappa+1, antigen_alpha, alpha+1, antigen_kappa, kappa+1, energy_model, N_epi, l, sera)

			sera_df = pd.DataFrame(sera)
			sera_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/DDG.csv', index = False)


		fig1, ax1 = plt.subplots(figsize=(5, 5), gridspec_kw={'left':0.1, 'right':.9, 'bottom':.1, 'top': 0.9})
		fig2a, ax2a = plt.subplots(figsize=(5, 5), gridspec_kw={'left':0.1, 'right':.9, 'bottom':.1, 'top': 0.9})
		fig2b, ax2b = plt.subplots(figsize=(5, 5), gridspec_kw={'left':0.1, 'right':.9, 'bottom':.1, 'top': 0.9})

		ax1.plot(np.linspace(-2, n_lines, 10), np.linspace(-2, n_lines, 10), color = 'k')
		ax2a.plot(np.linspace(-2, n_lines, 10), np.linspace(-2, n_lines, 10), color = 'k')
		ax2a.plot(np.linspace(-1, 4, 10), 3-1*np.linspace(-1, 4, 10), color = 'k', alpha = .9, ls = ':')
		ax2a.hlines(0, -2, 5, color = 'k', alpha = .9, ls = '--')
		ax2a.vlines(0, -2, 5, color = 'k', alpha = .9, ls = '--')

		ax2b.plot(np.linspace(-2, n_lines, 10), np.linspace(-2, n_lines, 10), color = 'k')
		ax2b.plot(np.linspace(-1, 4, 10), 3-1*np.linspace(-1, 4, 10), color = 'k', alpha = .9, ls = ':')
		ax2b.hlines(0, -2, 5, color = 'k', alpha = .9, ls = '--')
		ax2b.vlines(0, -2, 5, color = 'k', alpha = .9, ls = '--')

		sera_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/DDG.csv')
		sera1 = sera_df.loc[sera_df['past id']==0]
		for a1, antigen1 in enumerate(antigens):
			if a1 == 0:
				for a2, antigen2 in enumerate(antigens):
					d = []
					for i in range(len(antigen1)):
						if antigen1[i]==antigen2[i]:
							d.append(0)
						else:   
							d.append(1)
					if a2>=a1:
						DeltaDeltaGa1a2 = sera1.loc[(sera1['current id']==a1+1) & (sera1['test id']==a2+1),'DDG'].iloc[0]
						DeltaDeltaGa2a1 = sera1.loc[(sera1['current id']==a2+1) & (sera1['test id']==a1+1),'DDG'].iloc[0]

						var_DeltaDeltaGa1a2 = sera1.loc[(sera1['current id']==a1+1) & (sera1['test id']==a2+1),'var_DDG'].iloc[0]
						var_DeltaDeltaGa2a1 = sera1.loc[(sera1['current id']==a2+1) & (sera1['test id']==a1+1),'var_DDG'].iloc[0]

						ax1.scatter(x = DeltaDeltaGa1a2, y = DeltaDeltaGa2a1, color = my_blue, alpha = .8, s = 15*sum(d))
						# ax1.annotate(sum(d), (DeltaDeltaGa1a2, DeltaDeltaGa2a1), fontsize = 8)
		
		# Take colors at regular intervals spanning the colormap.
		cmap = mpl.colormaps['plasma']
		colors = cmap(np.linspace(0, 1, n_lines))
		d1max = 0
		d2amax = 0
		for a, antigen in enumerate(antigens):
			sera2 = sera_df.loc[sera_df['past id']==a+1]
			for a1, antigen1 in enumerate(antigens):
				d1 = []
				for i in range(len(antigen)):
					if antigen[i]==antigen1[i]:
						d1.append(0)
					else:   
						d1.append(1)
				if sum(d1)>d1max:
					d1max = sum(d1)
				for a2, antigen2 in enumerate(antigens):
					d2 = []
					for i in range(len(antigen1)):
						if antigen[i]==antigen2[i]:
							d2.append(0)
						else:   
							d2.append(1)

					d = []
					for i in range(len(antigen1)):
						if antigen1[i]==antigen2[i]:
							d.append(0)
						else:   
							d.append(1)
					if sum(d)>d2amax:
						d2amax = sum(d)
					

					if a2>=a1:
						DeltaDeltaGa1a2 = sera2.loc[(sera2['current id']==a1+1) & (sera2['test id']==a2+1),'DDG'].iloc[0]
						DeltaDeltaGa2a1 = sera2.loc[(sera2['current id']==a2+1) & (sera2['test id']==a1+1),'DDG'].iloc[0]

						var_DeltaDeltaGa1a2 = sera2.loc[(sera2['current id']==a1+1) & (sera2['test id']==a2+1),'var_DDG'].iloc[0]
						var_DeltaDeltaGa2a1 = sera2.loc[(sera2['current id']==a2+1) & (sera2['test id']==a1+1),'var_DDG'].iloc[0]

						# ax2a.errorbar(x = DeltaDeltaGa1a2, y = DeltaDeltaGa2a1, yerr = np.sqrt(var_DeltaDeltaGa2a1)/2, xerr = np.sqrt(var_DeltaDeltaGa1a2)/2, color = my_red, alpha = .4)
						if sum(d1)>=0:
							ax2a.scatter(x = DeltaDeltaGa1a2, y = DeltaDeltaGa2a1, color = colors[sum(d1)], s = 15*sum(d), alpha = .8)
							ax2b.scatter(x = DeltaDeltaGa1a2, y = DeltaDeltaGa2a1, color = colors[sum(d2)], s = 15*sum(d), alpha = .8)
							# if sum(d)==1:
							# 	ax2a.scatter(x = DeltaDeltaGa1a2, y = DeltaDeltaGa2a1, color = 'black', alpha = 1-.2*(sum(d1)), s = 6*sum(d), marker = '*')
							# 	ax2b.scatter(x = DeltaDeltaGa1a2, y = DeltaDeltaGa2a1, color = 'black', alpha = 1-.2*(sum(d2)), s = 6*sum(d), marker = '*')
		print(d1max, d2amax)		
		# df_pivot = sera_a1.pivot_table(index='test id', columns='current id', values='Z').sort_index(ascending=True)
		# diagonal = df_pivot.to_numpy().diagonal()
		# print(df_pivot, diagonal)

		# sns.heatmap(df_pivot.divide(diagonal, axis=1), ax = ax1, norm=LogNorm(),
		# center = 1e0, vmin = 1e-6, vmax = 1e2, cmap = 'seismic')

		legend_elements_2a = [Line2D([0], [0], ls = '', marker='o', color=colors[i], label='%d'%i, markerfacecolor=colors[i], markersize=10, alpha = .8) for i in range(n_lines-1)]
		legend_elements_2b= [Line2D([0], [0], ls = '', marker='o', color=colors[i], label='%d'%i, markerfacecolor=colors[i], markersize=10, alpha = .8) for i in range(n_lines-1)]


		my_plot_layout(ax = ax1, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax1.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
		# ax1.set_xlim(left = 2e-6, right = 2e2)
		# ax1.set_ylim(bottom = 2e-6, top = 2e2)
		ax1.set_xticks(range(-2, 7, 2))
		ax1.set_yticks(range(-2, 7, 2))
		# ax1.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig1.savefig('../../../Figures/memory_response/DeltaZ/DDG_Panel_1_'+energy_model+'_N_ens-'+str(N_ens)+'_N_evo-'+str(N_evo)+'.pdf')

		my_plot_layout(ax = ax2a, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax2a.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
		ax2a.legend(handles=legend_elements_2a, loc=0)
		# ax2a.set_xlim(left = 2e-6, right = 2e2)
		# ax2a.set_ylim(bottom = 2e-6, top = 2e2)
		ax2a.set_xticks(range(-2, 7, 2))
		ax2a.set_yticks(range(-2, 7, 2))
		# ax2a.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig2a.savefig('../../../Figures/memory_response/DeltaZ/DDG_Panel_2a_'+energy_model+'_N_ens-'+str(N_ens)+'_N_evo-'+str(N_evo)+'.pdf')

		my_plot_layout(ax = ax2b, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax2b.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
		ax2b.legend(handles=legend_elements_2b, loc=0)
		# ax2b.set_xlim(left = 2e-6, right = 2e2)
		# ax2b.set_ylim(bottom = 2e-6, top = 2e2)
		ax2b.set_xticks(range(-2, 7, 2))
		ax2b.set_yticks(range(-2, 7, 2))
		# axba.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig2b.savefig('../../../Figures/memory_response/DeltaZ/DDG_Panel_2b_'+energy_model+'_N_ens-'+str(N_ens)+'_N_evo-'+str(N_evo)+'.pdf')


if __name__ == "__main__":
    main()

