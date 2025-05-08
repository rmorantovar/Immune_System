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

	if N_evo == -1:
		N_evo = 'R'


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
	
	L0s = [int(1e4), int(1e5), int(1e6), int(1e7)]
	for l, L0 in enumerate(L0s):
		pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"

		
		fig_all, ax_all = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		fig_all_d, ax_all_d = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		Z_drop_d = np.zeros(13)
		Z_drop_d_counter = np.zeros(13)
		for N_evo in [1]:
			if N_evo == -1:
				N_evo = 'R'
			pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
			antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
			antigens = antigens_data['antigen']
			
			cmap = mpl.colormaps['viridis']
			colors = cmap(np.linspace(0, 1, 13))

			Z_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Z.csv')
			DDG_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/DDG.csv')
			
			Z_df = Z_df.loc[Z_df['past id']==0]

			DDG_df = DDG_df.loc[DDG_df['past id']==0]

			for kappa, antigen_kappa in enumerate(antigens):
				if kappa + 1 == 1:
					Z1_homo = Z_df.loc[(Z_df['current id']==kappa+1) & (Z_df['test id']==kappa+1),'Z'].iloc[0]

					for alpha, antigen_alpha in enumerate(antigens):
						d = []
						for i in range(len(antigen_kappa)):
							if antigen_kappa[i]==antigen_alpha[i]:
								d.append(0)
							else:   
								d.append(1)

						Z1_hete = Z_df.loc[(Z_df['current id']==kappa+1) & (Z_df['test id']==alpha+1),'Z'].iloc[0]

						DDG = DDG_df.loc[(DDG_df['current id']==kappa+1) & (DDG_df['test id']==alpha+1),'DDG'].iloc[0]
						var_DDG = DDG_df.loc[(DDG_df['current id']==kappa+1) & (DDG_df['test id']==alpha+1),'var_DDG'].iloc[0]

						# np.log2(Z1_hete)/np.log2(np.exp(1)) - var_DDG/4

						
						if alpha > kappa:
							if sum(d)<=12:
								Z_drop_d[sum(d)]+= np.log2(Z1_hete/Z1_homo)/np.log2(np.exp(1)) - var_DDG/4
								Z_drop_d_counter[sum(d)]+= 1

								ax_all.scatter(DDG, -np.log2(Z1_hete/Z1_homo)/np.log2(np.exp(1)) + var_DDG/4, color = colors[sum(d)], alpha = .8, marker = 's')
								ax_all_d.scatter(sum(d), np.log2(Z1_hete/Z1_homo)/np.log2(np.exp(1)) - var_DDG/4, color = colors[sum(d)], alpha = .8, marker = 'o')


		ax_all_d.plot(np.linspace(0, 12, 13),  Z_drop_d/Z_drop_d_counter, lw = 3, color = my_blue)

		ax_all.plot(np.linspace(-1.5, 14, 10),  (np.linspace(-1.5, 14, 10)), color = 'k')
		ax_all_d.plot(np.linspace(0, 10, 100),  -(np.linspace(0, 10, 100)), color = my_blue2)

		my_plot_layout(ax = ax_all, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax_all.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
		# ax_all.set_xlim(left = 7e-2, right = 3e1)
		# ax_all.set_ylim(bottom = 2e-6, top = 2e2)
		# ax_all.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax_all.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig_all.savefig('../../../Figures/memory_response/DeltaZ/Z_DDG_Panel_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ant-'+str(N_ant)+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

		my_plot_layout(ax = ax_all_d, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax_all_d.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
		# ax_all_d.set_xlim(left = 7e-2, right = 3e1)
		# ax_all_d.set_ylim(bottom = 2e-6, top = 2e2)
		# ax_all_d.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax_all_d.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig_all_d.savefig('../../../Figures/memory_response/DeltaZ/Z_d_Panel_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ant-'+str(N_ant)+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')




if __name__ == "__main__":
    main()

