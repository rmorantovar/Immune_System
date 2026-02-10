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
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	
	L0s = [int(1e4), int(1e5), int(1e6), int(1e7)]

	for l, L0 in enumerate(L0s):
		pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"

		fig_protection_null, ax_protection_null = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		fig_protection, ax_protection = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

		cmap = mpl.colormaps['viridis']
		colors = cmap(np.linspace(0, 1, 10))


		ys_null = []
		xs_null = []

		ys_null2 = []
		xs_null2 = []

		ys = []
		xs = []

		gauge_density = 2e15

		for N_evo in [1]:
			if N_evo == -1:
				N_evo = 'R'
			pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
			antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
			antigens = antigens_data['antigen']

			sera_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Z.csv')
			DDG_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/DDG.csv')
			DDG_df1 = DDG_df.loc[DDG_df['past id']==0]

			sera1 = sera_df.loc[sera_df['past id']==0]

			for kappa, antigen_kappa in enumerate(antigens):
				if kappa+1 == 1:
					Z1_homo = sera1.loc[(sera1['current id']==kappa+1) & (sera1['test id']==kappa+1),'Z'].iloc[0]
					sera2 = sera_df.loc[sera_df['past id']==kappa+1]
					Z2_homo_1 = sera2.loc[(sera2['current id']==kappa+1) & (sera2['test id']==kappa+1),'Z'].iloc[0]
					Z2_homo_1_memory = sera2.loc[(sera2['current id']==kappa+1) & (sera2['test id']==kappa+1),'Z memory'].iloc[0]
					
					ax_protection_null.scatter(-0, (1/(1+(gauge_density*Z1_homo/N_A)**(-1))), color = my_green, alpha = 1, marker = 'o', s = 100)
					ax_protection.scatter(-0, (Z2_homo_1_memory/Z2_homo_1), color = colors[0], alpha = 1, marker = 'o', s = 100)

					DDG_df2 = DDG_df.loc[DDG_df['past id']==kappa+1]

					for alpha, antigen_alpha in enumerate(tqdm(antigens)):
						d = []
						for i in range(len(antigen_kappa)):
							if antigen_kappa[i]==antigen_alpha[i]:
								d.append(0)
							else:   
								d.append(1)

						Z1_hete = sera1.loc[(sera1['current id']==kappa+1) & (sera1['test id']==alpha+1),'Z'].iloc[0]

						mean_DDG1 = DDG_df1.loc[(DDG_df1['current id']==kappa+1) & (DDG_df1['test id']==alpha+1),'DDG'].iloc[0]
						var_DDG1 = DDG_df1.loc[(DDG_df1['current id']==kappa+1) & (DDG_df1['test id']==alpha+1),'var_DDG'].iloc[0]
						DDG1 = mean_DDG1
						DDG1 = mean_DDG1 - var_DDG1/4
						T1 = DDG1#*np.log2(np.exp(1))
						# T1 = -np.log(Z1_hete/Z1_homo)

						mean_DDG2 = DDG_df2.loc[(DDG_df2['current id']==alpha+1) & (DDG_df2['test id']==kappa+1),'DDG'].iloc[0]
						var_DDG2 = DDG_df2.loc[(DDG_df2['current id']==alpha+1) & (DDG_df2['test id']==kappa+1),'var_DDG'].iloc[0]
						DDG2 = mean_DDG2
						DDG2 = mean_DDG2 - var_DDG2/4
						T2 = DDG2#*np.log2(np.exp(1))
						# T2 = -np.log(Z1_hete/Z1_homo)

						Z2_homo_2 = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==alpha+1),'Z'].iloc[0]
						Z2_homo_2_memory = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==alpha+1),'Z memory'].iloc[0]
						Z2_homo_2_naive = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==alpha+1),'Z naive'].iloc[0]


						if sum(d)<=10:
							ax_protection_null.scatter(T1, 1/(1+(gauge_density*Z1_hete/N_A)**(-1)) , color = my_green, alpha = .8, marker = 'o', s = 80)
							
							# ax_protection.scatter(T, (Z2_homo_2_memory/Z2_homo_2)/((Z2_homo_1_memory/Z2_homo_1)), color = colors[sum(d)], alpha = .8, marker = '^')
							ax_protection.scatter(T1, (Z2_homo_2_memory/Z2_homo_2)/1, color = 'mediumorchid', alpha = 1, marker = '*', s = 100)
							# ax_protection.scatter(T1, (Z2_homo_2_naive/Z2_homo_2)/1, color = 'aqua', alpha = .6, marker = '^', s = 60)

							if (1/(1+(gauge_density*Z1_hete/N_A)**(-1)))<1.0:
								xs_null.append(T1)
								ys_null.append(1/(1+(gauge_density*Z1_hete/N_A)**(-1)))

							if ((Z2_homo_2_memory/Z2_homo_2)/1)<1.0:
								xs.append(T1)
								ys.append((Z2_homo_2_memory/Z2_homo_2)/1)
					
		DDG_array = np.linspace(-1, 9, 100)
		titer_array = np.linspace(1, 12, 100)

		# ax_protection_null.plot(DDG_array, np.ones_like(DDG_array), color = 'k')
		
		ax_protection_null.plot(DDG_array, 1.0*(1/(1 + np.exp((DDG_array*np.log2(np.exp(1)) - (3))/1.4))), color = 'k', ls = '--', label = r'$ {\rm Heuristic: }E_0 = %.1f ; \gamma = %.1f$'%(3/np.log2(np.exp(1)), -1))
		popt2, pcov2 = curve_fit(my_linear_func, ydata = np.log(1/np.array(ys_null) - 1), xdata = np.array(xs_null))
		gamma = 1/popt2[1]
		E0 = popt2[0]*gamma
		ax_protection_null.plot(DDG_array, 1/(1+np.exp((1/gamma)*(DDG_array + E0))), color = 'k', ls = '-', label = r'$ {\rm Fit: } E_{0_1} = %.1f ; \gamma_1 = %.2f$'%(E0, gamma))

		my_plot_layout(ax = ax_protection_null, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax_protection_null.legend(fontsize = 16, title_fontsize = 18, loc = 1)
		ax_protection_null.set_xlim(left = -1.5, right = 9.5)
		ax_protection_null.set_ylim(bottom = -0.05, top = 1.05)
		# ax_protection_null.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax_protection_null.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig_protection_null.savefig('../../../Figures/memory_response/DeltaZ/Protection_null_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

		# my_plot_layout(ax = ax_protection_null, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax_protection_null.legend(fontsize = 18, title_fontsize = 22, loc = 'best')
		# ax_protection_null.set_xlim(left = -10.5, right = 1.5)
		# # ax_protection_null.set_ylim(bottom = 2e-6)
		# # ax_protection_null.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# # ax_protection_null.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		# fig_protection_null.savefig('../../../Figures/memory_response/DeltaZ/Protection_null_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_all_log.pdf')
		

		# ax_protection.plot(DDG_array, np.ones_like(DDG_array), color = 'k')
		ax_protection.plot(DDG_array, 1.0*(1/(1 + np.exp((DDG_array*np.log2(np.exp(1)) - (3))/1.4))), color = 'k', ls = '--', label = r'$ {\rm Heuristic: }E_0 = %.1f ; \gamma = %.1f$'%(3/np.log2(np.exp(1)), -1))
		popt2, pcov2 = curve_fit(my_linear_func, ydata = np.log(1/(np.array(ys)) - 1), xdata = np.array(xs))
		gamma = 1/popt2[1]
		E0 = popt2[0]*gamma
		ax_protection.plot(DDG_array, 1.0/(1+np.exp((1/gamma)*(DDG_array + E0))), color = 'k', ls = '-', label = r'$ {\rm Fit: } E_0 = %.1f ; \gamma = %.2f$'%(E0, gamma))
		ax_protection.plot(DDG_array, 1.0/(1+np.exp((1/0.5)*(DDG_array - 3))), color = my_blue2, ls = '-', label = r'$ {\rm Theory\, min: } E_0 = %.1f ; \gamma = %.2f$'%(1.7, 0.5))
		ax_protection.plot(DDG_array, 1.0/(1+np.exp((1/0.5)*(DDG_array - 1.7))), color = my_blue, ls = '-', label = r'$ {\rm Theory\, max: } E_0 = %.1f ; \gamma = %.2f$'%(3.0, 0.5))
		# ax_protection.plot(DDG_array, DDG_array, color = 'k', ls = '-', label = r'$ {\rm Fit: } E_0 = %.1f ; \gamma = %.3f$'%(E0, gamma))

		my_plot_layout(ax = ax_protection, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax_protection.legend(fontsize = 16, title_fontsize = 18, loc = 1)
		ax_protection.set_xlim(left = -1.5, right = 9.5)
		ax_protection.set_ylim(bottom = -0.05, top = 1.05)
		# ax_protection.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax_protection.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig_protection.savefig('../../../Figures/memory_response/DeltaZ/Protection_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

		my_plot_layout(ax = ax_protection, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax_protection.legend(fontsize = 16, title_fontsize = 18, loc = 1)
		ax_protection.set_xlim(left = -1.5, right = 9.5)
		ax_protection.set_ylim(bottom = -0.05, top = 1.05)
		# ax_protection.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax_protection.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig_protection.savefig('../../../Figures/memory_response/DeltaZ/Protection_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'_log.pdf')

		# ax_protection.plot(DDG_array[:-40], np.exp(popt1[0])*np.exp(popt1[1]*DDG_array[:-40]), color = 'k', ls = '-', label = r'$ \rm{ Fit: } \gamma = %f$'%(1/popt1[1]))

		# my_plot_layout(ax = ax_protection, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax_protection.legend(fontsize = 18, title_fontsize = 22, loc = 5)
		# ax_protection.set_xlim(left = -10.5 + 11, right = 1.5 + 11)
		# # ax_protection.set_ylim(bottom = 2e-6)
		# # ax_protection.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# # ax_protection.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		# fig_protection.savefig('../../../Figures/memory_response/DeltaZ/Protection_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_all_log.pdf')


if __name__ == "__main__":
    main()

