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
	subproject = 'PS'
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	
	fig_epsilon, ax_epsilon = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	L0s = [int(1e4), int(1e5), int(1e6), int(1e7)]


	for l, L0 in enumerate(L0s):
		pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"

		fig_DZ2, ax_DZ2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		fig_DZ1, ax_DZ1 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		# fig_DDG, ax_DDG = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		# fig_DDG_var, ax_DDG_var = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		fig_DZ21, ax_DZ21 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

		cmap = mpl.colormaps['viridis']
		colors = cmap(np.linspace(0, 1, 10))

		ys = []
		xs = []

		means_DDG = []
		ds = []

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

			for kappa, antigen_kappa in enumerate(tqdm(antigens)):
				if kappa+1 == 1:
					Z1_homo = sera1.loc[(sera1['current id']==kappa+1) & (sera1['test id']==kappa+1),'Z'].iloc[0]
					sera2 = sera_df.loc[sera_df['past id']==kappa+1]
					Z2_homo_1 = sera2.loc[(sera2['current id']==kappa+1) & (sera2['test id']==kappa+1),'Z'].iloc[0]
					Z2_homo_1_memory = sera2.loc[(sera2['current id']==kappa+1) & (sera2['test id']==kappa+1),'Z memory'].iloc[0]
					Z2_homo_1_naive = sera2.loc[(sera2['current id']==kappa+1) & (sera2['test id']==kappa+1),'Z naive'].iloc[0]
					
					ax_DZ2.scatter(0, np.log2(Z2_homo_1_memory/Z2_homo_1), color = my_blue, alpha = 1, marker = '*', s = 120, label = 'Memory')
					ax_DZ2.scatter(0, np.log2(Z2_homo_1_naive/Z2_homo_1), color = my_gold, alpha = 1, marker = '^', s = 100, label = 'Naive')
					ax_DZ2.scatter(0, np.log2(Z2_homo_1/Z2_homo_1), color = my_green, alpha = 1, marker = 'o', s = 100, label = 'Total')

					ax_DZ21.scatter(-np.log2(Z1_homo/Z1_homo), np.log2((Z2_homo_1_memory)/(Z1_homo)), color = 'mediumorchid', alpha = .8, marker = '*', s = 120, label = 'Memory')
					ax_DZ21.scatter(-np.log2(Z1_homo/Z1_homo), np.log2((Z2_homo_1_naive)/(Z1_homo)), color = 'aqua', alpha = .8, marker = '^', s = 100, label = 'Naive')
					ax_DZ21.scatter(-np.log2(Z1_homo/Z1_homo), np.log2((Z2_homo_1)/(Z1_homo)), color = colors[0], alpha = .8, marker = 'o', s = 100, label = 'Total')

					DDG_df2 = DDG_df.loc[DDG_df['past id']==kappa+1]

					for alpha, antigen_alpha in enumerate(antigens):
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

						# ax_DDG.scatter(mean_DDG1, mean_DDG2, color = colors[sum(d)], alpha = .8, marker = 'o', s = 60)
						# ax_DDG_var.scatter(var_DDG1, var_DDG2, color = colors[sum(d)], alpha = .8, marker = 'o', s = 60)

						Z2_homo_2 = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==alpha+1),'Z'].iloc[0]
						Z2_homo_2_memory = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==alpha+1),'Z memory'].iloc[0]
						Z2_homo_2_naive = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==alpha+1),'Z naive'].iloc[0]

						Z2_hete_2 = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==kappa+1),'Z'].iloc[0]
						Z2_hete_2_memory = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==kappa+1),'Z memory'].iloc[0]
						Z2_hete_2_naive = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==kappa+1),'Z naive'].iloc[0]

						ax_DZ1.scatter(DDG1, np.log2(Z1_hete/Z1_homo)/np.log2(np.exp(1)), alpha = .8, marker = 'o', s = 60, color = colors[2])
						# ax_DZ1.scatter(DDG2, np.log2(Z2_hete_2/Z2_homo_2)/np.log2(np.exp(1)), alpha = .8, marker = 'D', s = 60, color = colors[8])

						ax_DZ2.scatter(DDG1, np.log2(Z2_homo_2_memory/Z2_homo_1), color = my_blue, alpha = .8, marker = '*', s = 80)
						ax_DZ2.scatter(DDG1, np.log2(Z2_homo_2_naive/Z2_homo_1), color = my_gold, alpha = .8, marker = '^', s = 60)
						ax_DZ2.scatter(DDG1, np.log2(Z2_homo_2/Z2_homo_1), color = my_green, alpha = .8, marker = 'o', s = 60)

						ax_DZ21.scatter(DDG1, np.log2((Z2_homo_2_memory/Z2_homo_1)/(Z1_hete/Z1_homo)), color = 'mediumorchid', alpha = .8, marker = '*', s = 80)
						ax_DZ21.scatter(DDG1, np.log2((Z2_homo_2_naive/Z2_homo_1)/(Z1_hete/Z1_homo)), color = 'aqua', alpha = .8, marker = '^', s = 60)
						ax_DZ21.scatter(DDG1, np.log2((Z2_homo_2/Z2_homo_1)/(Z1_hete/Z1_homo)), color = colors[2], alpha = .8, marker = 'o', s = 60)

						ax_epsilon.scatter(sum(d), mean_DDG1, color = colors[3*l],  alpha = .8, marker = 'o')
						# ax_epsilon.scatter(sum(d), mean_DDG2, color = colors[3*l],  alpha = .8, marker = '^')

						means_DDG.append(mean_DDG1)
						ds.append(sum(d))

						if ((Z2_homo_2_memory/Z2_homo_2)/1)<1.0:
							xs.append(T1)
							ys.append((Z2_homo_2_memory/Z2_homo_2)/1)

		d_array = np.linspace(0, 8, 100)

		popt, pcov = curve_fit(my_linear_func, ydata = means_DDG, xdata = ds)
		ax_epsilon.plot(d_array, popt[0] + d_array*popt[1], color = colors[3*l], ls = '-', label = f"{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}")

		DDG_array = np.linspace(-2.5, 10, 100)
		titer_array = np.log2(np.linspace(2e-4, 3, 100))
		titer_array = -np.log2(np.exp(DDG_array - var_DDG2/4))

		ax_DZ1.plot(DDG_array, -DDG_array, color = 'k', ls = '--')

		my_plot_layout(ax = ax_DZ1, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax_DZ1.legend(fontsize = 18, title_fontsize = 22, loc = 'best')
		ax_DZ1.set_xlim(left = -2.8, right = 9.0)
		ax_DZ1.set_ylim(bottom = -9.0, top = 1.5)
		# ax_DZ1.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax_DZ1.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig_DZ1.savefig('../../../Figures/memory_response/DeltaZ/DZ1_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

		ax_DZ2.plot(DDG_array, -DDG_array, color = 'k', ls = '--')

		my_plot_layout(ax = ax_DZ2, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax_DZ2.legend(fontsize = 18, title_fontsize = 22, loc = 'best')
		ax_DZ2.set_xlim(left = -2.8, right = 9.0)
		ax_DZ2.set_ylim(bottom = -9.0, top = 1.8)
		# ax_DZ2.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax_DZ2.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig_DZ2.savefig('../../../Figures/memory_response/DeltaZ/DZ2_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

		ax_DZ21.plot(DDG_array, np.zeros_like(DDG_array), color = 'k', ls = '--')

		my_plot_layout(ax = ax_DZ21, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax_DZ21.legend(fontsize = 18, title_fontsize = 22, loc = 4)
		ax_DZ21.set_xlim(left = -2.8, right = 9.0)
		ax_DZ21.set_ylim(bottom = -6, top = 6)
		# ax_DZ21.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax_DZ21.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig_DZ21.savefig('../../../Figures/memory_response/DeltaZ/DZ21_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')


		# ax_DZ2.plot(DDG_array[:-40], np.exp(popt1[0])*np.exp(popt1[1]*DDG_array[:-40]), color = 'k', ls = '-', label = r'$ \rm{ Fit: } \gamma = %f$'%(1/popt1[1]))

		# my_plot_layout(ax = ax_DZ2, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# ax_DZ2.legend(fontsize = 18, title_fontsize = 22, loc = 5)
		# ax_DZ2.set_xlim(left = -10.5 + 11, right = 1.5 + 11)
		# # ax_DZ2.set_ylim(bottom = 2e-6)
		# # ax_DZ2.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# # ax_DZ2.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		# fig_DZ2.savefig('../../../Figures/memory_response/DeltaZ/Protection_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_all_log.pdf')


		# ax_DDG.plot(DDG_array, DDG_array, color = 'k', ls = '--')

		# my_plot_layout(ax = ax_DDG, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# # ax_DDG.legend(fontsize = 18, title_fontsize = 22, loc = 'best')
		# # ax_DDG.set_xlim(left = np.log2(1e-4), right = np.log2(8))
		# # ax_DDG.set_ylim(bottom = np.log2(1e-4), top = np.log2(8))
		# # ax_DDG.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# # ax_DDG.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		# fig_DDG.savefig('../../../Figures/memory_response/DeltaZ/DDG21_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

		# my_plot_layout(ax = ax_DDG_var, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# # ax_DDG_var.legend(fontsize = 18, title_fontsize = 22, loc = 'best')
		# # ax_DDG_var.set_xlim(left = np.log2(1e-4), right = np.log2(8))
		# # ax_DDG_var.set_ylim(bottom = np.log2(1e-4), top = np.log2(8))
		# # ax_DDG_var.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# # ax_DDG_var.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		# fig_DDG_var.savefig('../../../Figures/memory_response/DeltaZ/DDG_var_21_'+energy_model+'_L0-1e'+str(int(np.log10(L0)))+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

		

	my_plot_layout(ax = ax_epsilon, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	ax_epsilon.legend(fontsize = 18, title_fontsize = 22, loc = 0)
	ax_epsilon.set_xlim(left = -0.5, right = 10.5)
	ax_epsilon.set_ylim(bottom = -0.8, top = 12)
	# ax_epsilon.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax_epsilon.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig_epsilon.savefig('../../../Figures/memory_response/DeltaZ/epsilon_'+energy_model+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')



if __name__ == "__main__":
    main()

