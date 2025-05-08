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

	cmap = mpl.colormaps['viridis']
	colors = cmap(np.linspace(0, 1, 15))

	fig1, ax1 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	fig2, ax2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	fig3, ax3 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

	L0s = [int(1e4), int(1e5), int(1e6), int(1e7)]
	# L0s = [int(1e6), int(1e7)]

	gammas = []
	E0s = []

	gammas_null = []
	E0s_null = []

	gammas_null2 = []
	E0s_null2 = []

	vars_L0 = []

	gauge_density = 2e15

	homo_memory = []

	betas_star = []

	for L0 in L0s:
		pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
		ys = []
		xs = []
		ys_null = []
		xs_null = []
		ys_null2 = []
		xs_null2 = []

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
					antigen_seq = from_aa_to_i(antigen_kappa, energy_model, '../../')
					#--------------------------Energy Motif--------------------------
					motif = get_motif(antigen_seq, energy_model, '../../')*1.2

					#Change values by the minimum
					E_m = -3
					for i in np.arange(l):
						E_m+=np.min(motif[:,i], axis=0)
						motif[:,i]-=np.min(motif[:,i], axis=0)
					# if E_m<-22:
				
					# print('Em:', E_m)
					#--------------------------Entropy function--------------------------
					Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
					Kds = np.exp(Es[:-1])
					#--------------------------Repertoire properties--------------------------
					beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
					print('beta_r = %.1f'%beta_r, 'K_r = %.1e'%Kd_r, 'E_r = %.1f'%E_r)
					betas_star.append(beta_r)
					#--------------------------Proofreading properties--------------------------
					# beta_step, E_step, Kd_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
					# print('beta_step = %.2f'%beta_step)

					Z1_homo = sera1.loc[(sera1['current id']==kappa+1) & (sera1['test id']==kappa+1),'Z'].iloc[0]
					sera2 = sera_df.loc[sera_df['past id']==kappa+1]
					Z2_homo_1 = sera2.loc[(sera2['current id']==kappa+1) & (sera2['test id']==kappa+1),'Z'].iloc[0]					
					Z2_homo_1_memory = sera2.loc[(sera2['current id']==kappa+1) & (sera2['test id']==kappa+1),'Z memory'].iloc[0]

					homo_memory.append((Z2_homo_1_memory/Z2_homo_1))

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
						T2 = -DDG2#*np.log2(np.exp(1))
						# T2 = -np.log(Z1_hete/Z1_homo)

						Z2_homo_2 = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==alpha+1),'Z'].iloc[0]
						Z2_homo_2_memory = sera2.loc[(sera2['current id']==alpha+1) & (sera2['test id']==alpha+1),'Z memory'].iloc[0]

						if sum(d)<=14:

							if (1/(1+(gauge_density*Z1_hete/N_A)**(-1)))<1.0:
								xs_null.append(T1)
								ys_null.append(1/(1+(gauge_density*Z1_hete/N_A)**(-1)))

							if ((Z2_homo_2_memory/Z2_homo_2)/1)<1.0:
								xs.append(T1)
								ys.append((Z2_homo_2_memory/Z2_homo_2)/1)


		popt, pcov = curve_fit(my_linear_func, ydata = np.log(1/np.array(ys) - 1), xdata = np.array(xs))
		gamma = -1/popt[1]
		E0 = popt[0]*gamma

		gammas.append(gamma)
		E0s.append(E0)

		popt, pcov = curve_fit(my_linear_func, ydata = np.log(1/np.array(ys_null) - 1), xdata = np.array(xs_null))
		gamma_null = -1/popt[1]
		E0_null = popt[0]*gamma_null

		gammas_null.append(gamma_null)
		E0s_null.append(E0_null)

	print(gammas)
	ax1.plot(L0s, np.array(gammas), label = r'$\rm Noneq$', ls = '-', marker = 'D', lw = 3, color = 'mediumorchid')
	ax1.plot(L0s, np.array(gammas_null), label = r'$\rm Null$', ls = '-', marker = 'D', lw = 3, color = my_green)
	ax1.plot(L0s, np.ones_like(L0s)*(-1), label = r'$\rm Heuristic$', color = 'k', ls = '--')
	# ax1.plot(L0s, -1/(3-np.array(betas_star)), color = 'k')
	# ax1.plot(L0s, np.array(gammas_null2), label = r'$\rm Null 2$')
	
	ax2.plot(L0s, E0s, label = r'$\rm Noneq$', ls = '-', marker = 'D', lw = 3, color = 'mediumorchid')
	ax2.plot(L0s, E0s_null, label = r'$\rm Null$', ls = '-', marker = 'D', lw = 3, color = my_green)
	ax2.plot(L0s, np.ones_like(L0s)*(3/np.log2(np.exp(1))), label = r'$\rm Heuristic$', color = 'k', ls = '--')

	ax3.plot(L0s, homo_memory, ls = '-', marker = 'D', lw = 3, color = my_blue2)

	my_plot_layout(ax = ax1, xscale='log', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	ax1.legend(fontsize = 18, title_fontsize = 22, loc = 0)
	# ax1.set_xlim(left = -15 + 11, right = 2. + 11)
	# ax1.set_ylim(bottom = 2e-6)
	# ax1.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax1.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig1.savefig('../../../Figures/memory_response/DeltaZ/gamma_'+energy_model+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

	my_plot_layout(ax = ax2, xscale='log', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	ax2.legend(fontsize = 18, title_fontsize = 22, loc = 0)
	# ax2.set_xlim(left = -15 + 11, right = 2. + 11)
	# ax2.set_ylim(bottom = 2e-6)
	# ax2.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax2.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig2.savefig('../../../Figures/memory_response/DeltaZ/E0_'+energy_model+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')

	my_plot_layout(ax = ax3, xscale='log', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# ax3.legend(fontsize = 18, title_fontsize = 22, loc = 0)
	# ax3.set_xlim(left = -15 + 11, right = 2. + 11)
	# ax3.set_ylim(bottom = 2e-6)
	# ax3.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax3.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig3.savefig('../../../Figures/memory_response/DeltaZ/Memory_homo_'+energy_model+'_N_ens-'+str(N_ens)+'_'+str(experiment)+'.pdf')



if __name__ == "__main__":
    main()

