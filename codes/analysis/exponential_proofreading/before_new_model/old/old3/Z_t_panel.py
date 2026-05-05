import sys
sys.path.append('../../lib/')
from funcs import*
from matplotlib.colors import LogNorm

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=100, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=200, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = -1)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**7, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--new', type=int, default=0, help="run Z values again.")
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
	dT = 15./1500.
	C = 1e4
	time_array = np.linspace(0, 15, 1500)
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
	subproject = 'Z_dynamics'
	experiment = 0
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
	antigens = antigens_data['antigen']
	output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
	os.makedirs(output_plot, exist_ok=True)

	N_ant = len(antigens)
	
	if new:
		sera = defaultdict(list)
		for a1, antigen1 in enumerate(tqdm(antigens[:10])):
			antigen_seq = from_aa_to_i(antigen1, energy_model, '../../')
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
			try:
				data1 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/activated_repertoire.csv')
				# for a11, antigen11 in enumerate(antigens):
				sera = run_essay_time(data1, 'xxxx', 0, antigen1, a1+1, antigen1, a1+1, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, beta_r, time_array, sera)
				
				data2 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/%d'%(a1+1) + '/activated_repertoire.csv')
				# for a22, antigen22 in enumerate(antigens):
				sera = run_essay_time(data2, antigen1, a1+1, antigen1, a1+1, antigen1, a1+1, energy_model, N_epi, N_ens, l, p, lamA, lamB, 100*C, dT, beta_r, time_array, sera)

				# for a2, antigen2 in enumerate(antigens):
				# 	try:
				# 		data2 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/%d'%(a2+1) + '/activated_repertoire.csv')
				# 		# for a22, antigen22 in enumerate(antigens):
				# 		sera = run_essay_time(data2, antigen1, a1+1, antigen2, a2+1, antigen2, a2+1, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, time_array, sera)
				# 	except FileNotFoundError:
				# 		print(f'skipping {antigen2}; background {antigen1}')
				# 		continue
			except FileNotFoundError:
				print(f'skipping {antigen1}')
				continue
					
		sera_df = pd.DataFrame(sera)
		sera_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Z_time.csv', index = False)


	fig1, ax1 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	# fig2, ax2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

	sera_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Z_time.csv')
	range_betas = np.max(sera_df['beta'])-np.min(sera_df['beta'])
	sera_a1 = sera_df.loc[sera_df['past id']==0]
	cmap = mpl.colormaps['plasma']
	colors = cmap(np.linspace(0, 1, 51))
	
	Z1_homo = np.zeros_like(time_array)
	Z2_homo = np.zeros_like(time_array)
	counter = 0
	# exponents_time = defaultdict(list)
	for a1, antigen1 in enumerate(antigens[:2]):
			counter+=1
			sera_a2 = sera_df.loc[sera_df['past id']==a1+1]

		# 	antigen_seq = from_aa_to_i(antigen1, energy_model, '../../')
		# 	#--------------------------Energy Motif--------------------------
		# 	motif = get_motif(antigen_seq, energy_model, '../../')*1.2

		# 	#Change values by the minimum
		# 	E_m = -3
		# 	for i in np.arange(l):
		# 		E_m+=np.min(motif[:,i], axis=0)
		# 		motif[:,i]-=np.min(motif[:,i], axis=0)
		# 	# if E_m<-22:
		
		# 	# print('Em:', E_m)
		# 	#--------------------------Entropy function--------------------------
		# 	Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
		# 	Kds = np.exp(Es[:-1])
		# 	#--------------------------Repertoire properties--------------------------
		# 	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
			
			beta1 = sera_a1.loc[(sera_a1['current id']==a1+1) & (sera_a1['test id']==a1+1),'beta'].iloc[0]
			beta2 = sera_a2.loc[(sera_a2['current id']==a1+1) & (sera_a2['test id']==a1+1),'beta'].iloc[0]

			Z1_homo_a1 = np.array(sera_a1.loc[(sera_a1['current id']==a1+1) & (sera_a1['test id']==a1+1),'Z'].apply(literal_eval).iloc[0])
			Z2_homo_a1 = np.array(sera_a2.loc[(sera_a2['current id']==a1+1) & (sera_a2['test id']==a1+1),'Z'].apply(literal_eval).iloc[0])
			
			Z1_homo += np.log(Z1_homo_a1)
			Z2_homo += np.log(Z2_homo_a1)

		# 	# for a2, antigen2 in enumerate(antigens):
		# 	# 	Z1_hete = sera_a1.loc[(sera_a1['current id']==a1+1) & (sera_a1['test id']==a2+1),'Z'].apply(literal_eval).iloc[0]
				
		# 	# 	sera_a2 = sera_df.loc[sera_df['past id']==a1+1]
		# 	# 	Z2_homo = sera_a2.loc[(sera_a2['current id']==a1+1) & (sera_a2['test id']==a1+1),'Z'].apply(literal_eval).iloc[0]
		# 	# 	Z2_hete = sera_a2.loc[(sera_a2['current id']==a2+1) & (sera_a2['test id']==a2+1),'Z'].apply(literal_eval).iloc[0]
			ax1.plot(time_array, Z1_homo_a1, color = my_red, ls = '--', lw = 2, alpha = .8, marker = '')
			Z_i_min = Z1_homo[Z1_homo>0][0]
			Z_i_max = Z1_homo[-1]
			# ax1.plot(time_array[(Z1_homo>100*Z_i_min) & (Z1_homo<0.2*Z_i_max)], Z1_homo[(Z1_homo>100*Z_i_min) & (Z1_homo<0.2*Z_i_max)], color = colors[int(50*(beta1-0.5)/3)], alpha = .2, marker = '.', ls = '')

			ax1.plot(time_array, Z2_homo_a1, color = my_purple, ls = '-', lw = 2, alpha = 1, marker = '')
			Z_i_min = Z2_homo[Z2_homo>0][0]
			Z_i_max = Z2_homo[-1]
			min_time = 3.823823823823824
			# ax2.plot(time_array[(time_array>min_time*1.25) & (Z2_homo<0.5*Z_i_max)], Z2_homo[(time_array>min_time*1.25) & (Z2_homo<0.5*Z_i_max)], color = colors[int(50*(beta2-0.5)/3)], alpha = .2, marker = '.', ls = '')

	Z1_homo = np.exp(Z1_homo/counter)
	Z2_homo = np.exp(Z2_homo/counter)

	# df_pivot = sera_a1.pivot_table(index='test id', columns='current id', values='Z').sort_index(ascending=True)
	# diagonal = df_pivot.to_numpy().diagonal()
	# print(df_pivot, diagonal)

	# sns.heatmap(df_pivot.divide(diagonal, axis=1), ax = ax, norm=LogNorm(),
	# center = 1e0, vmin = 1e-6, vmax = 1e2, cmap = 'seismic')
	
	my_plot_layout(ax = ax1, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# ax1.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# ax1.set_xlim(left = 2e-6, right = 2e2)
	ax1.set_ylim(bottom = 1e5, top = 1e14)
	# ax1.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax1.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig1.savefig(output_plot + '/Z_t_Panel_'+energy_model+'_L0-'+str(L0)+'_N_ant-'+str(N_ant)+'_N_ens-'+str(N_ens)+'_N_evo-'+str(N_evo)+'.pdf')

	# my_plot_layout(ax = ax2, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# # ax2.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# # ax2.set_xlim(left = 2e-6, right = 2e2)
	# ax2.set_ylim(bottom = 1e5, top = 1e13)
	# # ax2.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# # ax2.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	# fig2.savefig(output_plot + '/Z2_t_Panel_'+energy_model+'_L0-'+str(L0)+'_N_ant-'+str(N_ant)+'_N_ens-'+str(N_ens)+'_N_evo-'+str(N_evo)+'.pdf')



if __name__ == "__main__":
    main()

