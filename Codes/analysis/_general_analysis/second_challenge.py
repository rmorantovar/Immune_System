import sys
sys.path.append('../../lib/')
from funcs import*
# from classes import*
#from functions_2 import*
from matplotlib.colors import LogNorm

def main():
	# Set up command-line argument parser
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")

	# Define command-line arguments with explanations
	parser.add_argument('--N_ant', type=int, default=4, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default=0, help="Evolution count.")
	parser.add_argument('--N_epi', type=int, default=2, help="Number of epitopes.")
	parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default=8.0, help="Activation time threshold.")
	parser.add_argument('--E_lim', type=float, default=-6.0, help="Threshold for the sum of entries.")
	parser.add_argument('--E_m', type=float, default=-24, help="Energy threshold.")
	parser.add_argument('--chunk_size', type=int, default=100000000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=3, help="Number of steps.")
	parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
	parser.add_argument('--lamA', type=float, default=6.0, help="Antigen growth rate A.")
	parser.add_argument('--lamB', type=float, default=2.0, help="Antigen growth rate B.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="Number of jobs for parallel processing.")
	parser.add_argument('--random_antigen', type=int, default=0, help="Random antigen flag.")
	parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR', help="Antigen sequence.")
	parser.add_argument('--energy_model', type=str, default='TCRen', help="Energy model.")
	parser.add_argument('--seqs', type=int, default=1, help="Number of sequences.")
	parser.add_argument('--one_Ref', type=int, default=1, help="Single Ref flag.")
	parser.add_argument('--secondary', type=int, default=1, help="Secondary infection flag.")
	parser.add_argument('--secondary_all', type=int, default=1, help="Secondary all infections flag.")
	parser.add_argument('--pro', type=str, default='epitope_complexity', help="Project name.")
	parser.add_argument('--subpro', type=str, default='Panel', help="Subproject name.")
	parser.add_argument('--exp', type=int, default=5, help="Experiment ID.")
    
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
	one_Ref = args.one_Ref
	secondary = args.secondary
	secondary_all = args.secondary_all

	if N_evo == -1:
		N_evo = 'R'

	dT = 0.001
	C = 1e4
	T = 12
	time_array = np.linspace(0, T, int((T-0)/dT))
	Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)

	project = args.pro
	subproject = args.subpro
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"/N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv", converters={"antigen": literal_eval})
	antigens = antigens_data['antigen']

	# ------------ ------------ ------------

	start_time = time.time()
	print('Starting simulation ...')

	if one_Ref:
		Refs = [antigens.iloc[0]]
	else:
		Refs = list(antigens)
	
	colors = [my_blue, my_red, my_green, my_cyan]
	cmap = mpl.colormaps['rainbow']
	colors = cmap(np.linspace(0, 1, 13))
	lss = ['-', '--']
	lws = [1, 2]
	fig0, ax0 = plt.subplots(figsize=(8, 8), gridspec_kw={'left':0.15, 'right':.95, 'bottom':.08, 'top': 0.95})
	fig, ax = plt.subplots(figsize=(8, 8), gridspec_kw={'left':0.15, 'right':.95, 'bottom':.08, 'top': 0.95})
	pre_drops = []
	post_drops = []
	markers_analysis = []
	colors_analysis = []
	alphas_analysis = []
	for a_ref, antigen_ref in enumerate(tqdm(Refs)):
		# fig_ref, ax_ref = plt.subplots(figsize=(8, 8), gridspec_kw={'left':0.15, 'right':.95, 'bottom':.08, 'top': 0.95})
		x = np.linspace(-8, 8, 100)
		# ax_ref.plot(x, x, ls = '--', color = 'k', alpha = .8)
		# ax_ref.vlines(0, -8, 8, ls = '--', color = 'k', alpha = .8)
		# ax_ref.hlines(0, -8, 8, ls = '--', color = 'k', alpha = .8)
		ax0.plot(x, x, ls = '--', color = 'k', alpha = .8)
		ax0.vlines(0, -8, 8, ls = '--', color = 'k', alpha = .8)
		ax0.hlines(0, -8, 8, ls = '--', color = 'k', alpha = .8)
		ax.plot(x, x, ls = '--', color = 'k', alpha = .8)
		ax.vlines(0, -8, 8, ls = '--', color = 'k', alpha = .8)
		ax.hlines(0, -8, 8, ls = '--', color = 'k', alpha = .8)
		
		for a_serum, antigen_serum in enumerate(antigens):
			#-------- calculate distance ---------
			d_ref_ser = 0
			for i in range(len(antigen_serum)):
				if antigen_serum[i]!=antigen_ref[i]:  
					d_ref_ser+=1
			#-------- ------------------ ---------
			output_dir1 = root_dir + pars_dir_1 + pars_dir_2 + "/%d"%(a_serum+1)
			output_dir2 = root_dir + pars_dir_1 + pars_dir_2 + "/%d"%(a_serum+1) + "/%d"%(a_ref+1)

			for a_test, antigen_test in enumerate((antigens)):
				#-------- calculate distance ---------
				d_ref_test = 0
				for i in range(len(antigen_test)):
					if antigen_test[i]!=antigen_ref[i]:  
						d_ref_test+=1
				#-------- ------------------ ---------

				output_file1 = os.path.join(output_dir1, 'potency_' + str(a_test+1) + '.csv')
				output_file1_ref = os.path.join(output_dir1, 'potency_' + str(a_ref+1) + '.csv')
				
				try:
					data_activation = pd.read_csv(output_file1, converters={"Z_t": literal_eval})
					data_activation_ref = pd.read_csv(output_file1_ref, converters={"Z_t": literal_eval})
					pre_test = 0
					pre_ref = 0
					for index, row in data_activation.iterrows():
						pre_test+=row['Z_t'][-1]
					for index, row in data_activation_ref.iterrows():
						pre_ref+=row['Z_t'][-1]

					pre_drop = np.log(pre_test/pre_ref)
					pre_drops = np.append(pre_drops, pre_drop)
				except FileNotFoundError:
					print(f'skipping primary infection with antigen # {a_serum+1}')
					continue
	
				output_file2 = os.path.join(output_dir2, 'potency_' + str(a_test+1) + '.csv')
				output_file2_ref = os.path.join(output_dir2, 'potency_' + str(a_ref+1) + '.csv')

				try:
					data_activation = pd.read_csv(output_file2, converters={"Z_t": literal_eval})
					data_activation_ref = pd.read_csv(output_file2_ref, converters={"Z_t": literal_eval})
					post_test = 0
					post_ref = 0
					for index, row in data_activation.iterrows():
						post_test+=row['Z_t'][-1]
					for index, row in data_activation_ref.iterrows():
						post_ref+=row['Z_t'][-1]			
					post_drop = np.log(post_test/post_ref)
					post_drops = np.append(post_drops, post_drop)
				except FileNotFoundError:
					print(f'skipping second infection with antigen # {a_test+1}')
					continue

				if (a_serum == a_test) and (a_serum != a_ref):
					colors_analysis.append(int(d_ref_ser/N_epi))
					markers_analysis.append('*')
					alphas_analysis.append(0.8)
					# ax_ref.scatter(pre_drop, post_drop, marker = '*', color = my_blue , s = (13 - d_ref_ser)*10, alpha = .8)
					# ax.scatter(pre_drop, post_drop, marker = '*', color = colors[int(d_ref_ser/N_epi)] , s = 50, alpha = .8)
				elif a_serum == a_ref:
					colors_analysis.append(int(d_ref_ser/N_epi))
					markers_analysis.append('s')
					alphas_analysis.append(0.8)
					# ax.scatter(pre_drop, post_drop, marker = 's', color = colors[int(d_ref_ser/N_epi)] , s = 30, alpha = .8)
				else:
					colors_analysis.append('grey')
					markers_analysis.append('.')
					alphas_analysis.append(0.08)
					# ax_ref.scatter(pre_drop, post_drop, marker = '.', color = my_green , s = (13 - d_ref_ser)*10, alpha = .5)
					# ax.scatter(pre_drop, post_drop, marker = '.', color = 'grey' , s = 55, alpha = .05)

				

			# ax.scatter(np.mean(pre_drops), np.mean(post_drops), marker = 'o', color = 'k' , s = 30, alpha = .3)
				
		# my_plot_layout(ax = ax_ref, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30, bottom = -8, left = -8, top = 8, right = 8)
		# output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
		# os.makedirs(output_plot, exist_ok=True)
		# fig_ref.savefig(output_plot +'/pre_post_ref-'+str(a_ref+1)+'_L0-'+str(L0)+'_N_epi-'+str(N_epi)+'_all.pdf')
		# plt.close(fig_ref)
		
	pre_drops = np.array(pre_drops)
	post_drops = np.array(post_drops)
	markers_analysis = np.array(markers_analysis)
	colors_analysis = np.array(colors_analysis)
	alphas_analysis = np.array(alphas_analysis)

	ax0.scatter(pre_drops, post_drops, marker = '.', c = 'grey', s = 55, alpha = .5)

	ax.scatter(pre_drops[markers_analysis=='.'], post_drops[markers_analysis=='.'], marker = '.', c = colors_analysis[markers_analysis=='.'] , s = 55, alpha = .08)
	ax.scatter(pre_drops[markers_analysis=='*'], post_drops[markers_analysis=='*'], marker = '*', c = [colors[int(i)] for i in colors_analysis[markers_analysis=='*']] , s = 50, alpha = .8)
	ax.scatter(pre_drops[markers_analysis=='s'], post_drops[markers_analysis=='s'], marker = 's', c = [colors[int(i)] for i in colors_analysis[markers_analysis=='s']] , s = 25, alpha = .7)

	my_plot_layout(ax = ax0, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30, bottom = -8, left = -8, top = 8, right = 8)
	output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
	os.makedirs(output_plot, exist_ok=True)
	fig0.savefig(output_plot +f'/pre_post_L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}'+'_N_epi-'+str(N_epi)+'_null.pdf')
	plt.close(fig0)

	my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30, bottom = -8, left = -8, top = 8, right = 8)
	output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
	os.makedirs(output_plot, exist_ok=True)
	fig.savefig(output_plot +f'/pre_post_L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}'+'_N_epi-'+str(N_epi)+'_clusters.pdf')
	plt.close(fig)

	# Print Final execution time
	end_time = time.time()
	print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")

	# # ---------------------Calculate motif---------------------
			# motif = get_motif(antigen_serum, energy_model, '../../')*1.2
			# E_ms = np.zeros(N_epi)
			# E_rs = np.zeros(N_epi)
			# for epi in range(N_epi):
			# 	E_m = -3
			# 	# Normalize motif
			# 	for i in range(l):
			# 		E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], ax_refis=0)
			# 		motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], ax_refis=0)
			# 	E_ms[epi] = E_m
			# 	# Calculate Q0, Es, dE, betas
			# 	Es, dE, Q0, betas = calculate_Q0(0.01, 50, 100000, motif[:, epi*l:(epi+1)*l], E_m, l)
			# 	#--------------------------Repertoire properties--------------------------
			# 	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
			# 	E_rs[epi] = E_r
			# print(np.exp(E_ms), np.exp(E_rs))

			# colors = [colors[c] for c in np.argsort(E_rs)]
			# Kstar_dom = np.exp(-np.min(E_rs))

if __name__ == "__main__":
    main()