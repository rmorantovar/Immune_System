import sys
sys.path.append('../../my_lib/')
from funcs import*
# from classes import*
#from functions_2 import*
from matplotlib.colors import LogNorm

def main():
	# Set up command-line argument parser
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")

	# Define command-line arguments with explanations
	parser.add_argument('--N_ant', type=int, default=100, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=40, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default=0, help="Evolution count.")
	parser.add_argument('--N_epi', type=int, default=3, help="Number of epitopes.")
	parser.add_argument('--L0', type=int, default=10**7, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default=8.0, help="Activation time threshold.")
	parser.add_argument('--E_lim', type=float, default=-6.0, help="Threshold for the sum of entries.")
	parser.add_argument('--E_m', type=float, default=-24, help="Energy threshold.")
	parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=4.0, help="Number of steps.")
	parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
	parser.add_argument('--lamA', type=float, default=6.0, help="Antigen growth rate A.")
	parser.add_argument('--lamB', type=float, default=2.0, help="Antigen growth rate B.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="Number of jobs for parallel processing.")
	parser.add_argument('--random_antigen', type=int, default=0, help="Random antigen flag.")
	parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR', help="Antigen sequence.")
	parser.add_argument('--energy_model', type=str, default='TCRen', help="Energy model.")
	parser.add_argument('--seqs', type=int, default=1, help="Number of sequences.")
	parser.add_argument('--one_WT', type=int, default=0, help="Single WT flag.")
	parser.add_argument('--secondary', type=int, default=0, help="Secondary infection flag.")
	parser.add_argument('--secondary_all', type=int, default=1, help="Secondary all infections flag.")
	parser.add_argument('--add_mutant', type=int, default=0, help="add mutant.")
	parser.add_argument('--pro', type=str, default='epitope_complexity', help="Project name.")
	parser.add_argument('--subpro', type=str, default='fluctuations', help="Subproject name.")
	parser.add_argument('--exp', type=int, default=0, help="Experiment ID.")
	parser.add_argument('--potency_all', type = int, default = 0)

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
	potency_all = args.potency_all
	add_mutant = args.add_mutant

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
	pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv", converters={"antigen": literal_eval})
	antigens = antigens_data['antigen']
	# print(antigens)

	# ------------ ------------ ------------

	start_time = time.time()
	print('Starting simulation ...')

	if one_WT:
		WTs = [antigens.iloc[0]]
	else:
		WTs = list(antigens)
	
	colors = [my_blue, my_red, my_green, my_cyan]
	lss = ['-', '--']
	lws = [1, 2]

	for kappa1, antigen_kappa1 in enumerate(WTs):
		print('primary infection')
		print(kappa1+1)
		output_dir1 = root_dir + pars_dir_1 + pars_dir_2 + "/%d"%(kappa1+1)
		input_file1 = os.path.join(output_dir1, 'activated_repertoire.csv')
		input_file1_DG = os.path.join(output_dir1, 'DG.csv')
		output_file1_pot = os.path.join(output_dir1, 'potency.csv')
		if os.path.isfile(input_file1) and not os.path.isfile(output_file1_pot):
			# ---------------------Calculate motif---------------------
			motif = get_motif(antigen_kappa1, energy_model, '../../')*1.2
			E_ms = np.zeros(N_epi)
			E_rs = np.zeros(N_epi)
			for epi in range(N_epi):
				E_m = -3
				# Normalize motif
				for i in range(l):
					E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
					motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
				# E_m = -23
				E_ms[epi] = E_m
				# Calculate Q0, Es, dE, betas
				Es, dE, Q0, betas = calculate_Q0(0.01, 50, 100000, motif[:, epi*l:(epi+1)*l], E_m, l)
				#--------------------------Repertoire properties--------------------------
				beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
				E_rs[epi] = E_r
				print(np.exp(E_m), np.exp(E_r))
				print(E_m, E_r)

			data_activation = pd.read_csv(input_file1)
			# data_activation['N_t'] = data_activation['N_t'].apply(lambda x: np.array(x, dtype=np.float32))
			if potency_all:
				data_DG = pd.read_csv(input_file1_DG)
				
				potency_dict = {}
				if not os.path.isfile(output_file1_pot):
					for alpha in data_DG.columns[7:].to_numpy():
						col_name = f"{alpha}"
						data_activation['E'] = data_DG[str(alpha)]
						data_activation['Z'] = data_activation['N']/np.exp(data_activation['E'])
						data_activation_mod = data_activation.groupby(['ens_id', 'epi', 'm']).agg({'E':'mean', 'Z':'sum'}).reset_index()
						potency_dict[col_name] = data_activation_mod['Z'].to_numpy()
						
					data_activation_mod = data_activation_mod.drop('Z', axis = 1)
					data_activation_mod = pd.concat([data_activation_mod, pd.DataFrame(potency_dict)], axis=1)
					data_activation_mod.to_csv(output_file1_pot, index = False)

			else:
				output_file1_pot = os.path.join(output_dir1, 'potency_' + str(kappa1+1) + '.csv')
				if not os.path.isfile(output_file1_pot):
					data_activation['Z'] = data_activation['N']/np.exp(data_activation['E'])
					data_activation_mod = data_activation.groupby(['ens_id', 'epi', 'm']).agg({'E':'mean', 
																't':'mean', 
																'Z':'sum'}).reset_index() 
					# data_activation_mod['Z'] = data_activation_mod['Z'].apply(lambda x: list(x))
					data_activation_mod.to_csv(output_file1_pot, index = False)


		if secondary: # Are there seconda infections?
			if secondary_all: # Do I want secondary infections with all the pathogens in the panel?
				for a2, antigen2_ in enumerate(tqdm(antigens)):
					# antigen2 = antigen2_.replace('-', '')
					# antigen2_seq = from_aa_to_i(antigen2, energy_model, '../../')
					output_dir2 = root_dir + pars_dir_1 + pars_dir_2 + "/%d"%(kappa1+1) + "/%d"%(a2+1)
					output_file2 = os.path.join(output_dir2, 'activated_repertoire.csv')
					output_file2_DG = os.path.join(output_dir2, 'DG.csv')

					data_activation = pd.read_csv(output_file2, converters={"N_t": literal_eval})
					data_activation['N_t'] = data_activation['N_t'].apply(lambda x: np.array(x, dtype=np.float32))
					if not potency_all:
						output_file2_pot = os.path.join(output_dir2, 'potency_' + str(a2+1) + '.csv')
						# if not os.path.isfile(output_file2_pot):
						try:
							data_activation['Z_t'] = data_activation['N_t']/np.exp(data_activation['E'])
							data_activation_mod = data_activation.groupby(['ens_id', 'epi', 'm']).agg({'E':'mean', 
																		't':'mean', 
																		'Z_t':'sum'}).reset_index()
							data_activation_mod['Z_t'] = data_activation_mod['Z_t'].apply(lambda x: list(x))
							data_activation_mod.to_csv(output_dir2 + '/potency_' + str(a2+1) + '.csv', index = False)
						except FileNotFoundError:
							print(f'skipping second infection with antigen # {a2+1}')
							continue
					else:
						data_DG = pd.read_csv(output_file2_DG)
						for alpha, antigen_alpha in enumerate((antigens)):
							output_file2_pot = os.path.join(output_dir2, 'potency_' + str(alpha+1) + '.csv')
							# if not os.path.isfile(output_file2_pot):
							try:
								data_activation['E'] = data_DG[str(alpha+1)]
								data_activation['Z_t'] = data_activation['N_t']/np.exp(data_activation['E'])
								data_activation_mod = data_activation.groupby(['ens_id', 'epi', 'm']).agg({'E':'mean', 
																			't':'mean', 
																			'Z_t':'sum'}).reset_index() 
								data_activation_mod['Z_t'] = data_activation_mod['Z_t'].apply(lambda x: list(x))
								data_activation_mod.to_csv(output_dir2 + '/potency_' + str(alpha+1) + '.csv', index = False)

							except FileNotFoundError:
								print(f'skipping primary infection with antigen # {a2+1}')
								continue

						
			else: # Are there secondary infections only with the WT pathogen?
				for inf in tqdm(range(N_inf-1)): # How many WT recurrent infections?
					# antigen2_seq = from_aa_to_i(antigen1, energy_model, '../../')
					
					input_file2 = os.path.join(output_dir1, 'activated_repertoire.csv')
					output_dir1 = output_dir1 + "/%d"%(kappa1+1)
					output_file2 = os.path.join(output_dir1, 'activated_repertoire.csv')

					# if os.path.isfile(input_file2):
					output_file1_pot = os.path.join(output_dir1, 'potency_' + str(kappa1+1) + '.csv')
					if not os.path.isfile(output_file1_pot):
						try:
							# if os.path.isfile(output_file2):
							data_activation = pd.read_csv(output_file2, converters={"N_t": literal_eval})
							data_activation['N_t'] = data_activation['N_t'].apply(lambda x: np.array(x, dtype=np.float32))
							data_activation['Z_t'] = data_activation['N_t']/np.exp(data_activation['E'])
							data_activation_mod = data_activation.groupby(['ens_id', 'epi', 'm']).agg({'E':'mean', 
																		't':'mean', 
																		'Z_t':'sum'}).reset_index()
							data_activation_mod['Z_t'] = data_activation_mod['Z_t'].apply(lambda x: list(x))
							data_activation_mod.to_csv(output_dir1 + '/potency_' + str(kappa1+1) + '.csv', index = False)
						except FileNotFoundError:
							print(f'skipping infection # {inf+2}')
							continue

	

	# Print Final execution time
	end_time = time.time()
	print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()