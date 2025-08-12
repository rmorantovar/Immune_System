from email import parser
import sys
sys.path.append('../../my_lib/')
from funcs import*

'''
------------ Instructions ------------

With this code you can run multiepitope recurrent infection events.
Use the variables one_WT, secondary and secondary_all to costumize your simulation.
For example, a recurrent infection with a single variant uses one_WT = 1, secondary = 1 and secondary_all = 0.
In this case you can use N_inf to set the number of recurrent infections
To run this simulations, you should have created a file antigens.txt in the corresponding folder.

'''
def main():
	# Setting up command-line argument parser
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=80, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = -1)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default=8., help="Threshold for activation time.") # Use 8 for L0>1e6
	parser.add_argument('--E_lim', type=float, default=-7., help="Threshold for the sum of entries.") # Use -6 for L0>1e6
	parser.add_argument('--E_m', type=float, default=-24, help="Threshold for the sum of entries.")
	parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=3.0, help="# steps.")
	parser.add_argument('--pmem', type=float, default=1.0, help="# steps for memory.")
	parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
	parser.add_argument('--lamA', type=float, default=6., help="Antigen growth rate.")
	parser.add_argument('--lamB', type=float, default=2., help="Antigen growth rate.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
	parser.add_argument('--random_antigen', type=int, default=0)
	# parser.add_argument('--antigen', type=str, default='TACNSEYPNTTRAKCGRWYR')
	parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR')
	parser.add_argument('--energy_model', type=str, default = 'TCRen')
	parser.add_argument('--use_seqs', type=int, default = 0)
	parser.add_argument('--one_WT', type=int, default = 1)
	parser.add_argument('--constant_K', type=int, default = 0)
	parser.add_argument('--reuse_repertoire', type=int, default = 0, help="Reuse the activated repertoire from a previous run.")
	parser.add_argument('--secondary', type=int, default = 0)
	parser.add_argument('--secondary_all', type=int, default = 0, help="Do I want secondary infections with all the pathogens in the panel?")
	parser.add_argument('--pro', type=str, default='memory_response', help="project.")
	parser.add_argument('--subpro', type=str, default='Z_NA_dynamics', help="subproject.")
	parser.add_argument('--exp', type=str, default='const_rep', help="experiment.")
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
	if L0>=1e7:
		chunk_size = args.chunk_size
	else:
		chunk_size = args.L0
	p = args.p
	pmem = args.pmem
	k_step = args.k_step
	lamA = args.lamA
	lamB = args.lamB
	n_jobs = args.n_jobs
	random_antigen = args.random_antigen
	antigen = args.antigen
	energy_model = args.energy_model
	
	use_seqs = args.use_seqs
	reuse_repertoire = args.reuse_repertoire
	one_WT = args.one_WT
	secondary = args.secondary
	secondary_all = args.secondary_all
	constant_K = args.constant_K

	if N_evo == -1:
		N_evo = 'R'

	dT = 0.002
	C = 2e4
	T = 12
	time_array = np.linspace(0, T, int((T-0)/dT))
	Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)

	project = args.pro
	subproject = args.subpro
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}"#/{experiment}"
	pars_dir_1 = f"L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
	antigens_path = os.path.join(root_dir, pars_dir_1, pars_dir_2, "antigens.csv")
	antigens = pd.read_csv(antigens_path, converters={"antigen": literal_eval})
	# antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
	# antigens = antigens_data['antigen']
	# print("PANEL OF ANTIGENS:")


	# ------------ ------------ ------------

	start_time = time.time()
	print('Starting simulation ...')

	if one_WT:
		WTs = antigens.iloc[[0]]
	else:
		WTs = antigens
	# 	# WTs = antigens.sample(n=5, replace = False)
	# 	WTs = antigens.iloc[[i*5 for i in range(20)]]

	print(WTs)
	for index, row in WTs.iterrows():
		antigen_kappa1 = row['antigen']
		kappa1 = index
		print('	primary infection with kappa ', kappa1+1)
		print(antigen_kappa1)
		output_dir1 = os.path.join(root_dir, pars_dir_1, pars_dir_2, str(kappa1+1))
		output_file1 = os.path.join(output_dir1, 'activated_repertoire.csv')
		input_file1 = ''

		# Calculate motif
		motif = get_motif(antigen_kappa1, energy_model, '../../')*1.2 # Change depending if antigens are strings or ints

		Q0s, Ess, dEs, Es_r, Es_r0, Es_ms, betas_r = get_Me_repertoire_properties(motif, N_epi, l, L0, constant_K)
		
		print(Es_r, Es_r0, Es_ms)
		sim_params = dict(
			Alphabet=Alphabet, motif=motif, Q0s=Q0s, Ess=Ess, dEs=dEs, Es_ms=Es_ms,
			time_array=time_array, dT=dT,
			N_ens=N_ens, N_epi=N_epi,
			t_lim=t_lim, E_lim=E_lim,
			L0=L0, l=l, p=p, pmem=pmem, k_step=k_step, lamA=lamA, lamB=lamB, C=C,
			infection=1,
			chunk_size=chunk_size,
			input_memory_file=input_file1,
			DDE=0,
			use_seqs=use_seqs,
			reuse_repertoire=reuse_repertoire,
			n_jobs=n_jobs
		)
		if not os.path.isfile(output_file1):
			# Execute process
			df_response = ensemble_of_responses(**sim_params)
			os.makedirs(output_dir1, exist_ok=True)
			df_response = df_response.sort_values(by=['ens_id', 'epi', 't'], ascending=[True, True, True])
			df_response.to_csv(output_file1, index=False)

			# Save parameters for reproducibility
			simple_params = {k: v for k, v in sim_params.items() if not isinstance(v, (np.ndarray, list, tuple))}
			with open(os.path.join(output_dir1, "params.json"), "w") as f:
				json.dump(simple_params, f, indent=2)

		print('	primary infection terminated')
		
		# Modify to work properly with antigens as a df !!!!!!
		if secondary: # Do I want secondary infections? 
			iterate_dict = {
				"DDE": [0.0, 1.0, 2.0, 3.0],
				"pmem": [1.0, 4.0]
			}
			param_names = list(iterate_dict.keys())
			param_values = list(iterate_dict.values())
			for values in itertools.product(*param_values):
				print('	secondary infection...')
				DDE, pmem = values  # Unpack explicitly
				folder_suffix = f"DDE_{DDE}_pmem_{pmem}"
				print(folder_suffix)
				input_file2 = os.path.join(output_dir1, 'activated_repertoire.csv')
				output_dir2 = os.path.join(output_dir1, folder_suffix)
				output_file2 = os.path.join(output_dir2, 'activated_repertoire.csv')
				sim_params.update({
					"input_memory_file": input_file2,
					"infection": 2,
					"DDE": DDE,
					"pmem": pmem
					})
				if os.path.isfile(input_file2):
					if not os.path.isfile(output_file2):
						# Execute process
						df_response = ensemble_of_responses(**sim_params)
						os.makedirs(output_dir2, exist_ok=True)
						df_response = df_response.sort_values(by=['ens_id', 'epi', 't'], ascending=[True, True, True])
						df_response.to_csv(output_file2, index=False)
						# Save parameters for reproducibility
						simple_params = {k: v for k, v in sim_params.items() if not isinstance(v, (np.ndarray, list, tuple))}
						with open(os.path.join(output_dir2, "params.json"), "w") as f:
							json.dump(simple_params, f, indent=2)

				print('	secondary infection terminated')

			# if secondary_all: # Do I want secondary infections with all the pathogens in the panel?
			# 	# for i_DDE, DDE in enumerate(tqdm(np.linspace(0., 8, 9))):				
			# 	for i_DDE, DDE in enumerate(tqdm([0.0])):
			# 		input_file2 = os.path.join(output_dir1, 'activated_repertoire.csv')
			# 		output_dir2 = os.path.join(output_dir1, folder_suffix)
			# 		output_file2 = os.path.join(output_dir2, 'activated_repertoire.csv')
			# 		json_file = os.path.join(output_dir2, "params.json")

			# 		if os.path.isfile(input_file2):
			# 			if not os.path.isfile(output_file2):
			# 				# Execute process
			# 				sim_params.update({
			# 					"input_memory_file": input_file2,
			# 					"infection": 2,
			# 					"DDE": DDE,  # or `inf + 1` in other block
			# 					"pmem": pmem 
			# 				})

			# 				df_response = ensemble_of_responses(**sim_params)
			# 				os.makedirs(output_dir2, exist_ok=True)
			# 				df_response = df_response.sort_values(by=['ens_id', 'epi', 't'], ascending=[True, True, True])
			# 				df_response.to_csv(output_file2, index=False)

			# 				# Save parameters for reproducibility
			# 				with open(os.path.join(output_dir2, "params.json"), "w") as f:
			# 					json.dump(sim_params, f, indent=2)

			# else: # Do I want secondary infections only with the WT pathogen?
			# 	for inf in tqdm(range(N_inf-1)): # How many WT recurrent infections?
			# 		# antigen2 = antigen2_.replace('-', '')
			# 		# antigen2_seq = from_aa_to_i(antigen1, energy_model, '../../')
					
			# 		input_file2 = os.path.join(output_dir1, 'activated_repertoire.csv')
			# 		output_dir2 = os.path.join(output_dir1, "%d"%(kappa1+1))
			# 		output_file2 = os.path.join(output_dir2, 'activated_repertoire.csv')
					
			# 		# THE SECONDARY INFECTION MUST BE MODIFY TO THE NEW ENSEMBLE_OF_RESPONSES() FUNCTION!!!!!!!
			# 		if os.path.isfile(input_file2):
			# 			if not os.path.isfile(output_file2):
			# 				# Execute process
			# 				df_response = ensemble_of_responses(Alphabet, motif, Q0s, Ess, dEs, time_array, dT, N_ens, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, lamB, C, inf+2, chunk_size, input_file2, N_epi)
			# 				os.makedirs(output_dir2, exist_ok=True)
			# 				df_response = df_response.sort_values(by=['ens_id', 'epi', 't'], ascending=[True, True, True])
			# 				df_response.to_csv(output_file2, index=False)
								

	# Print Final execution time
	end_time = time.time()
	print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
	main()
