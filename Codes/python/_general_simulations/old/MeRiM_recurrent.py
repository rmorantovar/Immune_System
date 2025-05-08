import sys
sys.path.append('../../lib/')
from funcs import*

def main():
	# Setting up command-line argument parser
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=40, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=10, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 0)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default=8., help="Threshold for activation time.") # Use 8 for L0>1e6
	parser.add_argument('--E_lim', type=float, default=-6., help="Threshold for the sum of entries.") # Use -6 for L0>1e6
	parser.add_argument('--E_m', type=float, default=-24, help="Threshold for the sum of entries.")
	parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=3, help="# steps.")
	parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
	parser.add_argument('--lamA', type=float, default=6., help="Antigen growth rate.")
	parser.add_argument('--lamB', type=float, default=2., help="Antigen growth rate.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
	parser.add_argument('--random_antigen', type=int, default=0)
	# parser.add_argument('--antigen', type=str, default='TACNSEYPNTTRAKCGRWYR')
	parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR')
	parser.add_argument('--energy_model', type=str, default = 'TCRen')
	parser.add_argument('--seqs', type=int, default = 1)
	parser.add_argument('--pro', type=str, default='epitope_complexity', help="project.")
	parser.add_argument('--subpro', type=str, default='panel', help="subproject.")
	parser.add_argument('--exp', type=int, default=1, help="experiment.")
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

	if N_evo == -1:
		N_evo = 'R'
		
	dT = 0.01
	C = 1e4
	T = 12
	time_array = np.linspace(0, T, int((T-0)/dT))
	Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)

	project = 'memory_response'
	project = args.pro
	subproject = 'multi-epitope'
	subproject = 'Recurrent'
	subproject = args.subpro
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
	antigens = antigens_data['antigen']
	print(antigens)

	# ------------ ------------ ------------
	
	start_time = time.time()
	print('Starting simulation ...')

	for a1, antigen1 in enumerate(antigens):
		if a1==0: # THIS LINE MAKES THIS CODE "RECURRENT"
			print(a1, antigen1)
			antigen1_seq = from_aa_to_i(antigen1, energy_model, '../../')
			output_dir1 = root_dir + pars_dir_1 + pars_dir_2 + "/%d"%(a1+1)

			input_file1 = ''
			output_file1 = os.path.join(output_dir1, 'activated_repertoire.csv')

			# Calculate motif
			motif = get_motif(antigen1_seq, energy_model, '../../')*1.2
			E_ms = np.ones(N_epi)
			for epi in range(N_epi):
				E_m = -3
				# Normalize motif
				for i in range(l):
					E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
					motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
				# print(E_m)
				E_ms[epi] = E_m
				# Calculate Q0, Es, dE, betas
				Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif[:, epi*l:(epi+1)*l], E_m, l)
			if not os.path.isfile(output_file1):
				# Execute process
				df_activation = ensemble_of_activations_Me(Alphabet, motif, np.cumsum(Q0*dE)[::2], Es[:-1][::2], N_ens, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, 1, chunk_size, input_file1, N_epi, seqs = seqs)
				# print(df_activation)
				if len(np.where(np.array([len(df_activation.loc[df_activation['ens_id']==k].index) for k in range(N_ens)]) == 0)[0])==0:	
					df_expansion = ensemble_of_expansions(df_activation, N_ens, p, time_array, lamB, C, dT)
					#df_expansion = df_expansion.groupby(['ens_id', 'E'], as_index=False).agg({'time': 'mean', 'm': 'mean', 'n': 'sum', 'sequence':'first'})
					os.makedirs(output_dir1, exist_ok=True)
					df_expansion.to_csv(output_file1, index=False)

			for inf in tqdm(range(N_inf-1)):

				antigen2_seq = from_aa_to_i(antigen1, energy_model, '../../')
				
				input_file2 = os.path.join(output_dir1, 'activated_repertoire.csv')
				output_dir1 = output_dir1 + "/%d"%(a1+1)
				output_file2 = os.path.join(output_dir1, 'activated_repertoire.csv')

				if os.path.isfile(input_file2):
					if not os.path.isfile(output_file2):
						# Execute process
						df_activation = ensemble_of_activations_Me(Alphabet, motif, np.cumsum(Q0*dE)[::2], Es[:-1][::2], N_ens, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, 2, chunk_size, input_file2, N_epi, seqs = seqs)
						if len(np.where(np.array([len(df_activation.loc[df_activation['ens_id']==k].index) for k in range(N_ens)]) == 0)[0])==0:
							df_expansion = ensemble_of_expansions(df_activation, N_ens, p, time_array, lamB, C, dT)
							#df_expansion = df_expansion.groupby(['ens_id', 'E'], as_index=False).agg({'time': 'mean', 'm': 'mean', 'n': 'sum', 'sequence':'first'})
							os.makedirs(output_dir1, exist_ok=True)
							df_expansion.to_csv(output_file2, index=False)

	# Print Final execution time
	end_time = time.time()
	print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()
