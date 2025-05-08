import sys
sys.path.append('../../lib/')
from functions_memory import*
# from classes import*

def main():
	# Setting up command-line argument parser
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 0)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=20, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default=8., help="Threshold for activation time.")
	parser.add_argument('--E_lim', type=float, default=-11., help="Threshold for the sum of entries.")
	parser.add_argument('--E_m', type=float, default=-24, help="Threshold for the sum of entries.")
	parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=3, help="# steps.")
	parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
	parser.add_argument('--lamA', type=float, default=6., help="Antigen growth rate.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
	parser.add_argument('--random_antigen', type=int, default=0)
	# parser.add_argument('--antigen', type=str, default='TACNSEYPNTTRAKCGRWYR')
	parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR')
	parser.add_argument('--energy_model', type=str, default = 'TCRen')
	parser.add_argument('--seqs', type=int, default = 1)
	args = parser.parse_args()


	# Retrieving arguments
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
	chunk_size = args.chunk_size
	p = args.p
	k_step = args.k_step
	lamA = args.lamA
	n_jobs = args.n_jobs
	random_antigen = args.random_antigen
	antigen = args.antigen
	energy_model = args.energy_model
	
	seqs = args.seqs

	lamB = 3 * np.log(2) #(days)^-1
	dT = 0.05
	C = 1e4
	time_array = np.linspace(0, 10, int((10-0)/dT))

	Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)

	for a in range(N_ant):
		if random_antigen:
			antigen = ''.join(np.random.choice(Alphabet[:-1], l*N_epi))
			antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
		else:
			antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
			l = len(antigen)

		# Create output directory if it doesn't exist
		project = 'memory'
		subproject = 'multi-epitope'
		subproject = 'PS'
		output_dir = (f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/N_ens_{N_ens}_N_evo_{N_evo}_N_inf_{N_inf}_N_epi_{N_epi}_L0_{L0}_l_{l}_p_{p}/" + "".join([antigen[epi*l:(epi+1)*l] + '-' for epi in range(N_epi)]))[:-1]
		os.makedirs(output_dir, exist_ok=True)

		# Save parameters as metadata in a JSON file
		metadata = {
			'N_ant': N_ant,
			'N_ens': N_ens,
			'N_inf': N_inf,
			'N_evo' : N_evo,
			'N_epi' : N_epi,
			'L0': L0,
			'l': l,
			't_lim': t_lim,
			'E_lim': E_lim,
			'chunk_size': chunk_size,
			'p': p,
			'k_step': k_step,
			'lamA': lamA,
			'n_jobs': n_jobs,
			'antigen': antigen,
			'energy_model': energy_model,
			'seqs' : seqs,
			'random_antigen': random_antigen
		}
		# Save metadata
		metadata_file = os.path.join(output_dir, 'metadata.json')
		with open(metadata_file, 'w') as f:
		    json.dump(metadata, f, indent=2)

		start_time = time.time()
		print('Starting simulation ...')
		# Iterate over antigens
		for infection in range(1, N_inf+1):
			start_time_inf = time.time()

			# Update CSV file paths
			input_file = os.path.join(output_dir, 'act_pop_' + str(infection-1) + '.csv')
			output_file = os.path.join(output_dir, 'act_pop_' + str(infection) + '.csv')
			
			fig, ax = plt.subplots()
			
			mutation = {}
			mutation['antigen_old'] = ("".join([antigen[epi*l:(epi+1)*l] + '-' for epi in range(N_epi)]))[:-1]

			# Simulate antigen evolution
			if (infection>1):
				for j in range(N_evo):
					pos_mut = np.random.randint(0, l*N_epi)
					temp_array = list(range(20))
					temp_array.remove(antigen_seq[pos_mut])
					antigen_seq[pos_mut] = np.random.choice(temp_array)
					antigen = from_i_to_aa(antigen_seq, energy_model, '../../')
					mutation['pos_mut_%d'%j] = pos_mut + 1
					mutation['epi_mut_%d'%j] = (pos_mut)//l + 1
					mutation['antigen_new_%d'%j] = ("".join([antigen[epi*l:(epi+1)*l] + '-' for epi in range(N_epi)]))[:-1]

			# Save mutation data
			mutation_file = os.path.join(output_dir, 'mutation_'+str(infection)+'.json')
			with open(mutation_file, 'w') as f:
			    json.dump(mutation, f, indent=2)

			# Calculate motif
			motif = get_motif(antigen_seq, energy_model, '../../')*1.2
			E_ms = np.ones(N_epi)

			for epi in range(N_epi):
				E_m = -3
				# Normalize motif
				for i in range(l):
					E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
					motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
				
				E_ms[epi] = E_m
				print('Master sequence binding energy = ', E_m, r'$%.1e$'%(np.exp(E_m)))

				# Calculate Q0, Es, dE, betas
				Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif[:, epi*l:(epi+1)*l], E_m, l)
				ax.plot(np.exp(Es), Q0*L0, label = '%d'%(epi+1))

			ax.hlines(1, 1e-8, 1e3, color = 'k', ls = '--')
			ax.set_yscale('log')
			ax.set_xscale('log')
			ax.legend()
			fig.savefig(output_dir+'/densities_%d.pdf'%infection)
			plt.close()

			# Execute process
			df_activation = ensemble_of_activations_Me(Alphabet, motif, np.cumsum(Q0[:-1]*dE)[::2], Es[:-1][::2], N_ens, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, infection, chunk_size, input_file, N_epi, seqs = seqs)
			df_expansion = ensemble_of_expansions(df_activation, N_ens, p, time_array, lamB, C, dT)
			#df_expansion = df_expansion.groupby(['ens_id', 'E'], as_index=False).agg({'time': 'mean', 'm': 'mean', 'n': 'sum', 'sequence':'first'})
			df_expansion.to_csv(output_file, index=False)

			# Print Final execution time
			end_time = time.time()
			print(f"execution time of infection {infection}: {(end_time - start_time_inf)/60:.2f} minutes")


		# Print Final execution time
		end_time = time.time()
		print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()
