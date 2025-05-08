import sys
sys.path.append('../../lib/')
from functions_1 import*
from classes import*
#from functions_2 import*

def main():
    # Setting up command-line argument parser
    parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
    parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
    parser.add_argument('--N_a', type=int, default=1, help="Number of antigens.")
    parser.add_argument('--L0', type=int, default=10**7, help="Number of random sequences.")
    parser.add_argument('--l', type=int, default=20, help="Length of the sequences.")
    parser.add_argument('--t_lim', type=float, default=6, help="Threshold for the sum of entries.")
    parser.add_argument('--E_lim', type=float, default=-10, help="Threshold for the sum of entries.")
    parser.add_argument('--E_m', type=float, default=-24, help="Threshold for the sum of entries.")
    parser.add_argument('--chunk_size', type=int, default=100000, help="Size of each chunk.")
    parser.add_argument('--p', type=float, default=4, help="# steps.")
    parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
    parser.add_argument('--lamA', type=float, default=6.0, help="Antigen growth rate.")
    parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
    parser.add_argument('--seqs', type=int, default=0)
    args = parser.parse_args()

    # Parameters
    energy_model = 'TCRen'
    # energy_model = 'MJ2'  # Uncomment this line if using the MJ2 energy model
    # Define antigen sequence
    antigen = 'TACNSEYPNTTRAKCGRWYC'  # L=20
    # antigen = 'TACNSEYTRAGRWYC'  # L=15
    antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
    l = len(antigen)

    # Retrieving arguments
    N_ens = args.N_ens
    N_a = args.N_a
    L0 = args.L0
    E_lim = args.E_lim
    t_lim = args.t_lim
    E_m = args.E_m
    chunk_size = args.chunk_size
    p = args.p
    k_step = args.k_step
    lamA = args.lamA
    n_jobs = args.n_jobs
    seqs = args.seqs

    # Create output directory if it doesn't exist
    project = 'memory'
    subproject = 'exploration'
    output_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/output_N_ens_{N_ens}_L0_{L0}_p_{p}_k_step_{k_step}_E_lim_{E_lim}_t_lim_{t_lim}_lamA_{lamA}"
    os.makedirs(output_dir, exist_ok=True)

    # Save parameters as metadata in a JSON file
    metadata = {
        'N_ens': N_ens,
        'N_a': N_a,
        'L0': L0,
        'l': l,
        't_lim': t_lim,
        'E_lim': E_lim,
        'E_m': E_m,
        'chunk_size': chunk_size,
        'p': p,
        'k_step': k_step,
        'lamA': lamA,
        'n_jobs': n_jobs,
        'seqs': seqs
    }

    # Generate antigens
    antigens = [antigen]
    for _ in range(N_a - 1):
        temp_antigen = np.random.choice(antigens)
        temp_antigen_seq = from_aa_to_i(temp_antigen, energy_model, '../../')
        for _ in range(2):
            temp_antigen_seq[np.random.randint(0, l)] = np.random.randint(0, 20)
        temp_antigen = from_i_to_aa(temp_antigen_seq, energy_model, '../../')
        antigens.append(temp_antigen)

    # Iterate over antigens
    for antigen in antigens:
        # Save metadata
        metadata_file = os.path.join(output_dir, 'metadata_' + antigen + '.json')
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        # Update output CSV file paths
        output_csv_file = os.path.join(output_dir, 'activated_population_with_memory_' + antigen + '.csv')
        input_memory_file = os.path.join(output_dir, 'activated_population_' + antigen + '_2.csv')

        # Calculate motif
        motif = get_motif(antigen_seq, energy_model, '../../')
        # Normalize motif
        for i in range(l):
            motif[:, i] -= np.min(motif[:, i], axis=0)

        # Calculate Q0, Es, dE, betas
        Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)

        # Execute process
        start_time = time.time()
        process_ensemble_2(motif, np.cumsum(Q0[:-1]*dE)[::100], Es[:-1][::100], N_ens, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, chunk_size, input_memory_file, output_csv_file)
        end_time = time.time()

        # Print execution time
        print(f"Execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()
