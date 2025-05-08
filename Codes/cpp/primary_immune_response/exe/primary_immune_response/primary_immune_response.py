import sys
sys.path.append('../library/')
from my_functions import*

def main():
    parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
    parser.add_argument('--N_ens', type=int, default=10, help="Number of times to execute the process.")
    parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
    parser.add_argument('--l', type=int, default=10, help="Length of the sequences.")
    parser.add_argument('--t_lim', type=float, default= 6, help="Threshold for the sum of entries.")
    parser.add_argument('--E_lim', type=float, default= -10, help="Threshold for the sum of entries.")
    parser.add_argument('--E_m', type=float, default= -24, help="Threshold for the sum of entries.")
    parser.add_argument('--chunk_size', type=int, default=1000, help="Size of each chunk.")
    parser.add_argument('--p', type=float, default=3, help="# steps.")
    parser.add_argument('--k_step', type=float, default=288, help="step rate.")
    parser.add_argument('--n_jobs', type=int, default=100, help="n_jobs.")
    args = parser.parse_args()


    # Parameters -----------------------------------------------------

    antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
    #antigen = 'TACNSEYPNT'
    l=len(antigen)
    energy_model = 'TCRen'
    #energy_model = 'MJ2'
    N_ens = args.N_ens  # Number of times to execute the process
    L0 = args.L0  # Number of random sequences
    E_lim = args.E_lim  # Threshold for the sum of entries
    t_lim = args.t_lim  # Threshold for the sum of entries
    E_m = args.E_m  #
    chunk_size = args.chunk_size  # Size of each chunk
    p = args.p
    k_step = args.k_step
    n_jobs = args.n_jobs
    #----------------------------------------------------------------
    # Create an output directory if it doesn't exist
    #output_dir = f"/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/primary_immune_response/output_N_ens_{N_ens}_L0_{L0}_p_{p}_k_step_{k_step}_E_lim_{E_lim}_t_lim_{t_lim}"
    output_dir = f"/data01/rmorantovar/text_files/primary_immune_response/output_N_ens_{N_ens}_L0_{L0}_p_{p}_k_step_{k_step}_E_lim_{E_lim}_t_lim_{t_lim}_E_m_{E_m}"
    os.makedirs(output_dir, exist_ok=True)

    # Save the parameters as metadata in a JSON file
    metadata = {
        'N_ens': N_ens,
        'L0': L0,
        'l': l,
        't_lim': t_lim,
        'E_lim': E_lim,
        'E_m': E_m,
        'chunk_size': chunk_size,
        'p': p,
        'k_step': k_step,
        'n_jobs': n_jobs,
    }
    metadata_file = os.path.join(output_dir, 'metadata.json')
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)

    # Update the output CSV file path to be inside the run folder
    output_csv_file = os.path.join(output_dir, 'filtered_sequence_properties.csv')

    PWM = get_motif(antigen, energy_model, '../')
    #Change values by the minimum
    for i in np.arange(l):
        PWM[:,i]-=np.min(PWM[:,i], axis=0)

    Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM, E_m, l)

    # Measure execution time and run the process
    start_time = time.time()
    #process_ensemble(PWM, N_ens, L0, l, K, p, k_step, chunk_size, output_csv_file)
    process_ensemble(np.cumsum(Q0[:-1]*dE)[::100], Es[:-1][::100], N_ens, L0, l, t_lim, E_lim, E_m, p, k_step, chunk_size, output_csv_file, n_jobs = n_jobs)
    end_time = time.time()

    print(f"Execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    # To avoid issues with multiprocessing on Windows, call the main() function within this block
    main()
