import sys
sys.path.append('../../lib/')
from func_memory import*
from classes import*
#from functions_2 import*

def main():
    parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
    parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
    parser.add_argument('--L0', type=int, default=10**7, help="Number of random sequences.")
    parser.add_argument('--l', type=int, default=20, help="Length of the sequences.")
    parser.add_argument('--t_lim', type=float, default= 6, help="Threshold for the sum of entries.")
    parser.add_argument('--E_lim', type=float, default= -8, help="Threshold for the sum of entries.")
    parser.add_argument('--E_m', type=float, default= -24, help="Threshold for the sum of entries.")
    parser.add_argument('--chunk_size', type=int, default=100000, help="Size of each chunk.")
    parser.add_argument('--p', type=float, default=4, help="# steps.")
    parser.add_argument('--k_step', type=float, default=288, help="step rate.")
    parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
    args = parser.parse_args()


    # Parameters -----------------------------------------------------

    N_ens = args.N_ens  # Number of times to execute the process
    L0 = args.L0  # Number of random sequences
    E_lim = args.E_lim  # Threshold for the sum of entries
    t_lim = args.t_lim  # Threshold for the sum of entries
    E_m = args.E_m  #
    chunk_size = args.chunk_size  # Size of each chunk
    p = args.p
    k_step = args.k_step
    n_jobs = args.n_jobs
    lambda_A = 6.0
    lambda_B = 3 * np.log(2) #(days)^-1
    dT = 0.05
    C = 1e4
    time_array = np.linspace(0, 10, int((10-0)/dT))
    #----------------------------------------------------------------
    # Create an output directory if it doesn't exist
    output_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/memory/exploration/output_N_ens_{N_ens}_L0_{L0}_p_{p}_k_step_{k_step}_E_lim_{E_lim}_t_lim_{t_lim}_lamA_{lambda_A}"
    #output_dir = f"/data01/rmorantovar/text_files/primary_immune_response/output_N_ens_{N_ens}_L0_{L0}_p_{p}_k_step_{k_step}_E_lim_{E_lim}_t_lim_{t_lim}_E_m_{E_m}"
    os.makedirs(output_dir, exist_ok=True)

    energy_model = 'TCRen'
    #energy_model = 'MJ2'
    antigens = ['TACNSEYPNTTRAKCGRWYC', 'TACNSEYPNTTFDKCGRWYC', 'MACNSEYPNTTRAKCGRWYC', 'MACNSEYPNTTRCKCLRWYC', 'YACNSEYPNTTFDKCGRWYC', 'TACNSTYPNTERAKCGRWYC', 'MACNSEYPNTTRCRKLRWYC',
    'TACNSKYPNDTTDKCGRWSC', 'MACNSECPNTTRCKWLRWYC', 'MACNSEYPNTTRCKEFRWYC'] #L=20
    antigens = ['TACNSEYPNTTRAKCGRWYC']
    
    for antigen in antigens:    
        antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
        l=len(antigen)

        # Update the output CSV file path to be inside the run folder
        output_csv_file = os.path.join(output_dir, 'activated_population_' + antigen + '_2.csv')

        motif = get_motif(antigen_seq, energy_model, '../../')
        #Change values by the minimum
        for i in np.arange(l):
            motif[:,i]-=np.min(motif[:,i], axis=0)

        Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)

        # Measure execution time and run the process
        start_time = time.time()

        #-----------------Loading data----------------------------
        #parameters_path = 'L-%d_Nbc-%d_Antigen-'%(l, L0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
        #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
        return_data_type = 0
        data, return_data_type = get_data(folder_path = '../../out/initial_expansion', data_type = 'ie')

        if(return_data_type==1):
            clone_size_total = data[0]
        else:
            data = pd.read_csv(output_dir + '/activated_population_' + antigen + '.csv')
            clone_size_total = []
            for i_ens in tqdm(np.arange(N_ens)):
                data_active = data.loc[data['ens_id']==i_ens]
                t_act_data = np.min(data_active['time'])
                data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
                activation_times = np.array(data_active['time'])
                energies  = np.array(data_active['E'])
                ids = data_active.index.tolist()

                #---------------------------- B cell linages ----------------------
                clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
                #--------------------------t_C filter-------------------------
                lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
                clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
                ids_C = np.array(ids)[filter_C]
                clone_size_total = np.concatenate((clone_size_total, clone_sizes_C[:,-1]/1))

        #process_ensemble(motif, np.cumsum(Q0[:-1]*dE)[::100], Es[:-1][::100], N_ens, L0, l, t_lim, E_lim, E_m, p, k_step, chunk_size, output_csv_file, n_jobs = n_jobs)
        
        
        final_clone_size = []
        i = 0
        for j in range(len(data)):
            if j in ids_C:
                final_clone_size.append(int(clone_size_total[i]))
                i+=1
            else:
                final_clone_size.append(0)
        
        data['n'] = final_clone_size
        data = data.loc[data['n']!=0]
        data.to_csv(output_csv_file, index=False)

        end_time = time.time()
        print(f"Execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    # To avoid issues with multiprocessing on Windows, call the main() function within this block
    main()

