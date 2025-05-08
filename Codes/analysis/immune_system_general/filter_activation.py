import sys
sys.path.append('../library/')
from functions import*

dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
k_pr = k_pr*24 #days^-1
Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

# Define the number and type of columns in the binary file
num_cols = 5
col_types = "diidi"

# Calculate the size of each row
row_size = struct.calcsize(col_types)

for p in tqdm([3.0]):
    for N_r in [1e8]:
        for N_ens in [1025]:
            parameters_path = 'L-%d_Nbc-%d_Antigen-'%(20, N_r)+'TACNSEYPNTTRAKCGRWYC'+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(6.0, 3.0, k_pr/24, p, 100000000, 0, N_ens)+'TCRen'

            with open(Text_files_path + 'Dynamics/Ensemble/L%d/'%20+parameters_path+'/energies_ensemble.bin', 'r+b') as f:
                with open(Text_files_path + 'Dynamics/Ensemble/L%d/'%20+parameters_path+'/energies_ensemble_active.bin', 'wb') as f2:
                    #dataraw = np.fromfile(f, dtype=dt)
                    #print('bin file has been read...')
                    #data = pd.DataFrame(dataraw)
                    #print(len(np.array((data.loc[data['active']==1]['active']))))
                    #print(len(np.array((data.loc[data['active']==0]['active']))))
                    #f.seek(0)

                    # Keep track of the current position in the file and the number of rows that have been removed
                    rows_removed = 0
                    # Iterate over all the rows in the file
                    while True:
                        # Read a row from the file
                        row = f.read(row_size)
                        # If we reached the end of the file, stop
                        if not row:
                            break
                        # Unpack the row into its columns
                        col1, col2, col3, col4, col5 = struct.unpack(col_types, row)
                        # Check if the value of the column we're interested in meets the condition
                        if col2 == 1:
                            f2.write(row)
                            pos_write = f2.tell()
                            
                        #else:
                            # If the condition is met, we want to remove this row from the file
                            # To do this, we will overwrite the current row with the next row in the file
                            #next_row = f.read(row_size)
                            #f.seek(pos)
                            #f.write(next_row)
                            #f.seek(pos + row_size)
                            # Increment the number of rows that have been removed
                            #rows_removed += 1                        

                    # Update the number of rows in the header to reflect the rows that were removed