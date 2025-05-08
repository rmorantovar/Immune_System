import sys
sys.path.append('../../lib/')
from functions_1 import*
from classes import*
import numpy as np
import pickle
from tqdm import tqdm
import pandas as pd
import math
import matplotlib.pyplot as plt
import ast

def main():
    parser = argparse.ArgumentParser(description="Simulate affinity maturation")
    parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
    parser.add_argument('--N', type=int, default=20, help="Number of random sequences.")
    parser.add_argument('--l', type=int, default=20, help="Length of the sequences.")
    parser.add_argument('--b', type=float, default=1, help="Birth rate.")
    parser.add_argument('--d', type=float, default=1, help="Death rate.")
    parser.add_argument('--m', type=float, default=0.01, help="Mutation rate.")
    parser.add_argument('--T', type=float, default=10, help="Simulation time.")
    parser.add_argument('--pop0', type=bool, default=False, help="Simulation time.")

    parser.add_argument('--chunk_size', type=int, default=1000, help="Size of each chunk.")
    parser.add_argument('--p', type=float, default=3, help="# steps.")
    parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")

    args = parser.parse_args()
    # Parameters for the simulation.
    N = args.N
    l = args.l
    b = args.b
    d = args.d
    m = args.m
    T = args.T
    pop0 = args.pop0
    sequences_pop0 = []
    if pop0:
        pop0_df = pd.read_csv('../../in/populations/filtered_sequence_properties_2.csv')
        #print(sequences_pop0_df)
        pop0_df['sequence'] = pop0_df['sequence'].apply(ast.literal_eval)
        sequences_pop0 = pop0_df['sequence'].tolist()
        ns_pop0 = pop0_df['n'].tolist()
        random_clones_index = np.random.choice(len(sequences_pop0), size = 50, replace = False, p = ns_pop0/np.sum(ns_pop0))
        sequences_pop0 = np.array(sequences_pop0)[random_clones_index]
        ns_pop0 = np.array(ns_pop0)[random_clones_index]
        N_clones = len(sequences_pop0)
        N_cells = np.sum(ns_pop0)

    energy_model = 'TCRen'
    antigen = 'TACNSEYPNTTRAKCGRWYC'
    antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
    l = len(antigen)    
    motif = get_motif(antigen_seq, energy_model, '../../')
    
    #Change values by the minimum
    for i in np.arange(l):
        motif[:,i]-=np.min(motif[:,i], axis=0)

    # Number of trajectories in the ensemble.
    N_ens = args.N_ens

    for i in range(N_ens):
        # Create and simulate a population.
        population = Population(N_clones, N_cells, l, b, d, m, motif, pop0 = pop0, sequences_pop0 = sequences_pop0, ns_pop0 = ns_pop0)
        clones, Ns_clones, Ns_cells, times, n_clones = population.simulate(T)
        # Convert states and times to a DataFrame.
        data = {}
        metadata = {}
        data2 = {}
        metadata2 = {}

        dtmin = np.min(np.diff(times))

        for n in range(n_clones-1):
            data[n] = np.array(clones[n].counts[:-1])/np.array(Ns_cells[:-1])
            metadata[n] = {'id': clones[n].id, 'sequence':[k for k in clones[n].sequence], 'fitness':clones[n].f, 'parent':clones[n].parent}
            if clones[n].parent != 0:
                metadata2[n] = {'Parent': clones[n].parent, 'Identity':clones[n].id}
        for i_t, t in enumerate(times):
            for n in range(n_clones-1):
                data2[i_t*(n_clones-1)+n] = {'Generation':int(t/dtmin), 'Identity':int(clones[n].id), 'Population': (clones[n].counts[i_t])}


        df = pd.DataFrame(data, index = times[:-1])
        metadf = pd.DataFrame(metadata)

        df2 = pd.DataFrame(data2)
        df2 = df2.T
        metadf2 = pd.DataFrame(metadata2)
        metadf2 = metadf2.T
        metadf2.index+=1
        df = df.T
        metadf = metadf.T
        # plt.plot(times, Ns)
        # plt.yscale('log')
        # plt.show()

        # # Save the DataFrame to a CSV file.
        df.to_csv(f"../../out/affinity_maturation/trajectory_{i}.csv", index = True)
        metadf.to_csv(f"../../out/affinity_maturation/metadata_{i}.csv", index = False)

        df2.to_csv(f"../../out/affinity_maturation/populations_df.csv", index = True)
        metadf2.to_csv(f"../../out/affinity_maturation/adjacency_df.csv", index = True)
    
if __name__ == "__main__":
    # To avoid issues with multiprocessing on Windows, call the main() function within this block
    main()
