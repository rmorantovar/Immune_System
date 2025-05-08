import sys
sys.path.append('../../lib/')
from functions_1 import*
from classes import*
#from functions_2 import*
import numpy as np
import pickle
from tqdm import tqdm
import pandas as pd
import math
import matplotlib.pyplot as plt
import ast
from scipy.spatial import distance_matrix

def main():
    parser = argparse.ArgumentParser(description="Simulate affinity maturation")
    parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
    parser.add_argument('--N', type=int, default=20, help="Number of random sequences.")
    parser.add_argument('--l', type=int, default=20, help="Length of the sequences.")
    parser.add_argument('--b', type=float, default=1, help="Birth rate.")
    parser.add_argument('--d', type=float, default=.1, help="Death rate.")
    parser.add_argument('--m', type=float, default=0.01, help="Mutation rate.")
    parser.add_argument('--T', type=float, default=10, help="Simulation time.")
    parser.add_argument('--steps', type=int, default=100, help="number of antigenic steps")
    parser.add_argument('--pop0', type=bool, default=False, help="add initial pop?")
    
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
    steps = args.steps
    pop0 = args.pop0

    if pop0:
        pop0_df = pd.read_csv('../../in/populations/filtered_sequence_properties_2.csv')
        pop0_df['sequence'] = pop0_df['sequence'].apply(ast.literal_eval)
        sequences_pop0 = pop0_df['sequence'].tolist()
        ns_pop0 = pop0_df['n'].tolist()
        N_clones = len(sequences_pop0)
        N_cells = np.sum(ns_pop0)

    energy_model = 'TCRen'
    antigens = ['TACNSEYPNTTRAKCGRWYC', 'TACNSEYPNTTFDKCGRWYC', 'MACNSEYPNTTRAKCGRWYC', 'MACNSEYPNTTRCKCLRWYC', 'YACNSEYPNTTFDKCGRWYC', 'TACNSTYPNTERAKCGRWYC', 'MACNSEYPNTTRCRKLRWYC',
    'TACNSKYPNDTTDKCGRWSC', 'MACNSECPNTTRCKWLRWYC', 'MACNSEYPNTTRCKEFRWYC'] #L=20

    # Number of trajectories in the ensemble.
    N_ens = args.N_ens

    letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']
    mean_e_s = {letters[d]:[] for d in range(len(antigens))}
    Zs = {letters[d]:[] for d in range(len(antigens))}
    print(Zs)
    for a, serum in enumerate(antigens):
        pop0_df = pd.read_csv('../../in/populations/activated_population_' + serum + '_2.csv')
        pop0_df['sequence'] = pop0_df['sequence'].apply(ast.literal_eval)
        sequences_pop0 = pop0_df['sequence'].tolist()
        ns_pop0 = pop0_df['n'].tolist()
        N_clones = len(sequences_pop0)
        N_cells = np.sum(ns_pop0)
        serum_seq = from_aa_to_i(serum, energy_model, '../../')
        l = len(serum)    
        motif = get_motif(serum_seq, energy_model, '../../')
        #Change values by the minimum
        for i in np.arange(l):
            motif[:,i]-=np.min(motif[:,i], axis=0)
        # Create and simulate a population.
        population = Population(N_clones, N_cells, l, b, d, m, motif, pop0 = pop0, sequences_pop0 = sequences_pop0, ns_pop0 = ns_pop0)

        for b, antigen in enumerate(antigens):
            antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
            motif = get_motif(antigen_seq, energy_model, '../../')
            #Change values by the minimum
            for i in np.arange(l):
                motif[:,i]-=np.min(motif[:,i], axis=0)
            population.motif = motif
            population.update_energies()
            energies = [population.clones[k].e for k in range(len(population.clones)) for _ in range(population.clones[k].counts[-1])]
            mean_energy = np.mean(energies)
            mean_e_s[letters[a]].append(np.exp(mean_energy))
            zs = [population.clones[k].counts[-1]/np.exp(population.clones[k].e) for k in range(len(population.clones))]
            Z = np.sum(zs)
            Zs[letters[a]].append(np.log10(Z))
    # for a in range(len(antigens)):
    for a in range(len(antigens)):
        homologous_Z = Zs[letters[a]][a]
        for b in range(len(antigens)):
            Zs[letters[a]][b] = Zs[letters[a]][b] - homologous_Z + 11
        
    Zs_df = pd.DataFrame(Zs)
    print(Zs_df)
    Zs_matrix = Zs_df.to_numpy()
    fig, ax = plt.subplots(figsize=(5*1.6, 5), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.96})
    sns.heatmap(Zs_df, ax = ax)
    ax.set_xlabel(r'Serum', fontsize = 22)
    ax.set_ylabel(r'Virus', fontsize = 22)
    fig.savefig('../../../Figures/coevolution/panel.png')

    fig, ax = plt.subplots(figsize=(5*1.6, 5), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.96})
    sns.heatmap(distance_matrix(Zs_matrix, Zs_matrix), ax = ax)
    ax.set_xlabel(r'Virus', fontsize = 22)
    ax.set_ylabel(r'Virus', fontsize = 22)
    fig.savefig('../../../Figures/coevolution/distance1.png')

    fig, ax = plt.subplots(figsize=(5*1.6, 5), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.96})
    sns.heatmap(distance_matrix(Zs_matrix.T, Zs_matrix.T), ax = ax)
    ax.set_xlabel(r'Serum', fontsize = 22)
    ax.set_ylabel(r'Serum', fontsize = 22)
    fig.savefig('../../../Figures/coevolution/distance2.png')
    
if __name__ == "__main__":
    # To avoid issues with multiprocessing on Windows, call the main() function within this block
    main()
