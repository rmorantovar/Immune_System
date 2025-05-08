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
    sequences_pop0 = []
    if pop0:
        pop0_df = pd.read_csv('../../in/populations/filtered_sequence_properties_2.csv')
        pop0_df['sequence'] = pop0_df['sequence'].apply(ast.literal_eval)
        sequences_pop0 = pop0_df['sequence'].tolist()
        ns_pop0 = pop0_df['n'].tolist()
        N_clones = len(sequences_pop0)
        N_cells = np.sum(ns_pop0)

    energy_model = 'TCRen'
    antigen = 'TACNSEYPNTTRAKCGRWYC'

    # Number of trajectories in the ensemble.
    N_ens = args.N_ens

    mean_e_s = {d:[] for d in range(steps)}
    Zs = {d:[] for d in range(steps)}

    for _ in range(N_ens):
        antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
        l = len(antigen)    
        motif = get_motif(antigen_seq, energy_model, '../../')
        #Change values by the minimum
        for i in np.arange(l):
            motif[:,i]-=np.min(motif[:,i], axis=0)
        # Create and simulate a population.
        population = Population(N_clones, N_cells, l, b, d, m, motif, pop0 = pop0, sequences_pop0 = sequences_pop0, ns_pop0 = ns_pop0)
        mean_e_s[0].append(np.exp(np.mean([population.clones[k].e for k in range(len(population.clones)) for _ in range(population.clones[k].counts[-1])])))
        Zs[0].append(np.sum([population.clones[k].counts[-1]/np.exp(population.clones[k].e) for k in range(len(population.clones))]))
        for j in range(1, steps):
            antigen_seq[np.random.randint(0, l)] = np.random.randint(0, 20)
            motif = get_motif(antigen_seq, energy_model, '../../')
            #Change values by the minimum
            for i in np.arange(l):
                motif[:,i]-=np.min(motif[:,i], axis=0)
            population.motif = motif
            population.update_energies()
            zs = [population.clones[k].counts[-1]/np.exp(population.clones[k].e) for k in range(len(population.clones))]
            energies = [population.clones[k].e for k in range(len(population.clones)) for _ in range(population.clones[k].counts[-1])]
            mean_energy = np.mean(energies)
            mean_e_s[j].append(np.exp(mean_energy))
            Z = np.sum(zs)
            Zs[j].append(Z)

    bins = np.logspace(-7, 2, 20)
    q = int(np.sqrt(steps))
    fig, ax = plt.subplots(q, q, figsize=(4*3*1.6, 4*3))
    means = []
    means2 = []
    stds = []
    stds2 = []
    for j in range(steps):
        ax[j//q, j%q].hist(mean_e_s[j], bins = bins, alpha = .6, label =  r'$%d$'%j, color = my_blue)
        ax[j//q, j%q].vlines(np.exp(np.mean(np.log(mean_e_s[j]))), 1, 100, color = my_red)
        means.append(np.exp(np.mean(np.log(mean_e_s[j]))))
        means2.append(np.exp(np.mean(np.log(Zs[j]))))
        stds.append(np.exp(np.std(np.log(mean_e_s[j]))))
        stds2.append(np.exp(np.std(np.log(Zs[j]))))
        ax[j//q, j%q].set_yscale('log')
        ax[j//q, j%q].set_xscale('log')
        ax[j//q, j%q].set_ylim(1, N_ens*2)
        ax[j//q, j%q].legend()
    fig.savefig('../../../Figures/coevolution/hist.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(5*1.6, 5))
    ax.plot(range(steps), means, ls = '', marker = 'x', color = my_red)
    ax.fill_between(range(steps), np.array(means)*np.array(stds), np.array(means)/np.array(stds), alpha = .2, color = my_blue)
    ax.set_yscale('log')
    fig.savefig('../../../Figures/coevolution/mean_Kd.png')

    fig, ax = plt.subplots(figsize=(5*1.6, 5))
    ax.plot(range(steps), means2, ls = '', marker = 'x', color = my_green)
    ax.fill_between(range(steps), np.array(means2)*np.array(stds2), np.array(means2)/np.array(stds2), alpha = .2, color = my_green)
    ax.set_yscale('log')
    fig.savefig('../../../Figures/coevolution/mean_Z.png')
    
if __name__ == "__main__":
    # To avoid issues with multiprocessing on Windows, call the main() function within this block
    main()
