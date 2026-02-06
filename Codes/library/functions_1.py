import struct
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy
import seaborn as sns
import scipy.special as sc
import time
from joblib import Parallel, delayed
import argparse
from scipy.integrate import odeint
import json
import warnings
from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
from tqdm import tqdm
from scipy.optimize import curve_fit
from matplotlib import style
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from itertools import product
from collections import deque
from collections import defaultdict
if(sys.version_info[1]<= 7):
    import pickle5 as pickle
else:
    import pickle
plt.rcParams['text.usetex'] = True

# ----- CONSTANTS ------
N_A = 6.02214076e23
k_BT = 1.380649e-23*293

my_red = np.array((228,75,41))/256.
my_purple = np.array((125,64,119))/256.
my_purple2 = np.array((116,97,164))/256.
my_green = np.array((125,165,38))/256.
my_blue = np.array((76,109,166))/256.
my_gold = np.array((215,139,45))/256.
my_brown = np.array((182,90,36))/256.
my_blue2 = np.array((80,141,188))/256.
my_yellow = np.array((246,181,56))/256.
my_green2 = np.array((158,248,72))/256.
my_cyan = 'tab:cyan'

# ------- functions -------

def Z_Motif(Motif, T):
    Z = 1
    for i in range(len(Motif[0,:])):
        Z_i = 0
        for j in range(len(Motif[:,0])):
            Z_i = Z_i + np.exp((-Motif[j, i]/T))
        Z = Z*Z_i
    return Z

def get_motif(antigen_seq, energy_model, Text_files_path):

    M = np.loadtxt(Text_files_path+'in/' + energy_model + '.txt', skiprows= 0, usecols=range(0,20))
    return M[:,antigen_seq]

def from_aa_to_i(aa_seq, energy_model, Text_files_path):
    Alphabet = np.loadtxt(Text_files_path+'in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    aa_list = [i for i in aa_seq]
    seq = []
    for i, aa in enumerate(aa_list):
        index = Alphabet_list.index(aa)
        seq.append(int(index))
    return seq

def from_aa_to_i_Alphabet(Alphabet, aa_seq):
    Alphabet_list = Alphabet.tolist()
    aa_list = [i for i in aa_seq]
    seq = []
    for i, aa in enumerate(aa_list):
        index = Alphabet_list.index(aa)
        seq.append(int(index))
    return seq

def from_aa_to_i_X(aa_seq, energy_model, Text_files_path):
    Alphabet = np.loadtxt(Text_files_path+'in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    aa_list = [i for i in aa_seq]
    seq = []
    for aa in (aa_list):
        if aa in ' +()-ZXBOUmpazgunyqk;tolfxhiebrdc0123456789/,':
            index = np.random.randint(0, 20)
        else:
            index = Alphabet_list.index(aa)
        seq.append(int(index))
    return seq

def from_i_to_aa(seq, energy_model, Text_files_path):
    Alphabet = np.loadtxt(Text_files_path+'in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    aa_seq = ''.join([Alphabet_list[i] for i in seq])
    return aa_seq

def from_i_to_aa_Alphabet(Alphabet, seq):
    Alphabet_list = Alphabet.tolist()
    aa_seq = ''.join([Alphabet_list[i] for i in seq])
    return aa_seq

def calculate_energy(motif, seq):
    return np.sum(motif[seq, np.arange(int(len(seq)))])

def calculate_Q0(Tmin, Tmax, n_T, Motif, E_m, L):

    Ts = np.linspace(Tmin, Tmax, n_T)
    betas = 1/Ts[:-1]

    F_PWM = -Ts*np.log(Z_Motif(Motif, Ts))
    S = -(np.diff(F_PWM)/np.diff(Ts))
    Es = F_PWM[:-1]+Ts[:-1]*S + E_m
    dE = np.diff(Es)

    #S2 = np.log(Z_Motif(Motif, Ts[:-2])) + betas[:-1]*(Es[:-1]-E_m)

    #Omega = 20**L
    #Omega1 = np.sum(np.exp(S1)*dE)
    Omega= np.sum(np.exp(S[:-1])*dE)
    #Omega = 2*np.sum(np.exp(S)*dE)

    #print('%.2e'%(20**L), '%.2e'%(np.sum(np.exp(S)*dE)))
    Q0 = np.exp(S)/Omega

    return Es, dE, Q0, betas

def calculate_QR(Q0, k_on, k_act, rho_A, Es, q, lambda_A, N_c, dE):

    p_a = (1/(1 + (k_on*np.exp(Es[:-1])/k_act))**q)
    u_on = rho_A*k_on
    R = 1-np.exp(-u_on*p_a*N_c/lambda_A)
    QR = Q0*R

    return u_on, p_a, R, QR

def get_clones_sizes_C(n_act, time_array, activation_times, lambda_B, C, dT):
    clone_sizes = np.ones((n_act, len(time_array)))
    for i_t, t in enumerate(time_array[:-1]):
        Nb = clone_sizes[:,i_t]
        N  = np.sum(Nb) - np.size(activation_times[activation_times>t])
        deltaNb = lambda_B*Nb*(1-(N/C))*np.heaviside(t-activation_times, 1)
        clone_sizes[:, i_t+1] = Nb + deltaNb*dT
        
    return clone_sizes

def apply_filter_C(clone_sizes, activation_times, energies, lim_size):
    filter_C = clone_sizes[:, -1] > lim_size
    n_C = np.sum(filter_C)
    clone_sizes_C = clone_sizes[filter_C, :]
    activation_times_C = activation_times[filter_C]
    energies_C = energies[filter_C]

    return clone_sizes_C, activation_times_C, energies_C, filter_C, n_C

def my_plot_layout(ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 24, xlabel = '', ylabel = '', title = '', x_fontsize=24, y_fontsize = 24, t_fontsize = 24):
    ax.tick_params(labelsize = ticks_labelsize)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize = x_fontsize)
    ax.set_ylabel(ylabel, fontsize = y_fontsize)
    ax.set_title(title, fontsize = t_fontsize)

# Function to generate random sequences and compute properties
def generate_sequences_and_properties(motif, cum_Omega_0, Es_avg, ensemble_id, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, chunk_size, seqs = False):
    T0 = 0
    Tf = 10
    Tf_sim = 7
    dT = 0.05
    lambda_A = lamA
    k_on = 1e6*24*3600; #(M*days)^-1
    b0 = 1e5
    N0=1
    N_A = 6.02*10**23
    times = np.linspace(0, Tf, 1000)
    properties = []
    for _ in range(L0 // chunk_size):
        # Sequences can be real integer sequences or random numbers between 0 and 1. In the second case,
        # they are the random numbers to be used to sample the random energies.
        if not seqs:
            sequences = np.random.rand(chunk_size)
            proto_Es = np.searchsorted(cum_Omega_0,sequences)-1
        else:
            sequences = np.random.randint(0, 20, size=(chunk_size, l))
        for proto_E in proto_Es:
            if not seqs:
                E = Es_avg[proto_E]
                if E < E_lim:
                    F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                    # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                    r1 = np.random.random()
                    #t1 = times[F1<r1][-1]
                    t1 = times[np.searchsorted(F1,r1)-1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'time': t1
                        })
            else:
                E = calculate_energy(motif, proto_E) + E_m
                if E < E_lim:
                    F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                    for n0 in range(int(N0)):
                        r1 = np.random.random()
                        #t1 = times[F1<r1][-1]
                        t1 = times[np.searchsorted(F1,r1)-1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'time': t1, 
                            'sequence': seq.tolist()
                        })
                
    return properties

def generate_sequences_and_properties_2(motif, cum_Omega_0, Es_avg, ensemble_id, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, chunk_size, dir_memory, seqs = False):
    T0 = 0
    Tf = 10
    Tf_sim = 7
    dT = 0.05
    lambda_A = lamA
    k_on = 1e6*24*3600; #(M*days)^-1
    b0 = 1e5
    N0=1
    N_A = 6.02*10**23
    times = np.linspace(0, Tf, 1000)
    properties = []
    for _ in range(L0 // chunk_size):
        # Sequences can be real integer sequences or random numbers between 0 and 1. In the second case,
        # they are the random numbers to be used to sample the random energies.
        if not seqs:
            sequences = np.random.rand(chunk_size)
            proto_Es = np.searchsorted(cum_Omega_0,sequences)-1
        else:
            sequences = np.random.randint(0, 20, size=(chunk_size, l))
        for proto_E in proto_Es:
            if not seqs:
                E = Es_avg[proto_E]
                if E < E_lim:
                    F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                    # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                    r1 = np.random.random()
                    #t1 = times[F1<r1][-1]
                    t1 = times[np.searchsorted(F1,r1)-1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'time': t1
                        })
            else:
                E = calculate_energy(motif, proto_E) + E_m
                if E < E_lim:
                    F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                    for n0 in range(int(N0)):
                        r1 = np.random.random()
                        #t1 = times[F1<r1][-1]
                        t1 = times[np.searchsorted(F1,r1)-1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'time': t1, 
                            'sequence': seq.tolist()
                        })
    
    memory_clones = pd.read_csv(dir_memory)
    memory_clones = memory_clones.loc[memory_clones['ens_id'] == ensemble_id]
    if len(memory_clones.index)/2 >0:
        for index, clone in memory_clones.sample(n=int(len(memory_clones.index)/2), replace=False).iterrows():
            E = float(clone['E'])
            for cell in range(int(clone['n'])):
                if E < E_lim:
                    F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                    # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                    r1 = np.random.random()
                    #t1 = times[F1<r1][-1]
                    t1 = times[np.searchsorted(F1,r1)-1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'time': t1
                        })

    return properties

def process_ensemble(motif, cum_Omega_0, Es_avg, N_ens, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, chunk_size, output_csv_file, n_jobs=-1, seqs = False):

    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        delayed(generate_sequences_and_properties)(motif, cum_Omega_0, Es_avg, i, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, chunk_size) for i in range(N_ens)
    )
    all_properties = [prop for sublist in results for prop in sublist]
    df = pd.DataFrame(all_properties)
    df.to_csv(output_csv_file, index=False)

def process_ensemble_2(motif, cum_Omega_0, Es_avg, N_ens, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, chunk_size, input_memory_file, output_csv_file, n_jobs=-1, seqs = False):
    
    # Create a progress bar using tqdm
    # progress_bar = tqdm(total=N_ens, desc="Processing", unit="tasks")

    # Define a callback function to update the progress bar
    # def update_progress_bar(result):
    #     progress_bar.update(n=1)

    # Wrap the generate_sequences_and_properties function to include the progress bar update
    # def wrapped_generate_sequences_and_properties(*args, **kwargs):
    #     result = generate_sequences_and_properties(*args, **kwargs)
    #     update_progress_bar()
    #     return result

    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        #delayed(generate_sequences_and_properties)(PWM, i, L0, l, E_lim, E_m, p, k_step, chunk_size) for i in range(N_ens)
        delayed(generate_sequences_and_properties_2)(motif, cum_Omega_0, Es_avg, i, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, chunk_size, input_memory_file) for i in range(N_ens)
    )
    # progress_bar.close()
    all_properties = [prop for sublist in results for prop in sublist]
    df = pd.DataFrame(all_properties)
    df.to_csv(output_csv_file, index=False)

def generate_repertoire(Alphabet, motif, cum_Omega_0, Es_avg, ensemble_id, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, infection, chunk_size, input_memory_file, seqs = False):
    T0 = 0
    Tf = 10
    Tf_sim = 7
    dT = 0.05
    lambda_A = lamA
    k_on = 1e6*24*3600; #(M*days)^-1
    b0 = 1e5
    N0=1
    times = np.linspace(0, Tf, 1000)
    properties = []
    for _ in range(L0 // chunk_size):
        # Sequences can be real integer sequences or random numbers between 0 and 1. In the second case,
        # they are the random numbers to be used to sample the random energies.
        if not seqs:
            random_numbers = np.random.rand(chunk_size)
            proto_Es = np.searchsorted(cum_Omega_0,random_numbers)-1
        else:
            proto_Es = np.random.randint(0, 20, size=(chunk_size, l))
            # proto_Es = np.array([calculate_energy(motif, seq) for seq in sequences])

        for proto_E in proto_Es:
            if not seqs:
                E = Es_avg[proto_E]
                if E < E_lim:
                    F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                    # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                    r1 = np.random.random()
                    #t1 = times[F1<r1][-1]
                    t1 = times[np.searchsorted(F1,r1)-1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'time': t1,
                            'm' : 0
                        })
            else:
                E = calculate_energy(motif, proto_E) + E_m
                if E < E_lim:
                    F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                    for n0 in range(int(N0)):
                        r1 = np.random.random()
                        #t1 = times[F1<r1][-1]
                        t1 = times[np.searchsorted(F1,r1)-1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'time': t1, 
                            'sequence': from_i_to_aa_Alphabet(Alphabet, proto_E),
                            'm' : 0
                        })

    if infection > 1:
        memory_clones = pd.read_csv(input_memory_file)
        memory_clones = memory_clones.loc[memory_clones['ens_id'] == ensemble_id]
        if len(memory_clones.index) >0:
            if not seqs:
                #for index, clone in memory_clones.sample(n=int(len(memory_clones.index)/2), replace=False).iterrows():
                for E in memory_clones['E'].sample(n = 100, replace = True, weights = np.array(memory_clones['n'])):
                    #E = float(clone['E'])
                    #for cell in range(int(clone['n'])):
                    if E < E_lim:
                        F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                        # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                        r1 = np.random.random()
                        #t1 = times[F1<r1][-1]
                        t1 = times[np.searchsorted(F1,r1)-1]
                        if t1 < t_lim:
                            properties.append({
                                'ens_id': ensemble_id,
                                'E': E,
                                'time': t1,
                                'm' : 1
                            })
            else:
                memory = memory_clones['sequence'].sample(n = 100, replace = True, weights = np.array(memory_clones['n']))
                for seq_aa in memory:
                    seq = from_aa_to_i_Alphabet(Alphabet, seq_aa)
                    E = calculate_energy(motif, np.array(seq)) + E_m
                    #for cell in range(int(clone['n'])):
                    if E < E_lim:
                        F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                        # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                        r1 = np.random.random()
                        #t1 = times[F1<r1][-1]
                        t1 = times[np.searchsorted(F1,r1)-1]
                        if t1 < t_lim:
                            properties.append({
                                'ens_id': ensemble_id,
                                'E': E,
                                'time': t1,
                                'sequence': seq_aa,
                                'm' : 1
                            })


    return properties

def generate_repertoire_seqs(Alphabet, motif, cum_Omega_0, Es_avg, ensemble_id, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, infection, chunk_size, input_memory_file, seqs = False):
    T0 = 0
    Tf = 10
    Tf_sim = 7
    dT = 0.05
    lambda_A = lamA
    k_on = 1e6*24*3600; #(M*days)^-1
    b0 = 1e5
    N0=1
    N_A = 6.02*10**23
    times = np.linspace(0, Tf, 1000)
    properties = []
    for _ in range(L0 // chunk_size):
        proto_Es = np.random.randint(0, 20, size=(chunk_size, l))
        for proto_E in proto_Es:
            E = calculate_energy(motif, proto_E) + E_m
            # print(E)
            if E < E_lim:
                F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                r1 = np.random.random()
                #t1 = times[F1<r1][-1]
                t1 = times[np.searchsorted(F1,r1)-1]
                # print(t1)
                if t1 < t_lim:
                    properties.append({
                        'ens_id': ensemble_id,
                        'E': E,
                        'time': t1, 
                        'sequence': from_i_to_aa_Alphabet(Alphabet, proto_E),
                        'm' : 0
                    })

    if infection > 1:
        memory_clones = pd.read_csv(input_memory_file)
        memory_clones = memory_clones.loc[memory_clones['ens_id'] == ensemble_id]
        if len(memory_clones.index) >0:
            memory = memory_clones['sequence'].sample(n = 1000, replace = True, weights = np.array(memory_clones['n']))
            for seq_aa in memory:
                seq = from_aa_to_i_Alphabet(Alphabet, seq_aa)
                E = calculate_energy(motif, np.array(seq)) + E_m
                #for cell in range(int(clone['n'])):
                if E < E_lim:
                    F1 = 1-np.exp(-(b0*10*k_on)/(lambda_A*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lambda_A*times)-1))
                    # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                    r1 = np.random.random()
                    #t1 = times[F1<r1][-1]
                    t1 = times[np.searchsorted(F1,r1)-1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'time': t1,
                            'sequence': seq_aa,
                            'm' : 1
                        })


    return properties

def ensemble_of_activations(Alphabet, motif, cum_Omega_0, Es_avg, N_ens, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, infection, chunk_size, input_memory_file, n_jobs=-1, seqs = False):

    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        delayed(generate_repertoire_seqs)(Alphabet, motif, cum_Omega_0, Es_avg, i, L0, l, t_lim, E_lim, E_m, p, k_step, lamA, infection, chunk_size, input_memory_file, seqs = seqs) for i in range(N_ens)
    )
    all_properties = [prop for sublist in results for prop in sublist]
    df = pd.DataFrame(all_properties)
    return df

def ensemble_of_expansions(data, N_ens, p, time_array, lambda_B, C, dT):
    data = data
    clone_size_total = []
    ids_C_total = []
    for i_ens in np.arange(N_ens):
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
        ids_C_total = np.concatenate((ids_C_total, ids_C))
    
    final_clone_size = []
    i = 0
    for j in range(len(data)):
        if j in ids_C_total:
            final_clone_size.append(int(clone_size_total[i]))
            i+=1
        else:
            final_clone_size.append(0)

    data['n'] = final_clone_size
    data = data.loc[data['n']!=0]
    return data

def ensemble_of_expansions_time(data, N_ens, p, time_array, lambda_B, C, dT):

    clone_size_total = []
    ids_C_total = []
    data_active = data
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

    data['n_t'] = [clone_sizes_C[i] for i in range(len(ids_C))]

    return data
















