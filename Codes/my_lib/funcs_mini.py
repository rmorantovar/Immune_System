# ----- imports -----
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
from matplotlib.patches import Patch
from datetime import datetime, timedelta
from tqdm import tqdm
from scipy.optimize import curve_fit
from matplotlib import style
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from itertools import product
from collections import deque
from collections import defaultdict
from ast import literal_eval
from numba import njit, prange
from scipy.integrate import solve_ivp

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
my_yellow2 = np.array((242, 192, 65))/256.
my_green2 = np.array((158,248,72))/256.
my_cyan = 'tab:cyan'

antigen_color = my_yellow

my_green_a = np.array((159, 206, 99))/256.
my_green_b = np.array((79, 173, 91))/256.
my_green_c = np.array((94, 129, 63))/256.

# ------- functions mini -------
def hamming_distance(chaine1, chaine2):

    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_complementary_seq_min(sequence, energy_model, Text_files_path):

    Alphabet = np.loadtxt(Text_files_path+'in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)
    M = np.loadtxt(Text_files_path+'in/' + energy_model + '.txt', skiprows= 0, usecols=range(0,20)).T
    #Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
    list_seq = list(sequence)
    list_comp_seq = []
    for i in list_seq:
        pos_i = np.where(np.isin(Alphabet,i))[0][0]
        list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],min(M[pos_i])))[0][0]])
    comp_seq = "".join(list_comp_seq)
    return comp_seq

def calculate_Es(motif, seqs_flat, R, l, d, L0, E_m): #For flat energies
    M = ((R==seqs_flat).T).reshape(int(L0), l, d)
    Es = np.sum(M*motif.T, axis = (1, 2)) + E_m
    return Es

def mcmc_sampling(l, energy_model, beta, num_samples, burn_in=0):
    """Run MCMC sampling to find low-energy sequences."""
    # Initialize a random sequence
    current_sequence = [np.argmin(energy_model[:, i]) for i in range(l)]
    current_energy = calculate_energy(energy_model, current_sequence)
    print(current_energy)
    samples = []
    samples.append((current_sequence.copy(), current_energy))
    # for k in range(1):
    #     # Propose a new sequence by mutating one position
    #     new_sequence = current_sequence.copy()
    #     idx = np.random.randint(0, l)
    #     new_sequence[idx] = np.random.randint(0, 20)
    #     new_energy = calculate_energy(energy_model, new_sequence)
    #     current_sequence, current_energy = new_sequence, new_energy
    # print(current_energy)

    for step in range(num_samples + burn_in):
        # Propose a new sequence by mutating one position
        new_sequence = current_sequence.copy()
        idx = np.random.randint(0, l)
        new_sequence[idx] = np.random.randint(0, 20)
        new_energy = calculate_energy(energy_model, new_sequence)

        # Metropolis acceptance criterion
        # if np.random.rand() < np.exp(-beta * (new_energy - current_energy)) and new_energy>4:
        if new_energy<18:
            current_sequence, current_energy = new_sequence, new_energy

        # Record samples after burn-in
        # if step >= burn_in:
        samples.append((current_sequence.copy(), current_energy))
    
    return samples

def my_linear_func(x, a, b):

    return a + b*x

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

def calculate_energy_repertoire(motif, seq, l, L0):
    array = motif[seq, [i for i in range(l)]*L0]
    reshaped_array = array.reshape(L0, l)
    print(reshaped_array[0])
    result = reshaped_array.sum(axis=1)
    return result

def calculate_Q0(Tmin, Tmax, n_T, E_matrix, E_ms, L):

    Ts = np.linspace(Tmin, Tmax, n_T)
    betas = 1/Ts[:-1]
    F_PWM = -Ts*np.log(Z_Motif(E_matrix, Ts))
    Es = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts)) + E_ms
    dE = np.diff(Es)
    
    S = np.cumsum(betas[:-1]*dE)

    #Omega = 20**L
    Omega = np.sum(np.exp(S)*dE)
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

def P_min_e_Q0(N, Q0, dE):

    return (N*(1-np.cumsum(Q0*dE))**(N-1)*(Q0))

def get_data(folder_path, data_type):
    return_data_type = 0 #default for returning processed data
    # if os.path.exists(folder_path+'/processed_data_'+data_type+'.pkl'):
    try:
        file = open(folder_path+'/processed_data_'+data_type+'.pkl', 'rb')
        return_data_type = 1
        # print('Data exists already and is proccesed.  Loading it ...')
        data = pickle.load(file)
        # print('Proccesed file has been read...')
        file.close()
        return data, return_data_type
    # else:
    except FileNotFoundError:
        # print(f'Processed data does not exist')
        return [], 0

def get_repertoire_properties(betas, Q0, Es, dE, N_r):
    beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
    E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
    Kd_r = np.exp(E_r)

    return beta_r, E_r, Kd_r

def get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on):
    E_pr = Es[:-1][Es[:-1]<np.log(k_pr/k_on)][-1]
    Kd_pr = np.exp(E_pr)
    beta_pr = betas[:-1][Es[:-1]<E_pr][-1]

    return beta_pr, E_pr, Kd_pr

def get_p_properties(betas, Q0, Es, dE, theta):

    beta_theta = betas[betas>theta][-1]
    E_theta = Es[betas>theta][-1]
    Kd_theta = np.exp(E_theta)

    return beta_theta, E_theta, Kd_theta

def get_clones_sizes_C(n_act, time_array, activation_times, lamB, C, dT):
    clone_sizes = np.ones((n_act, len(time_array)))
    for i_t, t in enumerate(time_array[:-1]):
        Nb = clone_sizes[:,i_t]
        N  = np.sum(Nb) - np.size(activation_times[activation_times>t])
        deltaNb = lamB*Nb*np.max([(1-(N/C)), 0])*np.heaviside(t-activation_times, 1)
        clone_sizes[:, i_t+1] = Nb + deltaNb*dT
    # clone_sizes[1:, :][clone_sizes[1:, :]==1] = 0
    clone_sizes[clone_sizes==1] = 0
        
    return clone_sizes

def get_clones_sizes_C_N0(n_act, time_array, activation_times, N0s, lamB, C, dT):
    clone_sizes = np.ones((n_act, len(time_array)))
    # clone_sizes[:, 0] = N0s
    for i_t, t in enumerate(time_array[:-1]):
        Nb = clone_sizes[:,i_t]
        N  = np.sum(Nb) - np.size(activation_times[activation_times>t])
        deltaNb = lamB*Nb*np.max([(1-(N/C)), 0])*np.heaviside(t-activation_times, N0s)
        clone_sizes[:, i_t+1] = Nb + deltaNb*dT
    # clone_sizes[1:, :][clone_sizes[1:, :]==1] = 0
    clone_sizes[clone_sizes==1] = 0
        
    return clone_sizes

def get_clones_sizes_C_new(n_act, time_array, activation_times, N0s, lamB, C, dT):
    
    t_span = (0, time_array[-1])

    N0 = N0s.copy()
    # print(N0)
    def dNdt(t, N):
        total = np.sum(N)
        # Mask entities that haven't started growing
        active_mask = t >= activation_times
        # Vectorized logistic growth only for active entities
        dN = np.zeros_like(N)
        dN[active_mask] = lamB * N[active_mask] * (1 - total / C)
        return dN

    # Solve the system
    sol = solve_ivp(dNdt, t_span, N0, t_eval=time_array, method='RK45')

    # Pre-fill the solution array with constant init_sizes until growth starts
    N_full = sol.y.copy()

    for i in range(n_act):
        t0 = activation_times[i]
        if t0 > 0:
            before_idx = sol.t < t0
            N_full[i, before_idx] = N0s[i]
    # print(N_full)
    return N_full

def apply_filter_C(clone_sizes, activation_times, energies, lim_size):
    filter_C = clone_sizes[:, -1] > lim_size
    n_C = np.sum(filter_C)
    clone_sizes_C = clone_sizes[filter_C, :]
    activation_times_C = activation_times[filter_C]
    energies_C = energies[filter_C]

    return clone_sizes_C, activation_times_C, energies_C, filter_C, n_C

def my_plot_layout(ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 24, xlabel = '', ylabel = '', title = '', x_fontsize=24, y_fontsize = 24, t_fontsize = 24, bottom = None, top = None, left = None, right = None):
    ax.tick_params(labelsize = ticks_labelsize)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlim(left = left, right = right)
    ax.set_ylim(bottom = bottom, top = top)
    ax.set_xlabel(xlabel, fontsize = x_fontsize)
    ax.set_ylabel(ylabel, fontsize = y_fontsize)
    ax.set_title(title, fontsize = t_fontsize)







