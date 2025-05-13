import struct
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.special as sc
import time
from scipy.integrate import odeint

if(sys.version_info[1]<= 7):
    import pickle5 as pickle
else:
    import pickle

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


def my_plot_layout(ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 24, xlabel = '', ylabel = '', title = '', x_fontsize=24, y_fontsize = 24, t_fontsize = 24):
    ax.tick_params(labelsize = ticks_labelsize)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize = x_fontsize)
    ax.set_ylabel(ylabel, fontsize = y_fontsize)
    ax.set_title(title, fontsize = t_fontsize)

def energy(E_matrix, W_matrix, peptide, tcr, l):

    #motif = E_matrix[:,peptide]
    E = 0
    for i in range(l):
        # E += np.sum(motif[tcr[i]]*W_matrix[i])
        for j in range(l):
            E += E_matrix[tcr[j], peptide[i]]*W_matrix[i, j]
    return E

def Z_PWM(PWM, T):
    Z = 1
    for i in range(len(PWM[0,:])):
        Z_i = 0
        for j in range(len(PWM[:,0])):
            Z_i = Z_i + np.exp((-PWM[j, i]/T))
        Z = Z*Z_i
    return Z

def plot_energy_matrix(Energy_Matrix, Alphabet, title, ax):

    M = Energy_Matrix
    
    Alphabet = Alphabet

    sns.heatmap(np.flip(M, axis = 0), ax = ax, cmap=plt.cm.seismic, center = 0, cbar = True)
    ax.set_title(title, fontsize = 22)
    ax.tick_params(labelsize = 20)
    ax.set_xticklabels(Alphabet)
    ax.set_yticklabels(np.flip(Alphabet));
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=18)

def calculate_Q0(Tmin, Tmax, n_T, E_matrix, E_ms, L):

    Ts = np.linspace(Tmin, Tmax, n_T)
    betas = 1/Ts[:-1]
    F_PWM = -Ts*np.log(Z_PWM(E_matrix, Ts))
    Es = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts)) + E_ms
    dE = np.diff(Es)
    
    # S1 = np.cumsum(betas[:-1]*dE)
    # S2 = np.log(Z_PWM(E_matrix, Ts[:-2])) + betas[:-1]*(Es[:-1]-E_ms)
    S2 = -(np.diff(F_PWM)/np.diff(Ts))
    #Omega = 20**L
    # Omega1 = np.sum(np.exp(S1)*dE)
    Omega2= np.sum(np.exp(S2[:-1])*dE)
    #Omega = 2*np.sum(np.exp(S)*dE)
    
    #print('%.2e'%(20**L), '%.2e'%(np.sum(np.exp(S)*dE)))
    # Q01 = np.exp(S1)/Omega1
    Q02 = np.exp(S2)/Omega2

    return Es, dE, Q02, betas


def my_linear_func(x, a, b):

    return a + b*x