import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as scipy
import pickle
import time
from tqdm import tqdm


a = -6
b = 3.0
N_epitopes = np.arange(1, 12, 1)
rho = 2.0
tau = 4
T = 2000*tau
R_bar = np.array([])
for g in N_epitopes:
    n_mutations = np.zeros(g) #Array with number of mutatios in each epitope
    neutralizations = binding_affinity(n_mutations, a, b) #Array with individuals neutralization probabilities
    time = np.array([0]) #Array with time
    t = time[-1]
    time_intervals = np.array([])
    exposure_times = np.array([0])
    sick_times = np.array([])
    R = np.array([1-np.min(1-neutralizations)]) #Array with total recognition function
    propensities = np.array(np.concatenate(([rho]*g, [1/tau]))) #Array with propensities
    i = 0
    while (t<T):
        i+=1
        cumsum = np.cumsum(propensities)
        alpha = np.sum(propensities)
        r1 = np.random.rand()
        dti = (1/alpha)*np.log(float(1/r1))
        time_intervals = np.append(time_intervals, dti)
        t = t+dti
        time = np.append(time, t)
        r2 = np.random.rand()
        idEvent = np.searchsorted(cumsum,r2*alpha)
        if(idEvent<g):#mutation
            n_mutations[idEvent]+=1
        else:#exposure
            exposure_times = np.append(exposure_times, t)
            r_sick = np.random.rand()
            if(r_sick>R[-1]): #sick
                n_mutations = np.zeros(g) #update mutations
                sick_times = np.append(sick_times, t)
        neutralizations = binding_affinity(n_mutations, a, b) #update neutralizations  
        R = np.append(R, 1-np.product(1-neutralizations)) #update R
    
    #R_bar = np.append(R_bar, np.sum(R[:-1]*time_intervals)/time[-1])
    R_bar = np.append(R_bar, 1-np.size(sick_times)/np.size(exposure_times))

df_R = pd.DataFrame(list(zip(N_epitopes, R_bar)), header = False)
df_R.to_csv('/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Complexity/R_fast.txt', sep = '\t', index = False)

