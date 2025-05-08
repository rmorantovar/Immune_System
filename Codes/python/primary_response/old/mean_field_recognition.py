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

rhos = np.logspace(np.log10(0.08), np.log10(2), 6)
N_epitopes = np.linspace(1, 6, 20)
g=4
n0 = 2
e0 = 2.5
dt = 0.001
T=100

time = np.linspace(0, T, int((T-0)/dt))

def R_MF(r):
	return np.array([(1-r_i)*np.cumsum(R[:np.where(time<(1/(1-r_i)))[0][-1]]*dt)[-1] for r_i in r])

R_bar = np.zeros((len(N_epitopes), len(rhos)))
r = np.linspace(0.01, 0.99, 500)

for i, rho in enumerate(rhos):
	for j, g in tqdm(enumerate(N_epitopes)):
		R = 1-(1-binding_affinity(time*rho, -n0*e0, e0))**g

		R_bar[j, i] = (R_MF(r))[np.where((R_MF(r)-r)>=0)[0][-1]]

df_R = pd.DataFrame(R_bar)
df_R.to_csv('/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Complexity/R_bar_MF_continuous.txt', sep = '\t', index = False, header = False)
