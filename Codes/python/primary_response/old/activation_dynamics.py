# Code to model the activation of BCCs as the antigen grows

import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import pickle
from tqdm import tqdm

L = 1000
V = L*L
N_R = 100

#Lattice for BCCs and antigens
lattice = np.ones(shape=(L,L))

#produce N_R random BBCs and locate them in the lattice
