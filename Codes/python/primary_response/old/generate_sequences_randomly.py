import sys
import numpy as np
import matplotlib.pyplot as plt
from Immuno_models import*
#from Bio import Phylo
from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
import scipy.special as sc
import os.path
import pickle
from matplotlib import style
from scipy.interpolate import interp1d
import time


M2 = np.loadtxt('../../../../Dropbox/Research/Evolution_Immune_System/Text_files/MJ2.txt', skiprows= 1, usecols=range(1,21)).tolist()
Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
L = int(input('What is the length of the sequences? '))

#antigen_sequence = "".join(np.random.choice(Alphabet, L))
antigen_sequence = 'kbbdsgesj'
#n_seq = input('How many sequences do you want to create? ')

n_seqs = np.array([1e4, 1e5, 1e6])

for n_seq in n_seqs:
	start_time = time.time()
	
	Sequences = generate_Sequences_randomly(int(n_seq), Energy_Matrix = M2, antigen_sequence = antigen_sequence, L = L, new_antigen = False)
	pickle.dump( Sequences, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Sequences_random_MJ2_L-%d_n_seq-%d.pkl"%(int(L), int(n_seq)), "wb" ) )

	print("--- %.1f minutes --- \n" %((time.time() - start_time)/60))