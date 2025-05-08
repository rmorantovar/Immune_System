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

antigen_sequence = "".join(np.random.choice(Alphabet, L))
#antigen_sequence = 'kbbdsgesj'
#n_seq = input('How many sequences do you want to create? ')

n_seqs = np.array([5e5])

master_sequence_energies = np.array([])
range_energies = np.array([])

for n_seq in n_seqs:

	print('n_seq = %e'%(n_seq))
	start_time = time.time()

	for i in range(int(n_seq)):

		antigen_sequence = "".join(np.random.choice(Alphabet, L))
		master_sequence = find_complementary_seq(antigen_sequence, M2)
		worst_sequence = find_complementary_seq_2(antigen_sequence, M2)
		master_sequence_energy = calculate_energy(M2, antigen_sequence, master_sequence)
		worst_sequence_energy = calculate_energy(M2, worst_sequence, master_sequence)

		master_sequence_energies = np.append(master_sequence_energies, master_sequence_energy)
		range_energies = np.append(range_energies, worst_sequence_energy - master_sequence_energy)

	pickle.dump( master_sequence_energies, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/master_sequence_energies_MJ2_L-%d_n_seq-%d.pkl"%(int(L), int(n_seq)), "wb" ) )
	pickle.dump( range_energies, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/range_energies_MJ2_L-%d_n_seq-%d.pkl"%(int(L), int(n_seq)), "wb" ) )

	print("--- %.1f minutes --- \n" %((time.time() - start_time)/60))