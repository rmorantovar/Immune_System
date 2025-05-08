import sys
sys.path.append('../lib/')
from functions_1 import*
from classes import*
import numpy as np
import pickle
from tqdm import tqdm
import pandas as pd
import math
import matplotlib.pyplot as plt
import ast

data = pd.read_csv('iedb_slim.csv')
data = data.loc[data['MHC.3']=='I']
print(np.unique(data['Epitope.2']))



potential = pd.read_csv('TCRen_potential.csv')
Alphabet_T = np.unique(potential['residue.aa.from'])
Alphabet_p = np.unique(potential['residue.aa.to'])

#
potential_dict = {}
for aa_p in Alphabet_p:
	#print(aa_p)
	potential_dict[aa_p] = {}
	A = potential[potential['residue.aa.to']==aa_p]
	for aa_T in Alphabet_T:
		#print(aa_T)
		B = A[A['residue.aa.from']==aa_T]
		potential_dict[aa_p][aa_T] = float(B['TCRen'])

potential_df = pd.DataFrame(potential_dict)
potential = np.array(potential_df)

fig, ax = plt.subplots()
sns.heatmap(potential_df, ax = ax)
fig.savefig('../../Figures/TCRen.pdf')
