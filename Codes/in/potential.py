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

fig_total, ax_total  = plt.subplots(1, 3, figsize = (30, 8))

#------- MJ ------
eir = [-3.57, -3.92, -4.76, -4.42, -4.81, -3.89, -3.81, -3.41,
-2.57, -2.19, -2.29, -1.98, -1.92, -2.00, -1.84, -1.79, -2.56, -2.11,
-1.52, -2.09]

potential_df = pd.read_csv('MJ.txt', header = 0, index_col = 0, delimiter = '\t')
alphabet_MJ = potential_df.columns.values.tolist()


for i, r in enumerate(potential_df.iterrows()):
	index = r[0]
	row = r[1]
	row.iloc[0:i] = potential_df[index].iloc[0:i]

for i, r in enumerate(potential_df.iterrows()):
	index = r[0]
	row = r[1]
	row.iloc[0:] = 3*row.iloc[0:] - 3*np.array(eir)

fig, ax = plt.subplots()
sns.heatmap(potential_df, ax = ax)
sns.heatmap(potential_df, ax = ax_total[0])
fig.savefig('../../Figures/MJ.pdf')

#------- BLOSUM62 ------
potential_df = pd.read_csv('BLOSUM62.txt', delimiter = '\t', header = 0, index_col = 0)
alphabet_BM = potential_df.columns.values.tolist()
potential_df.drop('*', axis=1, inplace = True)
potential_df.drop('*', axis=0, inplace = True)

potential_df = potential_df[alphabet_MJ]
potential_df = potential_df.reindex(alphabet_MJ)

fig, ax = plt.subplots()
sns.heatmap(potential_df, ax = ax)
sns.heatmap(potential_df, ax = ax_total[1])
fig.savefig('../../Figures/BLOSUM62.pdf')

#------- TCRen ------
potential_csv = pd.read_csv('TCRen_potential.csv')
Alphabet_T = np.unique(potential_csv['residue.aa.from'])
Alphabet_p = np.unique(potential_csv['residue.aa.to'])
#print(Alphabet_T, Alphabet_p)

potential_csv = pd.read_csv('TCRen_potential.csv')
Alphabet_T = np.unique(potential_csv['residue.aa.from'])
Alphabet_p = np.unique(potential_csv['residue.aa.to'])
#print(Alphabet_T, Alphabet_p)

potential_dict = {}
for aa_p in Alphabet_p:
	#print(aa_p)
	potential_dict[aa_p] = {}
	A = potential_csv[potential_csv['residue.aa.to']==aa_p]
	for aa_T in Alphabet_T:
		#print(aa_T)
		B = A[A['residue.aa.from']==aa_T]
		potential_dict[aa_p][aa_T] = float(B['TCRen'])

potential_df = pd.DataFrame(potential_dict)
potential_df = potential_df[alphabet_MJ]
alphabet_MJ.remove('C')
potential_df = potential_df.reindex(alphabet_MJ)

fig, ax = plt.subplots()
sns.heatmap(potential_df, ax = ax)
sns.heatmap(potential_df, ax = ax_total[2])
fig.savefig('../../Figures/TCRen.pdf')

fig_total.savefig('../../Figures/Potentials.pdf')

 