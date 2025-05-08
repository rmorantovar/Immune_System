import sys
sys.path.append('../../lib/')
from functions_1 import*
from classes import*
import numpy as np
import pickle
from tqdm import tqdm
import pandas as pd
import math
import matplotlib.pyplot as plt
import ast

#-----------
peptides = pd.read_csv('../../in/benchmark_candidate_epitopes_IEDB.txt')
peptides = np.array(peptides['peptide'])

data = pd.read_csv('../../in/iedb_slim.csv')
data = data.loc[data['MHC.3']=='I']
peptides2 = np.unique(data['Epitope.2'])

data = pd.read_csv('../../in/epitope_iedb_new_H.csv')
peptides3 = np.unique(data['Epitope'])

Homo_data = open('../../in/Homo_sapiens.faa', 'r')
Homo_data = Homo_data.readlines()
Homo_proteins = ''
for row in Homo_data:
    if list(row)[0]!= '>' :
        Homo_proteins = Homo_proteins+row
Homo_proteins = np.array(list(Homo_proteins))
Homo_proteins = Homo_proteins[Homo_proteins!='\n']
ninemers_Homo = [Homo_proteins[i: i+9] for i in range(len(Homo_proteins)-9)]

Influenza_A_data = open('../../in/Influenza_A.faa', 'r')
Influenza_A_data = Influenza_A_data.readlines()
Influenza_A_proteins = ''
for row in Influenza_A_data:
    if list(row)[0]!= '>' :
        Influenza_A_proteins = Influenza_A_proteins + row
Influenza_A_proteins = np.array(list(Influenza_A_proteins))
Influenza_A_proteins = Influenza_A_proteins[Influenza_A_proteins!='\n']
ninemers_Influenza_A = [Influenza_A_proteins[i: i+9] for i in range(len(Influenza_A_proteins)-9)]
#-----------

#-----------
potential = pd.read_csv('../../in/TCRen_potential.csv')
Alphabet_T = np.unique(potential['residue.aa.from'])
Alphabet_p = np.unique(potential['residue.aa.to'])
Alphabet_list = Alphabet_p.tolist()

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
#-----------

energy_models = ['TCRen']
for energy_model in energy_models:

	E_ms_dataH = []
	ls_dataH = []
	for p in tqdm(ninemers_Homo[::10]):
		l = len(p)
		if not(np.sum(np.isin(p, ['U', 'X']))):
			aa_list = [i for i in p]
			p_seq = []
			for i, aa in enumerate(aa_list):
				index = Alphabet_list.index(aa)
				p_seq.append(int(index))
			ls_dataH.append(l)
			motif = potential[:,p_seq]
			E_m = np.sum(np.min(motif, axis = 0))
			E_ms_dataH.append(E_m)

	E_ms_dataFlu = []
	ls_dataFlu = []
	for p in tqdm(ninemers_Influenza_A[::1]):
		l = len(p)
		if not(np.sum(np.isin(p, ['U', 'X']))):
			aa_list = [i for i in p]
			p_seq = []
			for i, aa in enumerate(aa_list):
				index = Alphabet_list.index(aa)
				p_seq.append(int(index))
			ls_dataFlu.append(l)
			motif = potential[:,p_seq]
			E_m = np.sum(np.min(motif, axis = 0))
			E_ms_dataFlu.append(E_m)

	E_ms_data = []
	ls_data = []
	for p in tqdm(peptides[::1]):
	    l = len(p)
	    if(l==9):
		    aa_list = [i for i in p]
		    p_seq = []
		    for i, aa in enumerate(aa_list):
		        index = Alphabet_list.index(aa)
		        p_seq.append(int(index))	       
		    ls_data.append(l)
		    motif = potential[:,p_seq]
		    E_m = np.sum(np.min(motif, axis = 0))
		    E_ms_data.append(E_m)
	
	E_ms_data2 = []
	ls_data2 = []
	for p in tqdm(peptides2[::1]):
	    l = len(p)
	    if(l==9):
		    aa_list = [i for i in p]
		    p_seq = []
		    for i, aa in enumerate(aa_list):
		        if aa in 'XZ':
		            break
		        index = Alphabet_list.index(aa)
		        p_seq.append(int(index))	       
		    ls_data2.append(l)
		    motif = potential[:,p_seq]
		    E_m = np.sum(np.min(motif, axis = 0))
		    E_ms_data2.append(E_m)

	E_ms_data3 = []
	ls_data3 = []
	for p in tqdm(peptides3[::1]):
	    l = len(p)
	    if(l==9):
		    aa_list = [i for i in p]
		    p_seq = []
		    for i, aa in enumerate(aa_list):
		        if aa in 'XZl':
		            break
		        index = Alphabet_list.index(aa)
		        p_seq.append(int(index))	       
		    ls_data3.append(l)
		    motif = potential[:,p_seq]
		    E_m = np.sum(np.min(motif, axis = 0))
		    E_ms_data3.append(E_m)
			

	E_ms_random = []
	for i in tqdm(range(1000000)):
		p_seq = np.random.randint(0, 20, 9)    
		motif = potential[:,p_seq]
		E_m = np.sum(np.min(motif, axis = 0))
		E_ms_random.append(E_m)
	special_epitopes = ['KQWLVWLFL', 'RLLHPHHPL', 'ALDAPLFGI', 'MRTOPHHPY', 'AFPYNMLRT', 'QTMNAPRST']
	special_epitopes = ['KQWLVWLFL', 'RLLHPHHPL', 'ALDAPLFGI']
	E_m_special = []
	for special_epitope in special_epitopes:
		p_seq = from_aa_to_i_X(special_epitope, energy_model, '../../')
		motif = potential[:,p_seq]
		E_m = np.sum(np.min(motif, axis = 0))
		E_m_special.append(E_m)

	print(len(E_ms_data), len(E_ms_data2), len(E_ms_data3), len(E_ms_dataH))

	bins = np.linspace(np.min(E_ms_data2)-1, np.max(E_ms_data2)+1, 80)
	#bins = 20
	# hist_data = np.histogram(E_ms_data, bins = bins, density = False )
	# hist_data2 = np.histogram(E_ms_data2, bins = bins, density = False )
	hist_data3 = np.histogram(E_ms_data3, bins = bins, density = False )
	hist_dataH = np.histogram(E_ms_dataH, bins = bins, density = False )
	hist_dataFlu = np.histogram(E_ms_dataFlu, bins = bins, density = False )
	hist_random = np.histogram(E_ms_random, bins = bins, density = False )

	fig, ax = plt.subplots(figsize = (12, 10))
	# ax.plot(hist_data[1][:-1], np.cumsum(hist_data[0])/len(E_ms_data), label = 'data', color = my_blue, ls = '', marker = 's', alpha = .6, ms = 8)
	# ax.plot(hist_data2[1][:-1], np.cumsum(hist_data2[0])/len(E_ms_data2), label = 'data2', color = my_blue, ls = '', marker = 'o', alpha = .6, ms = 8)
	ax.plot(hist_data3[1][:-1], np.cumsum(hist_data3[0])/len(E_ms_data3), label = 'data3', color = my_blue, ls = '', marker = '^', alpha = .6, ms = 8)
	ax.plot(hist_dataH[1][:-1], np.cumsum(hist_dataH[0])/len(E_ms_dataH), label = 'dataH', color = my_blue, ls = '', marker = '*', alpha = .6, ms = 8)
	ax.plot(hist_dataFlu[1][:-1], np.cumsum(hist_dataFlu[0])/len(E_ms_dataFlu), label = 'dataFlu', color = my_blue, ls = '', marker = 'o', alpha = .6, ms = 8)
	ax.plot(hist_random[1][:-1], np.cumsum(hist_random[0])/len(E_ms_random), label = 'random', color = my_red, alpha = .6)
	# ax.vlines(E_m_special, ax.get_ylim()[0],ax.get_ylim()[1], ls = '--', colors = [my_blue, my_blue2, my_blue2])
	# ax.vlines(np.min(E_ms_data), ax.get_ylim()[0],ax.get_ylim()[1], ls = '-', color = my_purple)
	#ax.set_yscale('log')
	ax.legend()
	fig.savefig('../../../Figures/antigenicity/F_E_m_' + energy_model + '0.pdf')

	fig, ax = plt.subplots(figsize = (12, 10))
	# ax.plot(hist_data[1][:-1], (hist_data[0])/len(E_ms_data), label = 'data', color = my_blue, ls = '', marker = 's', alpha = .6, ms = 8)
	# ax.plot(hist_data2[1][:-1], (hist_data2[0])/len(E_ms_data2), label = 'data2', color = my_blue, ls = '', marker = 'o', alpha = .6, ms = 8)
	ax.plot(hist_data3[1][:-1], (hist_data3[0])/len(E_ms_data3), label = 'data3', color = my_blue, ls = '', marker = '^', alpha = .6, ms = 8)
	ax.plot(hist_dataH[1][:-1], (hist_dataH[0])/len(E_ms_dataH), label = 'dataH', color = my_blue, ls = '', marker = '*', alpha = .6, ms = 8)
	ax.plot(hist_dataFlu[1][:-1], (hist_dataFlu[0])/len(E_ms_dataFlu), label = 'dataFlu', color = my_blue, ls = '', marker = 'o', alpha = .6, ms = 8)
	ax.plot(hist_random[1][:-1], (hist_random[0])/len(E_ms_random), label = 'random', color = my_red, alpha = .6)
	# ax.vlines(E_m_special, ax.get_ylim()[0],ax.get_ylim()[1], ls = '--', colors = [my_blue, my_blue2, my_blue2])
	# ax.vlines(np.min(E_ms_data), ax.get_ylim()[0],ax.get_ylim()[1], ls = '-', color = my_purple)
	#ax.set_yscale('log')
	ax.legend()
	fig.savefig('../../../Figures/antigenicity/p_E_m_' + energy_model + '0.pdf')

	fig, ax = plt.subplots(figsize = (12, 10))
	# ax.plot(hist_data[1][:-1], (hist_data[0])/len(E_ms_data), label = 'data', color = my_blue, ls = '', marker = 's', alpha = .6, ms = 8)
	# ax.plot(hist_data2[1][:-1], (hist_data2[0])/len(E_ms_data2), label = 'data2', color = my_blue, ls = '', marker = 'o', alpha = .6, ms = 8)
	ax.plot(hist_data3[1][:-1], (hist_data3[0])/len(E_ms_data3), label = 'data3', color = my_blue, ls = '', marker = '^', alpha = .6, ms = 8)
	ax.plot(hist_dataH[1][:-1], (hist_dataH[0])/len(E_ms_dataH), label = 'dataH', color = my_blue, ls = '', marker = '*', alpha = .6, ms = 8)
	ax.plot(hist_dataFlu[1][:-1], (hist_dataFlu[0])/len(E_ms_dataFlu), label = 'dataFlu', color = my_blue, ls = '', marker = 'o', alpha = .6, ms = 8)
	ax.plot(hist_random[1][:-1], (hist_random[0])/len(E_ms_random), label = 'random', color = my_red, alpha = .6)
	# ax.vlines(E_m_special, ax.get_ylim()[0],ax.get_ylim()[1], ls = '--', colors = [my_blue, my_blue2, my_blue2])
	# ax.vlines(np.min(E_ms_data), ax.get_ylim()[0],ax.get_ylim()[1], ls = '-', color = my_purple)
	ax.set_yscale('log')
	ax.legend()
	fig.savefig('../../../Figures/antigenicity/p_E_m_log_' + energy_model + '0.pdf')

	fig, ax = plt.subplots(figsize = (12, 10))
	# ax.plot(hist_data[1][:-1], np.log(((hist_data[0])/len(E_ms_data))/((hist_random[0])/len(E_ms_random))), color = my_green, label = 'data')
	# ax.plot(hist_data2[1][:-1], np.log(((hist_data2[0])/len(E_ms_data2))/((hist_random[0])/len(E_ms_random))), color = my_gold, label = 'data2')
	ax.plot(hist_data3[1][:-1], np.log(((hist_data3[0])/len(E_ms_data3))/((hist_random[0])/len(E_ms_random))), color = my_purple, label = 'data3/random')
	ax.plot(hist_dataH[1][:-1], np.log(((hist_dataH[0])/len(E_ms_dataH))/((hist_random[0])/len(E_ms_random))), color = my_red, label = 'dataH/random')
	ax.plot(hist_dataFlu[1][:-1], np.log(((hist_dataFlu[0])/len(E_ms_dataFlu))/((hist_random[0])/len(E_ms_random))), color = my_green, label = 'dataFlu/random')
	ax.plot(hist_dataH[1][:-1], np.log(((hist_data3[0])/len(E_ms_data3))/((hist_dataH[0])/len(E_ms_dataH))), color = my_purple, label = 'data3/dataH', ls = '--')
	ax.plot(hist_dataH[1][:-1], np.log(((hist_dataH[0])/len(E_ms_dataH))/((hist_dataH[0])/len(E_ms_dataH))), color = my_red, label = 'dataH/dataH', ls = '--')
	ax.plot(hist_dataH[1][:-1], np.log(((hist_dataFlu[0])/len(E_ms_dataFlu))/((hist_dataH[0])/len(E_ms_dataH))), color = my_green, label = 'dataFlu/dataH', ls = '--')
	# ax.vlines(E_m_special, ax.get_ylim()[0],ax.get_ylim()[1], ls = '--', colors = [my_blue, my_blue2, my_blue2])
	# ax.vlines(np.min(E_ms_data), ax.get_ylim()[0],ax.get_ylim()[1], ls = '-', color = my_purple)
	ax.legend()
	fig.savefig('../../../Figures/antigenicity/KL_d_' + energy_model + '0.pdf')

	fig, ax = plt.subplots(figsize = (12, 10))
	ax.hist(ls_data)
	ax.hist(ls_data2)
	fig.savefig('../../../Figures/antigenicity/ls_' + energy_model + '0.pdf')



