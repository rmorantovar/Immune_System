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

peptides = pd.read_csv('../../in/benchmark_candidate_epitopes_IEDB.txt')
peptides = np.array(peptides['peptide'])

data = pd.read_csv('../../in/iedb_slim.csv')
data = data.loc[data['MHC.3']=='I']
peptides2 = np.unique(data['Epitope.2'])

energy_models = ['TCRen', 'MJ2']
for energy_model in energy_models:
	E_ms_data = []
	ls_data = []
	for p in tqdm(peptides[::1]):
	    p_seq = from_aa_to_i_X(p, energy_model, '../../')
	    l = len(p_seq)
	    if l == 9:    
		    ls_data.append(l)
		    motif = get_motif(p_seq, energy_model, '../../')
		    E_m = np.sum(np.min(motif, axis = 0))
		    E_ms_data.append(E_m)
	print(l)

	E_ms_data2 = []
	ls_data2 = []
	for p in tqdm(peptides2[::1]):
	    p_seq = from_aa_to_i_X(p, energy_model, '../../')
	    l = len(p_seq)
	    if l == 9:    
		    ls_data2.append(l)
		    motif = get_motif(p_seq, energy_model, '../../')
		    E_m = np.sum(np.min(motif, axis = 0))
		    E_ms_data2.append(E_m)
	print(l)
	E_ms_random = []
	for i in tqdm(range(100000)):
		p_seq = np.random.randint(0, 20, 9)    
		motif = get_motif(p_seq, energy_model, '../../')
		E_m = np.sum(np.min(motif, axis = 0))
		E_ms_random.append(E_m)

	special_epitopes = ['KQWLVWLFL', 'RLLHPHHPL', 'ALDAPLFGI', 'MRTOPHHPY', 'AFPYNMLRT', 'QTMNAPRST']
	E_m_special = []
	for special_epitope in special_epitopes:
		p_seq = from_aa_to_i_X(special_epitope, energy_model, '../../')
		motif = get_motif(p_seq, energy_model, '../../')
		E_m = np.sum(np.min(motif, axis = 0))
		E_m_special.append(E_m)

	bins = np.linspace(np.min(E_ms_data2)-3, np.max(E_ms_data2)+2, 80)
	#bins = 20
	hist_data = np.histogram(E_ms_data, bins = bins, density = False )
	hist_data2 = np.histogram(E_ms_data2, bins = bins, density = False )
	hist_random = np.histogram(E_ms_random, bins = bins, density = False )
	fig, ax = plt.subplots()
	ax.plot(hist_data[1][:-1], np.cumsum(hist_data[0])/len(E_ms_data), label = 'data', color = my_blue)
	ax.plot(hist_data2[1][:-1], np.cumsum(hist_data2[0])/len(E_ms_data2), label = 'data2', color = my_green)
	ax.plot(hist_random[1][:-1], np.cumsum(hist_random[0])/len(E_ms_random), label = 'random', color = my_red)
	ax.vlines(E_m_special, ax.get_ylim()[0],ax.get_ylim()[1], ls = '--', colors = [my_blue, my_blue, my_blue, my_red, my_red, my_red])
	ax.vlines(np.min(E_ms_data), ax.get_ylim()[0],ax.get_ylim()[1], ls = '-', color = my_purple)
	#ax.set_yscale('log')
	ax.legend()
	fig.savefig('../../../Figures/antigenicity/F_E_m_' + energy_model + '.pdf')

	fig, ax = plt.subplots()
	ax.plot(hist_data[1][:-1], (hist_data[0])/len(E_ms_data), label = 'data', color = my_blue)
	ax.plot(hist_data2[1][:-1], (hist_data2[0])/len(E_ms_data2), label = 'data2', color = my_green)
	ax.plot(hist_random[1][:-1], (hist_random[0])/len(E_ms_random), label = 'random', color = my_red)
	ax.vlines(E_m_special, ax.get_ylim()[0],ax.get_ylim()[1], ls = '--', colors = [my_blue, my_blue, my_blue, my_red, my_red, my_red])
	ax.vlines(np.min(E_ms_data), ax.get_ylim()[0],ax.get_ylim()[1], ls = '-', color = my_purple)
	#ax.set_yscale('log')
	ax.legend()
	fig.savefig('../../../Figures/antigenicity/p_E_m_' + energy_model + '.pdf')

	fig, ax = plt.subplots()
	ax.plot(hist_data[1][:-1], (hist_data[0])/len(E_ms_data), label = 'data', color = my_blue)
	ax.plot(hist_data2[1][:-1], (hist_data2[0])/len(E_ms_data2), label = 'data2', color = my_green)
	ax.plot(hist_random[1][:-1], (hist_random[0])/len(E_ms_random), label = 'random', color = my_red)
	ax.vlines(E_m_special, ax.get_ylim()[0],ax.get_ylim()[1], ls = '--', colors = [my_blue, my_blue, my_blue, my_red, my_red, my_red])
	ax.vlines(np.min(E_ms_data), ax.get_ylim()[0],ax.get_ylim()[1], ls = '-', color = my_purple)
	ax.set_yscale('log')
	ax.legend()
	fig.savefig('../../../Figures/antigenicity/p_E_m_log_' + energy_model + '.pdf')

	fig, ax = plt.subplots()
	ax.plot(hist_data[1][:-1], np.log(((hist_data[0])/len(E_ms_data))/((hist_random[0])/len(E_ms_random))), color = my_cyan, label = 'data')
	ax.plot(hist_data2[1][:-1], np.log(((hist_data2[0])/len(E_ms_data2))/((hist_random[0])/len(E_ms_random))), color = my_gold, label = 'data2')
	ax.vlines(E_m_special, ax.get_ylim()[0],ax.get_ylim()[1], ls = '--', colors = [my_blue, my_blue, my_blue, my_red, my_red, my_red])
	ax.vlines(np.min(E_ms_data), ax.get_ylim()[0],ax.get_ylim()[1], ls = '-', color = my_purple)
	fig.savefig('../../../Figures/antigenicity/KL_d_' + energy_model + '.pdf')

	fig, ax = plt.subplots()
	ax.hist(ls_data)
	ax.hist(ls_data2)
	fig.savefig('../../../Figures/antigenicity/ls_' + energy_model + '.pdf')



