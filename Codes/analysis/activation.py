import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
from Immuno_models import*
import scipy.special as sc
import pickle
from tqdm import tqdm

Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/"


fig_p_a, ax_p_a = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15, 'bottom': .15, 'right':.8})
# Paramters:
N_engagements = 1e2
k_on = 1e6
m_off = 0
m_on_array = np.array([1e-4, 1])
m_on_array = np.array([1e-2, 1, 1e2])
markers = ['^', 's', 'o']
colors = ['darkred', 'darkblue', 'green']
labels = [r'$m_{\mathrm{on}}\ll \gamma$', r'$m_{\mathrm{on}}\approx \gamma$', r'$m_{\mathrm{on}}\gg \gamma$']
#l_off = k_off
l_on = 0
gamma = 1
k_off_array2 = np.logspace(-6, 2, 100)
k_off_array2 = np.logspace(-5, 4, 100)

for j, m_on in enumerate(m_on_array):
	file_p_a = open(Text_files_path+'recognition_gamma-%.0e_m_on-%.0e.pkl'%(gamma, m_on),'rb')
	k_off_array, p_a_array = pickle.load(file_p_a)

	ax_p_a.plot(k_off_array/k_on, p_a_array, marker = markers[j], linestyle = '', ms = 8, color= colors[j], label = labels[j], alpha = .8)
	ax_p_a.plot(k_off_array2/k_on, (m_on_array[j]*gamma)/((m_on_array[j]*gamma) + (m_on_array[j] + gamma)*k_off_array2 +k_off_array2**2), marker = '', linestyle = '-', ms = 8, color = colors[j])
	#ax_p_a.plot(k_off_array2/k_on, np.ones_like(k_off_array2), marker = '', linestyle = '--', ms = 8, color = colors[j], alpha = .6)
	#ax_p_a.plot(k_off_array2/k_on, (m_on_array[j]*gamma)/((m_on_array[j]*gamma)+k_off_array2**2), marker = '', linestyle = '--', ms = 8, color = colors[j], alpha = .6)

	#if(m_on<gamma):
	#	ax_p_a.plot(k_off_array2/k_on, m_on_array[j]/(m_on_array[j]+k_off_array2), marker = '', linestyle = '--', ms = 8, color = colors[j], alpha = .6)
	#if(m_on>gamma):
	#	ax_p_a.plot(k_off_array2/k_on, gamma/(gamma+k_off_array2), marker = '', linestyle = '--', ms = 8, color = colors[j], alpha = .6)

ax_p_a.vlines([np.concatenate((m_on_array/k_on, [0]))], 1e-7, 1, linestyle=':', alpha = .4, color='grey' )
my_plot_layout(ax=ax_p_a, yscale = 'log', xscale = 'log', ticks_labelsize = 24, xlabel = r'$k_{\mathrm{off}}$', ylabel = r'$p_a$', title = '', x_fontsize=24, y_fontsize = 24, t_fontsize = 24)
#ax_p_a.set_ylim(bottom = 1e-7)
ax_p_a.legend(loc = 3, fontsize = 24, title = r'$\gamma=%.1f$'%(gamma), title_fontsize = 24)
fig_p_a.savefig('../../Figures/7_Recognition/K-PR_model.pdf')


