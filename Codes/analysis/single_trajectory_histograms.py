import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
from Immuno_models import*
import scipy.special as sc
import pickle
from matplotlib import style
from scipy.optimize import curve_fit

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'


N_A = 6.02214076e23
k_BT = 1.380649e-23*293
style.use('seaborn-paper')

M1 = np.loadtxt(Text_files_path+'MJ.txt', skiprows= 1, usecols=range(1,21)).tolist()
M2 = (np.loadtxt(Text_files_path+'MJ2.txt', skiprows= 1, usecols=range(1,21))).tolist()
M3 = np.loadtxt(Text_files_path+'BLOSUM62.txt', skiprows= 1, max_rows = 23, usecols=range(1,24)).tolist()
Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
Alphabet2 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w']
Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)


NC = 2e4
T0 = 0
Tf = 25
dT = 0.005
alpha = 1
gamma = 0.0000
beta = 0.5

antigen = 'TACNSEYPNTTK'
antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FTSENAYCGR'
antigen = 'FMLFMAVFVMTSWYC'

L=len(antigen)

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
models_name = ['exponential', 'linear', ]
colors = ['tab:blue', 'tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = [0, 1]

#fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
#fig3, ax3 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})

colors_fate = [['darkblue', 'royalblue'], ['darkred', 'indianred']]

min_e = -76
max_e = -58
gauge_e = 50

lambd = 0.5976041516921407
lambd - 0.7313511981449714

for energy_model in energy_models:
    
    fig, ax = plt.subplots(2,4,figsize=(40,18), gridspec_kw={'hspace':0.25})

    for j, linear in enumerate(growth_models):

        fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
        fig3, ax3 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})

        colors_activation = []
        data_bcells = np.loadtxt(Text_files_path + 'Dynamics/Single_trajectory/bcells_L-%d_N-%d_Antigen-'%(L, NC)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_linear-%d_'%(alpha, beta, gamma, linear)+energy_model+'.txt')
        data_energies = np.loadtxt(Text_files_path + 'Dynamics/Single_trajectory/energies_L-%d_N-%d_Antigen-'%(L, NC)+antigen+'_linear-%d_'%(linear)+energy_model+'.txt')

        #min_e = np.min(data_energies[:,0])
        #max_e = np.max(data_energies[:,0])
        print(min_e, max_e)

        energies_array = np.linspace(min_e, max_e-4, 50)

        data_bcells_active = np.transpose(data_bcells[:,np.where(data_energies[:,1]==1)[0]])

        energies_active = data_energies[np.where(data_energies[:,1]==1)[0],:]

        #print(np.min(data_energies))

        #---- DISTRIBUTION ENERGIES ----
        ax[j,0].hist(data_energies[:,0], bins = np.linspace(min_e, max_e, 20), color = colors[j], alpha = .8, histtype = 'bar')
        my_plot_layout(ax = ax[j,0], yscale = 'log', xlabel = 'Energies', ylabel = 'counts')
        ax[j,0].set_ylim(bottom = .8)
        if(linear == 0):
            fig_seq, ax_seq = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
            expfit = np.exp(lambd*(energies_array+74.5))
            ax_seq.plot(energies_array, expfit, color = 'indigo', linestyle = '--', linewidth = 4, alpha = .7)
            hist = np.histogram(data_energies[:,0], bins = np.linspace(min_e, max_e, 20))
            ax_seq.plot(hist[1][:-1][:-3], hist[0][:-3], color = 'indigo', alpha = 1, linestyle = '', marker = '^', ms = 10)
            ax_seq.text(x=-74, y=8e1, s = r'$\sim e^{\lambda \epsilon}$', fontsize=48, color = 'indigo')
            my_plot_layout(ax = ax_seq, yscale = 'log', xlabel = r'$\epsilon$', ylabel = r'$\Lambda(\epsilon)$')
            ax_seq.set_ylim(bottom = .8)
            fig_seq.savefig('../../Figures/1_Dynamics/Sequences_expansion_'+energy_model+'.pdf')

        #---- DISTRIBUTION ACTIVATED ENERGIES ----
        ax[j,1].hist(energies_active[:,0], bins = np.linspace(min_e, max_e, 20), color = colors[j], alpha = .8, histtype = 'bar')
        my_plot_layout(ax = ax[j,1], yscale = 'log', xlabel = 'Energies', ylabel = 'counts')
        ax[j,1].set_ylim(bottom = ax[j,0].get_ylim()[0], top=ax[j,0].get_ylim()[1])

        #---- DISTRIBUTION GC VS PLASMA ----
        ax[j,2].hist([energies_active[np.where(energies_active[:,2]==1)[0],0], energies_active[np.where(energies_active[:,2]==0)[0],0]], bins = np.linspace(min_e, max_e, 20), color = colors_fate[j], alpha = .6, histtype = 'barstacked', label = ['GC', 'Plasma'])
        my_plot_layout(ax = ax[j,2], yscale = 'log', xlabel = 'Energies', ylabel = 'counts')
        ax[j,2].set_ylim(bottom = ax[j,0].get_ylim()[0], top=ax[j,0].get_ylim()[1])
        ax[j,2].legend(fontsize = 24)

        ax2.hist([energies_active[np.where(energies_active[:,2]==1)[0],0], energies_active[np.where(energies_active[:,2]==0)[0],0]], bins = np.linspace(min_e, max_e, 20), color = colors_fate[j], alpha = .6, histtype = 'barstacked', label = ['GC', 'Plasma'])
        my_plot_layout(ax = ax2, yscale = 'log', xlabel = 'Energies', ylabel = 'counts')
        ax2.set_ylim(bottom = ax[j,0].get_ylim()[0], top=ax[j,0].get_ylim()[1])
        ax2.legend(fontsize = 24)

        #---- SERA ----
        clone_sizes = data_bcells_active[np.where(energies_active[:,2]==0)[0],-1]
        Kds = np.exp(energies_active[np.where(energies_active[:,2]==0)[0],0]+gauge_e)
        data_Kds = ax[j,3].hist(Kds, bins = np.logspace(np.log10(np.exp(min_e+gauge_e)), np.log10(np.exp(max_e+gauge_e)), 20), density = False, color = colors_fate[j][1], alpha = .5, histtype = 'bar')
        Kds_array = (data_Kds[1][:-1]+data_Kds[1][1:])/2
        ax[j,3].hist(Kds, bins = np.logspace(np.log10(np.exp(min_e+gauge_e)), np.log10(np.exp(max_e+gauge_e)), 20), density = False, color = 'silver', alpha = .4, histtype = 'bar', weights = data_bcells_active[np.where(energies_active[:,2]==0)[0],-1])
        clone_sizes_binned = np.array([])
        for i in range(int(len(data_Kds[0]))):
            clone_sizes_binned = np.append(clone_sizes_binned, np.mean(clone_sizes[(Kds<data_Kds[1][i+1]) & (Kds>data_Kds[1][i])]))
        ax[j,3].plot(Kds_array, clone_sizes_binned, color = 'dimgray', linewidth =3, linestyle = '', marker = '*', ms = 15, label = 'Clone size')
        if(linear==0):
        	ax[j,3].plot(Kds_array[:], 2e6*Kds_array[:]**(lambd), color = colors_fate[j][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
        	ax[j,3].plot(Kds_array[:], np.exp(beta*Tf)*np.exp(beta*25/alpha)*N_A**(-beta/alpha)*Kds_array[:]**(-beta/alpha), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
        	ax[j,3].plot(Kds_array[:], 2e6*np.exp(beta*Tf)*np.exp(beta*25/alpha)*N_A**(-beta/alpha)*Kds_array[:]**(lambd-beta/alpha), color = 'silver', linewidth =3, linestyle = '--', marker = '', ms = 15)
        if(linear==1):
        	ax[j,3].plot(Kds_array[:-11], 2e6*Kds_array[:-11]**(lambd), color = colors_fate[j][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
        	ax[j,3].plot(Kds_array[:-11], np.exp(beta*Tf)*np.exp(beta*1e3/(alpha*2000))*np.exp(-beta/(alpha*2000)*np.exp(-25)*N_A*Kds_array[:-11]), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
        ax[j,3].vlines(np.exp(np.average(np.log(Kds))), 0, ax[j,3].get_ylim()[1], color = colors_fate[j][1], linewidth = 3, linestyle = '-')
        #ax[j,3].vlines(np.exp(np.average(np.log(Kds), weights = data_bcells_active[np.where(energies_active[:,2]==0)[0],-1])), 0, ax[j,3].get_ylim()[1], color = colors_fate[j][1], linewidth = 3, linestyle = '--')
        my_plot_layout(ax = ax[j,3], yscale = 'log', xscale = 'log', xlabel = r'$K_D$', ylabel = 'counts')
        ax[j,3].legend(fontsize = 24)

        ax3.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e+gauge_e)), np.log10(np.exp(max_e+gauge_e)), 20), density = False, color = colors_fate[j][1], alpha = .5, histtype = 'bar')
        ax3.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e+gauge_e)), np.log10(np.exp(max_e+gauge_e)), 20), density = False, color = 'silver', alpha = .4, histtype = 'bar', weights = data_bcells_active[np.where(energies_active[:,2]==0)[0],-1])
        ax3.plot(Kds_array, clone_sizes_binned, color = 'dimgray', linewidth =3, linestyle = '', marker = '*', ms = 15, label = 'Clone size')
        if(linear==0):
        	ax3.plot(Kds_array[:], 2e6*Kds_array[:]**(lambd), color = colors_fate[j][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
        	ax3.plot(Kds_array[:], np.exp(beta*Tf)*np.exp(beta*25/alpha)*N_A**(-beta/alpha)*Kds_array[:]**(-beta/alpha), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
        	ax3.plot(Kds_array[:], 2e6*np.exp(beta*Tf)*np.exp(beta*25/alpha)*N_A**(-beta/alpha)*Kds_array[:]**(lambd-beta/alpha), color = 'silver', linewidth =3, linestyle = '--', marker = '', ms = 15)
        if(linear==1):
        	ax3.plot(Kds_array[:-11], 2e6*Kds_array[:-11]**(lambd), color = colors_fate[j][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
        	ax3.plot(Kds_array[:-11], np.exp(beta*Tf)*np.exp(beta*1e3/(alpha*2000))*np.exp(-beta/(alpha*2000)*np.exp(-25)*N_A*Kds_array[:-11]), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
        #ax3.vlines(np.exp(np.average(np.log(Kds))), 0, ax[j,3].get_ylim()[1], color = colors_fate[j][1], linewidth = 3, linestyle = '-')
        #ax3.vlines(np.exp(np.average(np.log(Kds), weights = data_bcells_active[np.where(energies_active[:,2]==0)[0],-1])), 0, ax[j,3].get_ylim()[1], color = colors_fate[j][1], linewidth = 3, linestyle = '--')
        my_plot_layout(ax = ax3, yscale = 'log', xscale = 'log', xlabel = r'$K_D$', ylabel = 'counts')
        ax3.legend(fontsize = 24, loc = 4)

        fig2.savefig("../../Figures/1_Dynamics/Plasma_vs_GC_linear-%d.pdf"%(linear))
        fig3.savefig("../../Figures/1_Dynamics/Sera_linear-%d.pdf"%(linear))


    fig.savefig('../../Figures/1_Dynamics/summary_2_single_trajectory_'+energy_model+'.png')
    
