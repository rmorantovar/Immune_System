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
#style.use('seaborn-paper')

M1 = np.loadtxt(Text_files_path+'MJ.txt', skiprows= 1, usecols=range(1,21)).tolist()
M2 = (np.loadtxt(Text_files_path+'MJ2.txt', skiprows= 1, usecols=range(1,21))).tolist()
M3 = np.loadtxt(Text_files_path+'BLOSUM62.txt', skiprows= 1, max_rows = 23, usecols=range(1,24)).tolist()
Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
Alphabet2 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w']
Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)


NC = 2e5
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

lambd = 0.65

for energy_model in energy_models:
    
    fig, ax = plt.subplots(2,4,figsize=(40,18), gridspec_kw={'hspace':0.25, 'wspace':.25})
    fig_b, ax_b = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
    fig_a, ax_a = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
    for j, linear in enumerate(growth_models):

        colors_activation = []
        data_antigen = np.loadtxt(Text_files_path + 'Dynamics/Single_trajectory/antigen_L-%d_N-%d_Antigen-'%(L, NC)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_linear-%d_'%(alpha, beta, gamma, linear)+energy_model+'.txt')
        data_bcells = np.loadtxt(Text_files_path + 'Dynamics/Single_trajectory/bcells_L-%d_N-%d_Antigen-'%(L, NC)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_linear-%d_'%(alpha, beta, gamma, linear)+energy_model+'.txt')
        data_N_active_linages = np.loadtxt(Text_files_path + 'Dynamics/Single_trajectory/N_active_linages_L-%d_N-%d_Antigen-'%(L, NC)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_linear-%d_'%(alpha, beta, gamma, linear)+energy_model+'.txt')
        data_bcells_active = np.transpose(data_bcells[:,np.where(data_bcells[-1,:]>1)[0]])
        print('Activated clones:',data_N_active_linages[-1], np.shape(data_bcells_active))
        
        final_clone_sizes = data_bcells_active[:, -1]
        final_clone_sizes_sorted = np.flip(np.sort(final_clone_sizes))

        t0 = np.where(data_N_active_linages==0)[0][-1]*dT
        a0 = data_antigen[0]

        for i in data_bcells[-1,:]:
            if(i==1):
                colors_activation.append('silver')
            else:
                colors_activation.append(colors[j])

        #---- Antigen ----
        ax[0,0].plot(time, data_antigen, color = colors[j], label = models_name[j] + ' growth', linewidth = 4, linestyle  = '-', marker = '', ms = 10)
        if(linear == 0):
            expfit = 2*np.exp(alpha*time)
            ax_a.plot(time[::500], data_antigen[::500], color = 'cornflowerblue', linewidth = 2, linestyle  = '-', marker = '*', ms = 9);
            ax_a.plot(time[::500], expfit[::500], color = colors_fit[j], linestyle = '--', linewidth = 4)
            ax_a.text(x=8, y=1e8, s = r'$\sim e^{\alpha t}$', fontsize=48, color = colors_fit[j])
            my_plot_layout(ax = ax_a, yscale = 'log', xlabel = 'Time', ylabel = 'Antigen')
            ax_a.set_ylim(bottom = .8)
            fig_a.savefig('../../Figures/1_Dynamics/Antigen_expansion_'+energy_model+'.png')
        if(linear == 1):
            ax[0,0].plot(time, 1e3 + time*2000, color = colors[j], label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
        my_plot_layout(ax = ax[0,0], yscale = 'log', xlabel = 'Time', ylabel = 'Antigen')
        #ax[0,0].set_ylim(bottom = 1)

        #---- B cell linages ----
        for i in np.arange(int(data_N_active_linages[-1])):
            ax[0,1].plot(time, data_bcells_active[i,:], color = colors[j], linewidth = 4, linestyle  = '-', marker = '', ms = 12);
        my_plot_layout(ax = ax[0,1], yscale = 'log', xlabel = 'Time', ylabel = 'Clone size')
        if(linear == 0):
            expfit = 1e-1*np.exp(beta*time)
            for i in range(int(data_N_active_linages[-1])):
                ax_b.plot(time[::500], data_bcells_active[i,:][::500], color = 'royalblue', linewidth = 1, linestyle  = '-', marker = 'o', ms = 5);
            ax_b.plot(time[::500], expfit[::500], color = colors_fit[j], linestyle = '--', linewidth = 4)
            ax_b.text(x=8, y=1e3, s = r'$\sim e^{\beta t}$', fontsize=48, color = colors_fit[j])
            my_plot_layout(ax = ax_b, yscale = 'log', xlabel = 'Time', ylabel = 'Clone size')
            ax_b.set_ylim(bottom = .8)
            fig_b.savefig('../../Figures/1_Dynamics/B_cells_expansion_'+energy_model+'.png')
        #---- Activation rate ----
        ax[1,0].plot(time[::100], data_N_active_linages[::100], linestyle = '-', marker = 'o', ms = 5, linewidth = 2, label = 'simulation', color = colors[j])
        if(linear == 0):
            theory = (1/(alpha*lambd))*(np.exp(alpha*lambd*(time-t0))-1)
            ax[1,0].plot(time, theory, color = colors[j], label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
        if(linear == 1):
            theory = ((1/(alpha*2000*(lambd)*(alpha*2000*t0+a0)**(lambd-1)))*((alpha*2000*time+a0)**(lambd)-(alpha*2000*t0+a0)**(lambd)))
            #theory = time**(lambd)
            ax[1,0].plot(time, theory, color = colors[j], label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
        my_plot_layout(ax = ax[1,0], xscale = 'linear', yscale = 'log', xlabel = 'Time', ylabel = r'$m(t)$')

        #---- Entropy ----
        total_pop = np.sum(np.transpose(data_bcells), axis = 0)
        bcell_freqs = np.transpose(data_bcells)/total_pop
        entropy = -np.sum(bcell_freqs*np.log(bcell_freqs), axis = 0)
        ax[1,1].plot(time, entropy/entropy[0], color = colors[j], linewidth = 3)
        my_plot_layout(ax = ax[1,1], yscale = 'log', xlabel = 'Time', ylabel = 'Entropy')

        #---- Total Pop size ----
        ax[0,2].plot(time, total_pop, color = colors[j], linewidth = 3)
        my_plot_layout(ax = ax[0,2], xlabel = 'Time', ylabel = 'Total Pop. size', yscale = 'log')

        #---- Relative clone sizes ----
        #ax[1,2].plot(range(1,int(8+1)), final_clone_sizes_sorted[:8]/final_clone_sizes_sorted[0], color = colors[j], linewidth = .5, linestyle = '--', marker = '^', ms = 15)
        my_plot_layout(ax = ax[1,2], xlabel = 'Clone', ylabel = 'Relative size', yscale = 'log')
        #ax[1,2].set_xlim(right = 8)

        #---- Stackplots ----
        ax[j, 3].stackplot(time, bcell_freqs, colors = colors_activation);
        my_plot_layout(ax = ax[j, 3], ticks_labelsize=24, title=models_name[j], xlabel = 'time', ylabel = 'Clone frequency')

    fig.savefig('../../Figures/1_Dynamics/summary_1_single_trajectory_'+energy_model+'.png')
    

