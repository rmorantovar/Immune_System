import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['text.usetex'] = True
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

Matrix = 'MJ2'
#Matrix = 'MM'

if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()


N_ens = 1
N_r = 5e4
N_r = 1e5
T0 = 0
Tf = 6
dT = 0.01
lambda_A = 6
k_pr = 0.1 # hour^-1
k_pr = k_pr*24 #days^-1

#k_pr= 0.000277
qs = [1, 2]
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e3
E_ms = -27

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'FTSENAYCGR'
antigen = 'TACNSEYPNTTK'
#antigen = 'TANSEYPNTK'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'


L=len(antigen)

#----------------------------------------------------------------
antigen_list = [i for i in antigen]
antigen_seq = np.array([], dtype = int)
for i, aa in enumerate(antigen_list):
    index = Alphabet_list.index(aa)
    antigen_seq = np.append(antigen_seq, int(index))
PWM_data = M2[:,antigen_seq]

#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

Es, dE, Q0, lambdas = calculate_Q0(0.01, 50, PWM_data, E_ms, L)
#----------------------------------------------------------------

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
models_name = ['exponential']#, 'linear',]
colors = ['tab:blue', 'tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = [0]

lambd = 1.4
for q in qs:
    for energy_model in energy_models:
        
        fig, ax = plt.subplots(2,4,figsize=(40,18), gridspec_kw={'hspace':0.25, 'wspace':.25})
        fig_b, ax_b = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
        fig_a, ax_a = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
        fig_act, ax_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
        fig_clones, ax_clones = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.06, 'right':.98, 'bottom':.1, 'top': 0.9})
        for j, linear in enumerate(growth_models):

            colors_activation = []

            parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_q-%d_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, q, j, N_ens)+energy_model
            data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies.txt', sep = '\t', header=None)
            #data_antigen = np.loadtxt(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/antigen.txt')
            #data_bcells = np.loadtxt(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/bcells.txt')

            #data_N_active_linages = np.loadtxt(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/m_bar.txt')
            #data_bcells_active = np.transpose(data_bcells[:,np.where(data_bcells[-1,:]>1)[0]])

            data_active = data.loc[data[1]==1]

            activations_times = np.array(data_active[3])

            ar1, ar2 = np.histogram(activations_times, bins = time)

            data_N_active_linages = np.cumsum(ar1)

            print('Activated clones:',data_N_active_linages[-1], np.shape(data_active))
            
            #final_clone_sizes = data_bcells_active[:, -1]
            #final_clone_sizes_sorted = np.flip(np.sort(final_clone_sizes))

            t0 = np.where(data_N_active_linages==0)[0][-1]*dT

            for i in data[1]:
                if(i==0):
                    colors_activation.append('silver')
                else:
                    colors_activation.append(colors[j])

            #---- Antigen ----

            if(linear == 0):
                expfit = np.exp(lambda_A*time)
                ax[0,0].plot(time, expfit, color = colors[j], label = models_name[j] + ' growth', linewidth = 4, linestyle  = '-', marker = '', ms = 10)
                #ax_a.plot(time, expfit, color = colors_fit[j], linestyle = '--', linewidth = 4)
                #ax_a.text(x=7, y=5e7, s = r'$\sim e^{\lambda_A t}$', fontsize=48, color = colors_fit[j])
            if(linear == 1):
                ax[0,0].plot(time, 1 + time*alpha_lin, color = colors[j], label = models_name[j] + ' growth', linewidth = 4, linestyle  = '-', marker = '', ms = 10)
                #ax_a.plot(time, 1 + time*alpha_lin, color = colors[j], label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
                #ax_a.text(x=16, y=1e2, s = r'$\sim \lambda_A t$', fontsize=48, color = colors_fit[j])
            my_plot_layout(ax = ax[0,0], yscale = 'log', xlabel = 'Time', ylabel = 'Antigen')
            ax[0,0].set_ylim(bottom = 1)
            ax[0,0].set_xlim(3, Tf)

            #---- B cell linages ----
            expfit = 1e-5*np.exp(lambda_B*time)
            clone_sizes = np.ones((int(data_N_active_linages[-1]), len(time)))
            for i in np.arange(0, int(data_N_active_linages[-1])):
                clone_sizes[i, int(activations_times[i]/dT):] = np.exp(lambda_B*(time[int(activations_times[i]/dT):] - activations_times[i] ))
                ax[0,1].plot(time, expfit, color = 'indigo', linestyle = '--', linewidth = 4)
                if(i%5==0):
                    ax[0,1].plot(time, clone_sizes[i,:], color = colors[j], linewidth = 2, linestyle  = '-', marker = '', ms = 12, alpha = .6);
                    ax_b.plot(time, clone_sizes[i,:], color = colors[j], linewidth = 1, linestyle  = '-', marker = '', ms = 5, alpha = .6);
            if(linear == 0):
                #ax_b.plot(time, expfit, color = 'indigo', linestyle = '--', linewidth = 4)
                ax_b.text(x=5, y=3e3, s = r'$\sim e^{\lambda_B t}$', fontsize=48, color = 'indigo')

            my_plot_layout(ax = ax[0,1], yscale = 'log', xlabel = 'Time [days]', ylabel = 'Clone size')
            ax[0,1].set_ylim(bottom = .8)
            ax[0,1].set_xlim(left = 3.5, right = Tf)

            #---- Activation rate ----
            #----------------------------------------------------------------
            u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf)/N_A, Es, q, lambda_A, N_c, dE)
            psi = ((k_on*N_c)/(N_A))*p_a
            M_r = N_r*N_c*np.sum(Q0*p_a*dE)
            ax[1, 2].plot(Es[:-1], p_a)
            ax[1, 2].plot(Es[:-1], N_r*Q0)
            ax[1, 2].plot(Es[:-1], N_r*Q0*p_a)
            my_plot_layout(ax = ax[1,2], xlabel = 'Energy', yscale = 'log')
            print(M_r)
            #----------------------------------------------------------------
            ax[1,0].plot(time[:-1], data_N_active_linages[:], linestyle = '', marker = 'o', ms = 5, linewidth = 2, label = 'simulation', color = colors[j], alpha = .6)
            ax_a.plot(time[:-1], data_N_active_linages[:], linestyle = '', marker = 'o', ms = 5, linewidth = 2, label = 'simulation', color = colors[j], alpha = .6)
            if(linear == 0):
                m_bar = np.array([N_r*(1-np.sum(np.exp(-(psi/lambda_A)*(np.exp(lambda_A*t)-1))*Q0*dE)) for t in time])
                m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))
                ax[1,0].plot(time, m_bar, color = colors[j], label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
                ax[1,0].plot(time, m_bar_approx, color = 'grey', label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
                ax[1,0].plot(time, m_bar_approx, color = 'grey', label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
                ax_a.plot(time, np.ones_like(time), color = 'black', linestyle = '-', linewidth = 1)
                ax_a.plot(time, m_bar, color = colors[j], label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
                ax_a.plot(time, m_bar_approx, color = 'grey', label = models_name[j] + ' growth', linestyle = ':', linewidth = 3)
                #ax[1,0].bar(time[:-1], ar1, alpha = .5)
            if(linear == 1):
                theory = ((k_on*M_r)/(N_A))*(1e8*time + 0.5*lambda_A*time**2)
                ax[1,0].plot(time[time>=t0], m_bar_approx[time>=t0], color = colors[j], label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
            my_plot_layout(ax = ax[1,0], xscale = 'linear', yscale = 'log', xlabel = 'Time', ylabel = r'$m(t)$')
            ax[1,0].set_ylim(top=.5*N_r, bottom = 1e-6)
            ax[1,0].set_xlim(3, Tf)
            #ax_a.set_ylim(top=.5*N_r)


            #---- Entropy ----
            size_lim = 0
            if q==1:
                size_lim=40
            if q==2:
                size_lim=2

            clone_sizes = clone_sizes[np.where(clone_sizes[:,-1]>(np.min(clone_sizes[:,-1])*size_lim))[0],:]
            total_pop = np.sum(clone_sizes, axis = 0)
            bcell_freqs = clone_sizes/total_pop
            entropy = -np.sum(bcell_freqs*np.log(bcell_freqs), axis = 0)
            ax[1,1].plot(time, entropy/entropy[0], color = colors[j], linewidth = 3)
            my_plot_layout(ax = ax[1,1], yscale = 'log', xlabel = 'Time', ylabel = 'Entropy')
            ax[1,1].set_xlim(left = 3.5, right = Tf)

            #---- Total Pop size ----
            ax[0,2].plot(time, total_pop, color = colors[j], linewidth = 3)
            my_plot_layout(ax = ax[0,2], xlabel = 'Time', ylabel = 'Total Pop. size', yscale = 'log')
            ax[0,2].set_xlim(left = 3.5, right = Tf)

            #---- Relative clone sizes ----
            #ax[1,2].plot(range(1,int(8+1)), final_clone_sizes_sorted[:8]/final_clone_sizes_sorted[0], color = colors[j], linewidth = .5, linestyle = '--', marker = '^', ms = 15)
            #my_plot_layout(ax = ax[1,2], xlabel = 'Clone', ylabel = 'Relative size', yscale = 'log')
            #ax[1,2].set_xlim(right = 8)

            #---- Stackplots ----
            #ax[j, 3].stackplot(time, bcell_freqs);
            #ax_clones.stackplot(time[400:], bcell_freqs[np.where(bcell_freqs[:, -1]>.01)[0], 400:]);
            colors = []
            min_bell_freq = np.min(bcell_freqs[:,-1])
            for c in range(int(len(clone_sizes[:,0]))):
                if bcell_freqs[c, -1]>(50*min_bell_freq):
                    if q==1:
                        colors.append('indianred')
                    if q==2:
                        colors.append('dodgerblue')
                else:
                    colors.append('silver')

            ax_clones.stackplot(time, bcell_freqs, colors = colors);

            my_plot_layout(ax = ax[j, 3], ticks_labelsize=30, title=models_name[j], xlabel = 'time', ylabel = 'Clone frequency')
            my_plot_layout(ax = ax_clones, ticks_labelsize=34)
            ax_clones.set_xlim(left = 3.5, right = Tf)
            ax_clones.set_ylim(0, 1)
            ax_clones.set_yticks([])
            ax_clones.set_xticks([3, 4, 5, 6])
        
        my_plot_layout(ax = ax_a, yscale = 'log', xlabel = 'Time', ylabel = r'$\bar m$')
        ax_a.set_xlim(left = 3.5, right = Tf)
        ax_a.set_ylim(top=.5*N_r, bottom = 1e-6)
        fig_a.savefig('../../Figures/1_Dynamics/Trajectories/activation_q-%d.pdf'%q)

        my_plot_layout(ax = ax_b, yscale = 'log', xlabel = 'Time', ylabel = 'Clone size')
        ax_b.set_xlim(left = 3.5, right = Tf)
        ax_b.set_ylim(bottom = .8)
        fig_b.savefig('../../Figures/1_Dynamics/Trajectories/B_cells_expansion_q-%d.pdf'%q)
        fig_clones.savefig('../../Figures/1_Dynamics/Trajectories/B_cell_clones_q-%.d.pdf'%(q))

        fig.savefig('../../Figures/1_Dynamics/Trajectories/summary_1_single_trajectory_q-%d.png'%q)
    

