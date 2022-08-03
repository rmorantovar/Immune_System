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
N_r = 5e5
N_r = 1e8
T0 = 0
Tf = 8
Tf = 6
dT = 0.01
lambda_A = 6
k_pr = 0.1 # hour^-1
k_pr = k_pr*24 #days^-1

#k_pr= 0.000277
thetas = [2, 1.5, 1]
#thetas = [2, 1.5]
colors_q = ['darkblue', 'olive', 'darkred']
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
E_ms = -27.63

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'FTSENAYCGR'
antigen = 'TACNSEYPNTTK'
antigen = 'TACNSEYPNTTKCGRWYC'
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

Es, dE, Q0, lambdas = calculate_Q0(0.01, 50, 200000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

beta_r = lambdas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
Kd_r = np.exp(E_r)
#----------------------------------------------------------------

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
models_name = ['exponential']#, 'linear',]
colors = ['tab:blue', 'tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = [0]

lambd = 1.4
for i_theta, theta in enumerate(thetas):

    beta_theta = lambdas[lambdas>theta][-1]
    E_theta = Es[lambdas>theta][-1]
    Kd_theta = np.exp(E_theta)
    for energy_model in energy_models:
        
        fig_b, ax_b = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
        fig_a, ax_a = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
        fig_act, ax_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
        fig_clones, ax_clones = plt.subplots(1, 3, figsize=(30,10), gridspec_kw={'left':0.06, 'right':.98, 'bottom':.1, 'top': 0.9})
        for j, linear in enumerate(growth_models):

            colors_activation = []

            parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, theta, j, N_ens)+energy_model
            data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies.txt', sep = '\t', header=None)

            min_e_data = np.min(data[0])
            max_e_data = np.max(data[0])

            data_active = data.loc[data[1]==1]

            activations_times = np.array(data_active[3])
            energies  = np.array(data_active[0])

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

            #---- B cell linages ----
            expfit = 1e-5*np.exp(lambda_B*time)
            clone_sizes = np.ones((int(data_N_active_linages[-1]), len(time)))
            for i in np.arange(0, int(data_N_active_linages[-1])):
                clone_sizes[i, int(activations_times[i]/dT):] = np.exp(lambda_B*(time[int(activations_times[i]/dT):] - activations_times[i] ))
                if(i%100==0):
                    ax_b.plot(time, clone_sizes[i,:], color = colors[j], linewidth = 1, linestyle  = '-', marker = '', ms = 5, alpha = .6);
            if(linear == 0):
                #ax_b.plot(time, expfit, color = 'indigo', linestyle = '--', linewidth = 4)
                ax_b.text(x=5, y=3e3, s = r'$\sim e^{\lambda_B t}$', fontsize=48, color = 'indigo')

            #---- Activation rate ----
            #----------------------------------------------------------------
            u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf)/N_A, Es, theta, lambda_A, N_c, dE)
            psi = ((k_on*N_c)/(N_A))*p_a
            M_r = N_r*N_c*np.sum(Q0*p_a*dE)
            print(M_r)
            #----------------------------------------------------------------
            ax_a.plot(time[:-1], data_N_active_linages[:], linestyle = '', marker = 'o', ms = 5, linewidth = 2, label = 'simulation', color = colors[j], alpha = .6)
            if(linear == 0):
                m_bar = np.array([N_r*(1-np.sum(np.exp(-(psi/lambda_A)*(np.exp(lambda_A*t)-1))*Q0*dE)) for t in time])
                m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))
                ax_a.plot(time, np.ones_like(time), color = 'black', linestyle = '-', linewidth = 1)
                ax_a.plot(time, m_bar, color = colors[j], label = models_name[j] + ' growth', linestyle = '--', linewidth = 3)
                ax_a.plot(time, m_bar_approx, color = 'grey', label = models_name[j] + ' growth', linestyle = ':', linewidth = 3)
            if(linear == 1):
                theory = ((k_on*M_r)/(N_A))*(1e8*time + 0.5*lambda_A*time**2)
            #ax_a.set_ylim(top=.5*N_r)


            #---- Entropy ----
            size_lim = 0
            if theta==1:
                size_lim=1000
            if theta==2:
                size_lim=1

            #energies = energies[np.where((clone_sizes[:,-1]>=(np.min(clone_sizes[:,-1])*size_lim)))[0]]
            #activations_times = activations_times[np.where((clone_sizes[:,-1]>=(np.min(clone_sizes[:,-1])*size_lim)))[0]]
            #clone_sizes = clone_sizes[np.where((clone_sizes[:,-1]>=(np.min(clone_sizes[:,-1])*size_lim)))[0],:]

            total_pop = np.sum(clone_sizes, axis = 0)
            bcell_freqs = clone_sizes/total_pop
            entropy = -np.sum(bcell_freqs*np.log(bcell_freqs), axis = 0)

            #---- Total Pop size ----

            #---- Relative clone sizes ----

            #---- Stackplots ----
            #ax_clones.stackplot(time[400:], bcell_freqs[np.where(bcell_freqs[:, -1]>.01)[0], 400:]);
            # colors = []
            # min_bell_freq = np.min(bcell_freqs[:,-1])
            # for c in range(int(len(clone_sizes[:,0]))):
            #     if bcell_freqs[c, -1]>(50*min_bell_freq/theta):
            #         if theta==1:
            #             colors.append('indianred')
            #         if theta==2:
            #             colors.append('dodgerblue')
            #     else:
            #         colors.append('silver')

            days_plot = np.linspace(Tf-1.5, Tf, 3)
            filter_size = clone_sizes[:, -1].argsort()
            filter_factor = 10**(2.2*(2-theta))
            bcell_freqs = bcell_freqs[filter_size, :][-int(100*filter_factor):, :]
            energies = energies[filter_size][-int(100*filter_factor):]
            activations_times = activations_times[filter_size][-int(100*filter_factor):]

            #positions = np.random.random((len(energies), 2))
            angles = np.random.random(len(energies))*2*np.pi
            energies = energies - (E_ms)
            radious = ((energies/(30)))*4

            for i_plot in range(len(days_plot)):
                for i_c in range(len(energies)):
                    if (activations_times[i_c]<=days_plot[i_plot]):
                        #print(i_c, positions[i_c], int(days_plot[i_plot]*len(time)/8)-1, time[int(days_plot[i_plot]*len(time)/8)-1], clone_sizes[i_c, int(days_plot[i_plot]*len(time)/8)-1], activations_times[i_c])
                        #circle = plt.Circle(positions[i_c], np.sqrt(bcell_freqs[i_c, int(days_plot[i_plot]*len(time)/Tf)-1]/(np.pi*Tf)), color = colors_q[i_theta], alpha = 1-(energies[i_c]/np.max(energies)))
                        circle = plt.Circle((radious[i_c]*np.cos(angles[i_c]), radious[i_c]*np.sin(angles[i_c])), np.sqrt(bcell_freqs[i_c, int(days_plot[i_plot]*len(time)/Tf)-1]/(np.pi)), color = colors_q[i_theta], alpha = .8)
                        ax_clones[i_plot].add_patch(circle)
                circle = plt.Circle((0, 0), 3, edgecolor = colors_q[i_theta], facecolor = 'none')
                ax_clones[i_plot].add_patch(circle)
                circle = plt.Circle((0, 0), 4*((E_r-E_ms)/30), edgecolor = 'grey', facecolor = 'none', linestyle = 'dashed', linewidth = 4)
                ax_clones[i_plot].add_patch(circle)
                circle = plt.Circle((0, 0), 4*((E_theta-E_ms)/30), edgecolor = 'grey', facecolor = 'none', linestyle = 'dotted', linewidth = 4)
                ax_clones[i_plot].add_patch(circle)

        for i_plot in range(len(days_plot)): 
            my_plot_layout(ax = ax_clones[i_plot], ticks_labelsize=34)
            ax_clones[i_plot].set_xlim(-3.01, 3.01)
            ax_clones[i_plot].set_ylim(-3.01, 3.01)
            ax_clones[i_plot].set_xticks([])
            ax_clones[i_plot].set_yticks([])
        
        my_plot_layout(ax = ax_a, yscale = 'log', xlabel = 'Time', ylabel = r'$\bar m$')
        ax_a.set_xlim(left = 3.5, right = Tf)
        ax_a.set_ylim(top=2*N_r, bottom = 1e-6)
        fig_a.savefig('../../Figures/1_Dynamics/Trajectories/activation_theta-%.1f.pdf'%theta)

        my_plot_layout(ax = ax_b, yscale = 'log', xlabel = 'Time', ylabel = 'Clone size')
        ax_b.set_xlim(left = 3.5, right = Tf)
        ax_b.set_ylim(bottom = .8)
        fig_b.savefig('../../Figures/1_Dynamics/Trajectories/B_cells_expansion_theta-%.1f.pdf'%theta)
        fig_clones.savefig('../../Figures/1_Dynamics/Trajectories/B_cell_clones_theta-%.1f.pdf'%(theta), dpi = 10)
    

