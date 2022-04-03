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
style.use('seaborn-paper')

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

N_r = 1e4
T0 = 0
Tf = 35
dT = 0.001
alpha = 1
alpha_lin = alpha*1e8
gamma = 0.0333
gamma = 0.000277
beta = 0.5
k_on = 1e6 #(M*s)^-1
N_c = 1e0
e_MS = -27


antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'FTSENAYCGR'
antigen = 'TACNSEYPNTTK'

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

gauge_e = 0

lambd = 0.5976041516921407
lambd = 0.8
lambd = 1.2

Tmin = .01
Tmax = 50

antigen_list = [i for i in antigen]
antigen_seq = np.array([], dtype = int)
for i, aa in enumerate(antigen_list):
    index = Alphabet_list.index(aa)
    antigen_seq = np.append(antigen_seq, int(index))
PWM_data = M2[:,antigen_seq]

#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

Ts = np.linspace(Tmin, Tmax, 20000)
lambdas = 1/Ts[:-1]
F_PWM = -Ts*np.log(Z_PWM(PWM_data, Ts))
Us = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts)) + e_MS
print(Us[0])
dU = np.diff(Us)

for energy_model in energy_models:
    
    fig, ax = plt.subplots(2,4,figsize=(40,18), gridspec_kw={'hspace':0.25})

    for j, linear in enumerate(growth_models):

        fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
        fig3, ax3 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})

        colors_activation = []

        parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_linear-%d_'%(alpha, beta, gamma, linear)+energy_model
        #data_bcells = np.loadtxt(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/bcells.txt')
        data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies.txt', sep = '\t', header=None)

        N_data = len(data[0])/10

        min_e_data = np.min(data[0])
        max_e_data = np.max(data[0])

        print(min_e_data, max_e_data)

        energies_array = np.linspace(min_e_data-1, max_e_data, 50)

        data_active = data.loc[data[1]==1]

        data_plasma = data_active.loc[data_active[2]==1]
        data_GC = data_active.loc[data_active[2]==0]


        #print(np.min(data))

        #---- DISTRIBUTION ENERGIES ----
        hist = ax[j,0].hist(data[0], bins = np.linspace(min_e_data-1, max_e_data+1, 10), color = colors[j], alpha = .8, histtype = 'bar')
        my_plot_layout(ax = ax[j,0], yscale = 'log', xlabel = 'Energies', ylabel = 'counts')
        ax[j,0].set_ylim(bottom = .8)
        if(linear == 0):
            fig_seq, ax_seq = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
            expfit = np.exp(lambd*(energies_array-(min_e_data-1)))
            ax_seq.plot(energies_array, expfit, color = 'indigo', linestyle = '--', linewidth = 4, alpha = .5)
            #hist = np.histogram(data[:,0], bins = np.linspace(min_e_data-1, max_e_data, 10))
            ax_seq.plot(hist[1][:-1][:], hist[0][:], color = 'indigo', alpha = 1, linestyle = '', marker = '^', ms = 10)
            #ax_seq.text(x=min_e_data+1, y=1e2, s = r'$\sim e^{\lambda \epsilon}$', fontsize=48, color = 'indigo')
            my_plot_layout(ax = ax_seq, yscale = 'log', xlabel = r'$\epsilon$', ylabel = r'$\Lambda(\epsilon)$')
            ax_seq.set_ylim(bottom = .8)
            fig_seq.savefig('../../Figures/1_Dynamics/Trajectories/Sequences_expansion_'+energy_model+'.pdf')

        #---- DISTRIBUTION ACTIVATED ENERGIES ----
        Omega = 20**L
        S = np.cumsum(lambdas[:-1]*dU)
        P0 = 1/Omega
        Q0 = np.exp(S)*P0
        p_a = (gamma/(gamma + (k_on*np.exp(Us[:-1]))**2 ) )
        p_e = ((np.exp(alpha*Tf)/N_A)*k_on)/(alpha/(24*60*60)) #probability of engagement up to time t
        Q_R = N_r*Q0*(1-np.exp(-p_e*p_a*N_c))
        

        ax[j,1].hist(data_active[0], bins = np.linspace(min_e_data-1, max_e_data+1, 10), color = colors[j], alpha = .8, histtype = 'bar')
        ax[j,1].plot(Us[:-1], Q_R, linestyle = '-', linewidth = 2, color = colors[j])
        my_plot_layout(ax = ax[j,1], yscale = 'log', xlabel = 'Energies', ylabel = 'counts')
        ax[j,1].set_ylim(bottom = ax[j,0].get_ylim()[0], top=ax[j,0].get_ylim()[1])

        #---- DISTRIBUTION GC VS PLASMA ----
        #ax[j,2].hist([data_active[np.where(data_active[:,2]==1)[0],0], data_active[np.where(data_active[:,2]==0)[0],0]], bins = np.linspace(min_e, max_e, 20), color = colors_fate[j], alpha = .6, histtype = 'barstacked', label = ['GC', 'Plasma'])
        ax[j,2].hist(data_GC[0], bins = np.linspace(min_e_data-1, max_e_data+1, 15), color = colors_fate[j][0], alpha = .6, histtype = 'barstacked', label = ['GC'])
        my_plot_layout(ax = ax[j,2], yscale = 'linear', xlabel = 'Energies', ylabel = 'counts')
        #ax[j,2].set_ylim(bottom = ax[j,0].get_ylim()[0], top=ax[j,0].get_ylim()[1])
        ax[j,2].set_ylim(bottom = ax[j,2].get_ylim()[0])
        ax[j,2].legend(fontsize = 24)
        ax2.hist([data_GC[0], data_plasma[0]], bins = np.linspace(min_e_data-1, max_e_data+1, 15), color = colors_fate[j], alpha = .6, histtype = 'barstacked', label = ['GC', 'Plasma'])
        ax2.plot(Us[:-1], Q_R, linestyle = '-', linewidth = 2, color = colors[j])
        my_plot_layout(ax = ax2, yscale = 'linear', xlabel = 'Energies', ylabel = 'counts')
        ax2.set_ylim(bottom = ax2.get_ylim()[0], top=ax2.get_ylim()[1])
        ax2.legend(fontsize = 24)

        #---- SERA ----
        clone_sizes = np.exp(beta*(Tf-np.array(data_plasma[3])))
        Kds = np.exp(data_plasma[0]+gauge_e)
        data_Kds = ax[j,3].hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data-1+gauge_e)), np.log10(np.exp(max_e_data+1+gauge_e)), 8), density = False, color = colors_fate[j][1], alpha = .5, histtype = 'bar', zorder=10)
        Kds_array_data = (data_Kds[1][:-1]+data_Kds[1][1:])/2
        Kds_array = np.logspace(np.log10(np.exp(min_e_data-1+gauge_e)), np.log10(np.exp(max_e_data+1+gauge_e)), 50)
        Kds_array = (Kds_array[:-1]+Kds_array[1:])/2
        ax[j,3].hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data-1+gauge_e)), np.log10(np.exp(max_e_data+1+gauge_e)), 8), density = False, color = 'silver', alpha = .4, histtype = 'bar', weights = clone_sizes, zorder=0)
        clone_sizes_binned = np.array([])
        var_clone_sizes_binned = np.array([])
        max_clone_sizes_binned = np.array([])
        for i in np.arange(int(len(data_Kds[0]))):
            clone_sizes_binned = np.append(clone_sizes_binned, np.mean( np.concatenate((clone_sizes[(Kds<=data_Kds[1][i+1]) & (Kds>data_Kds[1][i])], np.array([1]) )) )  )
            #var_clone_sizes_binned = np.append(var_clone_sizes_binned, np.var(clone_sizes[(Kds<data_Kds[1][i+1]) & (Kds>data_Kds[1][i])]))
            max_clone_sizes_binned = np.append(max_clone_sizes_binned, np.max(clone_sizes[(Kds<data_Kds[1][i+1]) & (Kds>data_Kds[1][i])],initial=1))
        ax[j,3].plot(Kds_array_data, clone_sizes_binned, color = 'dimgray', linewidth =3, linestyle = '', marker = 's', ms = 12, alpha = .8)
        ax[j,3].plot(Kds_array_data, max_clone_sizes_binned, color = 'dimgray', linewidth =3, linestyle = '', marker = '*', ms = 10, alpha = .8)
        #ax[j,3].errorbar(x=Kds_array_data, y=clone_sizes_binned, yerr = np.sqrt(var_clone_sizes_binned), capsize = 10 , linestyle = '', color = 'dimgray', linewidth =2)
        if(linear==0):
        	ax[j,3].plot(Kds_array[:], (Kds_array[:]/np.min(Kds_array))**(lambd), color = colors_fate[j][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
        	#ax[j,3].plot(Kds_array[:], np.exp(beta*Tf)*np.exp(beta*(gauge_e-25)/alpha)*N_A**(-beta/alpha)*(Kds_array[:]/np.min(Kds_array))**(0), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
        	#ax[j,3].plot(Kds_array[:], np.exp(beta*Tf)*np.exp(beta*(gauge_e-25)/alpha)*N_A**(-beta/alpha)*(Kds_array[:]/np.min(Kds_array))**(lambd+0), color = 'silver', linewidth =3, linestyle = '--', marker = '', ms = 15)
        if(linear==1):
            ax[j,3].plot(Kds_array[:], (Kds_array[:]/np.min(Kds_array))**(lambd), color = colors_fate[j][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
            #ax[j,3].plot(Kds_array[:-17], np.exp(beta*Tf)*np.exp(beta*1/(alpha_lin))*np.exp(-beta/(alpha_lin)*np.exp(-(gauge_e-25))*N_A*Kds_array[:-17]), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
            ax[j,3].plot(Kds_array[:-17], 5e11*Kds_array[:-17]**(lambd)*np.exp(beta*Tf)*np.exp(beta*1/(alpha_lin))*np.exp(-beta/(alpha_lin)*np.exp(-(gauge_e-25))*N_A*Kds_array[:-17]), color = 'silver', linewidth =3, linestyle = '--', marker = '', ms = 15)

        ax[j,3].vlines(np.exp(np.average(np.log(Kds))), 0, ax[j,3].get_ylim()[1], color = colors_fate[j][1], linewidth = 3, linestyle = '-')
        #ax[j,3].vlines(np.exp(np.average(np.log(Kds), weights = data_bcells_active[np.where(data_active[:,2]==0)[0],-1])), 0, ax[j,3].get_ylim()[1], color = colors_fate[j][1], linewidth = 3, linestyle = '--')
        my_plot_layout(ax = ax[j,3], yscale = 'log', xscale = 'log', xlabel = r'$K_D$', ylabel = 'counts')
        #ax[j,3].legend(loc = 4, fontsize = 24)

        ax3.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data-1+gauge_e)), np.log10(np.exp(max_e_data+1+gauge_e)), 8), density = False, color = colors_fate[j][1], alpha = .5, histtype = 'bar', zorder=10)
        ax3.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data-1+gauge_e)), np.log10(np.exp(max_e_data+1+gauge_e)), 8), density = False, color = 'silver', alpha = .4, histtype = 'bar', weights = clone_sizes, zorder=0)
        ax3.plot(Kds_array_data, clone_sizes_binned, color = 'dimgray', linewidth =3, linestyle = '', marker = '*', ms = 12, alpha = .8)
        #ax3.plot(Kds_array_data, max_clone_sizes_binned, color = colors_fate[j][1], linewidth =3, linestyle = '', marker = '*', ms = 14, alpha = .8)
        #ax3.errorbar(x=Kds_array_data, y=clone_sizes_binned, yerr = np.sqrt(var_clone_sizes_binned) , capsize = 10, linestyle = '', color = 'dimgray', linewidth =2)

        if(linear==0):
            ax3.plot(Kds_array[:], (Kds_array[:]/np.min(Kds_array))**(lambd), color = colors_fate[j][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
            #ax3.plot(Kds_array[:], 5e2*np.exp(beta*Tf)*np.exp(beta*(gauge_e-25)/alpha)*N_A**(-beta/alpha)*Kds_array[:]**(0), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
            ax3.plot(Kds_array[:], np.exp(beta*Tf)*(lambd*1e8/N_data)**(-beta/alpha)*(Kds_array[:]/np.min(Kds_array))**(0), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
            ax3.plot(Kds_array[:], np.exp(beta*Tf)*(lambd*1e8/N_data)**(-beta/alpha)*(Kds_array[:]/np.min(Kds_array))**(beta*lambd/alpha), color = colors_fate[j][1], linewidth = 3, linestyle = ':', marker = '', ms = 15)
            ax3.plot(Kds_array[:], np.exp(beta*Tf)*(lambd*1e8/N_data)**(-beta/alpha)*(Kds_array[:]/np.min(Kds_array))**(0 + lambd), color = 'silver', linewidth =3, linestyle = '--', marker = '', ms = 15)
        if(linear==1):
            ax3.plot(Kds_array[:], 5e11*Kds_array[:]**(lambd), color = colors_fate[j][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
            #ax3.plot(Kds_array[:-17], np.exp(beta*Tf)*np.exp(beta*1/(alpha_lin))*np.exp(-beta/(alpha_lin)*np.exp(-(gauge_e-25))*N_A*Kds_array[:-17]), color = 'dimgray', linewidth =3, linestyle = '--', marker = '', ms = 15)
            #ax3.plot(Kds_array[:-17], 1e7*Kds_array[:-17]**(lambd)*np.exp(beta*Tf)*np.exp(beta*1/(alpha_lin))*np.exp(-beta/(alpha_lin)*np.exp(-(gauge_e-25))*N_A*Kds_array[:-17]), color = 'silver', linewidth =3, linestyle = '--', marker = '', ms = 15)
        #ax3.vlines(np.exp(np.average(np.log(Kds))), 0, ax[j,3].get_ylim()[1], color = colors_fate[j][1], linewidth = 3, linestyle = '-')
        #ax3.vlines(np.exp(np.average(np.log(Kds), weights = data_bcells_active[np.where(data_active[:,2]==0)[0],-1])), 0, ax[j,3].get_ylim()[1], color = colors_fate[j][1], linewidth = 3, linestyle = '--')
        my_plot_layout(ax = ax3, yscale = 'log', xscale = 'log', xlabel = r'$K_D$', ylabel = 'counts')
        #ax3.legend(fontsize = 24, loc = 4)

        fig2.savefig("../../Figures/1_Dynamics/Trajectories/Plasma_vs_GC_linear-%d.pdf"%(linear))
        fig3.savefig("../../Figures/1_Dynamics/Trajectories/Sera_linear-%d.pdf"%(linear))


    fig.savefig('../../Figures/1_Dynamics/Trajectories/summary_2_single_trajectory_'+energy_model+'.png')
    
