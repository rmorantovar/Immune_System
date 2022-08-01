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

N_ens = 1
N_r = 1e5
T0 = 0
Tf = 8
dT = 0.05
lambda_A = 6
k_pr = 0.1 # hour^-1
k_pr = k_pr*24 #days^-1

#k_pr = 0.000277
qs = [1, 2]
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
E_ms = -28

print('k_on/k_pr = %.1e'%(k_on/k_pr))

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'FTSENAYCGR'
antigen = 'TACNSEYPNTTK'
antigen = 'TACNSEYPNTTKCGRWYC'
#antigen = 'TANSEYPNTK'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'


L=len(antigen)

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
models_name = ['exponential', 'linear', ]
colors = ['darkred', 'olive']
colors_fit = ['tab:red', 'tab:olive']
growth_models = [0]#, 1]

#fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
#fig3, ax3 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})

colors_fate = [['darkred', 'tab:red'], ['olive', 'tab:olive']]

gauge_e = 0

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
#----------------------------------------------------------------

for energy_model in energy_models:
    
    for l, linear in enumerate(growth_models):

        for n_q, q in enumerate(qs):

            fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
            fig3, ax3 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
            fig4, ax4 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

            parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_q-%d_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, q, l, N_ens)+energy_model
            data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies.txt', sep = '\t', header=None)

            min_e_data = np.min(data[0])
            max_e_data = np.max(data[0])

            energies_array = np.linspace(min_e_data-1, max_e_data, 50)

            data_active = data.loc[data[1]==1]

            data_plasma = data_active.loc[data_active[2]==1]
            data_GC = data_active.loc[data_active[2]==0]

            #---- DISTRIBUTION ENERGIES ----
            
            if(linear == 0):
                fig_seq, ax_seq = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
                hist = ax_seq.hist(np.exp(data[0]), bins = np.logspace(np.log10(np.exp(min_e_data-1)), np.log10(np.exp(max_e_data)), 20), color = colors[n_q], alpha = 0, histtype = 'bar')
                counts = hist[0][np.where(hist[0]!=0)]
                energies = np.log(hist[1][np.where(hist[0]!=0)])
                popt, pcov = curve_fit(f = my_linear_func , xdata = energies[:4], ydata= np.log(counts)[:4] )
                print('beta = %.2f'%(popt[1]))
                lambd_act = popt[1]

                expfit = (np.exp(lambd_act*energies)/np.exp(lambd_act*energies[0]))*counts[0]
                ax_seq.plot(np.exp(energies), counts, color = 'indigo', alpha = 1, linestyle = '', marker = '^', ms = 10)
                #ax_seq.plot(energies, expfit, color = 'indigo', linestyle = '--', linewidth = 4, alpha = .5)
                ax_seq.plot(np.exp(energies[:8]), np.exp(my_linear_func(energies[:8], *popt)), color = 'indigo', linestyle = '--', linewidth = 4, alpha = .5, label = '%.2f'%(popt[1]))
                #ax_seq.text(x=min_e_data+1, y=1e2, s = r'$\sim e^{\beta \epsilon}$', fontsize=48, color = 'indigo')
                my_plot_layout(ax = ax_seq, xscale = 'log', yscale = 'log', xlabel = r'$K_D$', ylabel = r'$\Lambda(\epsilon)$')
                #ax_seq.set_ylim(bottom = .8)
                ax_seq.legend(title = r'$\beta$', fontsize = 30, title_fontsize = 33)
                

            #---- DISTRIBUTION ACTIVATED ENERGIES ----
            #----------------------------------------------------------------
            u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf)/N_A, Es, q, lambda_A, N_c, dE)
            QR2 = Q0*(1-np.exp(-u_on*p_a*N_c/lambda_A))
            #----------------------------------------------------------------
            #data_Kds = ax2.hist( data_active[0], bins = np.logspace(np.log10(np.exp(min_e_data-1)), np.log10(np.exp(max_e_data)), 22), color = colors_fate[n_q][0], alpha = .2, density = False)
            data_Kds = ax2.hist([np.exp(data_GC[0]), np.exp(data_plasma[0])], bins = np.logspace(np.log10(np.exp(min_e_data-1)), np.log10(np.exp(max_e_data)), 12), color = colors_fate[n_q], alpha = .6, histtype = 'barstacked', label = ['GC', 'Plasma'], density = False)
            ax2.plot(np.exp(Es[:-1]), QR2*N_r, linestyle = '--', linewidth = 2, marker = '', color = colors_fate[n_q][0])
            ax2.plot(np.exp(Es[:-1]), QR*N_r, linestyle = '-', linewidth = 2, marker = '', color = colors_fate[n_q][0])
            lambd_peak = lambdas[:-1][QR == np.max(QR)][0]
            print('beta = %.2f'%(lambd_peak))
            #ax2.vlines(np.exp(Es)[lambdas[:] < q][0], ax2.get_ylim()[0], np.max(QR), color = colors[n_q], linestyle = ':', linewidth = 2)
            my_plot_layout(ax = ax2, xscale = 'log', yscale = 'log', ticks_labelsize = 38)
            ax2.set_ylim(bottom = 1e-1)
            #ax2.legend(fontsize = 24)

            #---- SERA ----
            clone_sizes = np.exp(lambda_B*(Tf-np.array(data_plasma[3])))
            Kds = np.exp(data_plasma[0]+gauge_e)
           
            data_Kds = ax3.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data-1+gauge_e)), np.log10(np.exp(max_e_data+gauge_e)), 12), density = False, color = colors_fate[n_q][1],histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
            counts = data_Kds[0][np.where(data_Kds[0]!=0)]
            Kds_array_data = (data_Kds[1][np.where(data_Kds[0]!=0)])#+data_Kds[1][1:])/2
            ax3.plot(Kds_array_data, counts, color = colors_fate[n_q][0], alpha = .8, marker = 'o', ms = 12, linestyle = '-')
            ax3.hist(Kds, bins = np.logspace(np.log10(np.exp(min_e_data-1+gauge_e)), np.log10(np.exp(max_e_data+gauge_e)), 12), density = False, color = colors_fate[n_q][0], alpha = .2, histtype = 'step', weights = clone_sizes, zorder=0, align = 'left', linewidth = 2)
            clone_sizes_binned = np.array([])
            var_clone_sizes_binned = np.array([])
            max_clone_sizes_binned = np.array([])
            for i in np.arange(int(len(counts))):
                clone_sizes_binned = np.append(clone_sizes_binned, np.mean( np.concatenate((clone_sizes[(Kds>=data_Kds[1][i]) & (Kds<data_Kds[1][i+1])], np.array([1]) )) )  )
                #var_clone_sizes_binned = np.append(var_clone_sizes_binned, np.var(clone_sizes[(Kds<data_Kds[1][i+1]) & (Kds>data_Kds[1][i])]))
                max_clone_sizes_binned = np.append(max_clone_sizes_binned, np.max(clone_sizes[(Kds>=data_Kds[1][i]) & (Kds<data_Kds[1][i+1]) ], initial=1))
            ax4.plot(Kds_array_data, clone_sizes_binned, color = colors_fate[n_q][1], linewidth =3, linestyle = '', marker = 's', ms = 10, alpha = .4)
            print(-(q*lambda_B/lambda_A))
            ax4.plot(Kds_array_data[:-2], clone_sizes_binned[0]*(Kds_array_data[:-2]/Kds_array_data[0])**(-(q*lambda_B/lambda_A)), color = colors_fate[n_q][1], linewidth =2, linestyle = '--', marker = '', ms = 15, alpha = .8)
            #ax4.plot(Kds_array_data[:], max_clone_sizes_binned[0]*(Kds_array_data[:]/np.min(Kds_array_data))**((lambda_B/lambda_A)*(-1/q)), color = colors_fate[n_q][1], linewidth =2, linestyle = ':', marker = '', ms = 15, alpha = .8)
            
            #ax4.plot(Kds_array_data[:-5], clone_sizes_binned[0]*(Kds_array_data[:-5]/np.min(Kds_array_data))**((lambda_B/alpha)*(q)), color = colors_fate[n_q][1], linewidth =2, linestyle = '--', marker = '', ms = 15, alpha = .8)
            #ax4.plot(Kds_array_data[:], max_clone_sizes_binned[0]*(Kds_array_data[:]/np.min(Kds_array_data))**((lambda_B/alpha)*(q)), color = colors_fate[n_q][1], linewidth =2, linestyle = ':', marker = '', ms = 15, alpha = .8)
                
            if(linear==1):
                ax3.plot(Kds_array_data[:], 5e11*Kds_array_data[:]**(lambd_act), color = colors_fate[n_q][1], linewidth =3, linestyle = '--', marker = '', ms = 15)
            
            #ax3.vlines(np.exp(np.average(np.log(Kds))), 0, ax[n_q,3].get_ylim()[1], color = colors_fate[n_q][1], linewidth = 3, linestyle = '-')
            #ax3.vlines(np.exp(np.average(np.log(Kds), weights = data_bcells_active[np.where(data_active[:,2]==0)[0],-1])), 0, ax[n_q,3].get_ylim()[1], color = colors_fate[n_q][1], linewidth = 3, linestyle = '--')
            my_plot_layout(ax = ax3, yscale = 'log', xscale = 'log', xlabel = r'$K_D$', ylabel = 'counts', ticks_labelsize = 38)
            my_plot_layout(ax = ax4, yscale = 'log', xscale = 'log', xlabel = r'$K_D$', ylabel = 'Clone size')
            #ax3.legend(fontsize = 24, loc = 4)

            fig_seq.savefig('../../Figures/1_Dynamics/Trajectories/Sequences_expansion_q-%d.pdf'%q)
            fig2.savefig("../../Figures/1_Dynamics/Trajectories/Plasma_vs_GC_q-%d.pdf"%(q))
            fig3.savefig("../../Figures/1_Dynamics/Trajectories/Sera_q-%d.pdf"%(q))
            fig4.savefig("../../Figures/1_Dynamics/Trajectories/Sera_2_q-%d.pdf"%(q))


