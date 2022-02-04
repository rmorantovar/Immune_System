import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
from Immuno_models import*
import scipy.special as sc
import pickle
from matplotlib import style
from scipy.optimize import curve_fit
import time

start = time.time()
print('Starting...\n')

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

N_A = 6.02214076e23
k_BT = 1.380649e-23*293
style.use('seaborn-paper')
#antigen = 'TACNSFVMTSATNLFSEYPN'
#antigen = 'FMLFMAVFVMTSWYC'
antigen = 'TACNSEYPNTTK'
#antigen = 'NTKTAATNLF'
antigens = ['NTKTAATNLF', 'TACNSEYPNTTK', 'FMLFMAVFVMTSWYC', 'TACNSFVMTSATNLFSEYPN', 'TACNSFVMTSSEYPNPNEFYMAWYC']
antigens = ['TACNSFVMTSSEYPNPNEFYMAWYC']

#-------------------------RELATIVE CLONE SIZE _ N----------------------------------
for antigen in antigens:
    L=len(antigen)
    T=.5
    T0 = 0
    Tf = 25
    dT = 0.5
    alpha = 1
    gamma = 0
    beta = 0.5
    Ns = np.array([2e2, 2e3, 2e4, 2e5, 2e6])
    Ns = np.array([2e2, 2e3, 2e4, 2e5])
    markers = ['o', '^', '*', 's', 'X']
    exponents_rcs = np.array([])
    lambdas_simulation = np.array([])
    vars_rcs = np.array([])
    energy_model='MJ'
    linear = 0
    growth_models = ['exponential', 'linear']

    lambda_Ns_file = open(Text_files_path+'Dynamics/Ensemble/lambda_Ns_alpha-%.2f_beta-%.2f_L-%d_'%(alpha, beta, L)+growth_models[linear]+'_'+energy_model+'.pkl','wb')

    fig, ax = plt.subplots(figsize = (10,8))

    for n, N in enumerate(Ns):
        # ------- from relative clone size -------
        data_bcells_ensemble = np.loadtxt(Text_files_path + 'Dynamics/Ensemble/bcells_ensemble_L-%d_N-%d_Antigen-'%(L, N)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_linear-%d_'%(alpha, beta, gamma, linear)+energy_model+'.txt')
        N_final_active = np.loadtxt(Text_files_path + "Dynamics/Ensemble/N_final_active_L-%d_N-%d_Antigen-"%(L, N)+antigen+"_alpha-%.6f_beta-%.6f_gamma-%.6f_Linear-%d_"%(alpha, beta, gamma, linear)+energy_model+".txt")
        N_final_active = np.concatenate(([0], N_final_active))
        N_final_active_cum = np.cumsum(N_final_active)
        N_clones = 10 #Maybe reduce this so that the exponent is closer to \epsilon_m
        Clone_relative_sizes = np.zeros(N_clones)
        N_ens = 0
        for j, N_final in enumerate(N_final_active_cum[:-1]):
            if(N_final_active[j+1]>=N_clones):
                N_ens +=1
                temp_array = np.flip(np.sort(data_bcells_ensemble[int(N_final):int(N_final+N_final_active[j+1])]))
                for k in range(N_clones):
                    Clone_relative_sizes[k] += (temp_array[k]/temp_array[0])
                    
        n_array = np.linspace(1, N_clones, N_clones)
        ax.plot(n_array, Clone_relative_sizes/N_ens, marker = markers[n], linestyle = '', ms = 11, linewidth = 4, alpha = .5, label = 'N=%.e'%(N), color = plt.cm.Set1(n))
        popt_csd, pcov_csd = curve_fit(my_linear_func, np.log(n_array), np.log(Clone_relative_sizes/N_ens))
        print(-beta/(alpha*(popt_csd[1])))
        exponents_rcs = np.append(exponents_rcs, popt_csd[1])
        lambdas_simulation = np.append(lambdas_simulation, -beta/(alpha*(popt_csd[1])))
        vars_rcs = np.append(vars_rcs, pcov_csd[1,1])
        ax.plot(n_array, np.exp(my_linear_func(np.log(n_array), *popt_csd)), linestyle = '-', marker = '', ms = '10', linewidth = 2, color=plt.cm.Set1(n), alpha= .8)

    pickle.dump(np.array([lambdas_simulation, vars_rcs], dtype=object), lambda_Ns_file)
    lambda_Ns_file.close()

    my_plot_layout(ax = ax, xscale='log', yscale= 'log', xlabel=r'Clone $n$', ylabel=r'Relative clone size $x_n$',
                      ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax.legend(fontsize = 22, loc = 0)
    fig.savefig('../../Figures/1_Dynamics/relative_clone_size_n_N_L-%d.png'%(L))



