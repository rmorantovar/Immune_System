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

Matrix = 'BLOSUM62'
Matrix = 'MJ2'
#Matrix = 'MM'

fig, ax = plt.subplots(1, 4, figsize=(40,8))

#-------------------------MATRIX----------------------------------
#w, v = LA.eig(M2_list)
if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    M2_list = M2.tolist()
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'MJ Model', ax = ax[0])
if(Matrix == 'MM'):
    M2 = (np.loadtxt(Text_files_path + Matrix + '.txt')+1)*e0
    M2_list = M2.tolist()
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'MJ Model', ax = ax[0])
if(Matrix == 'BLOSUM62'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,25))
    M2_list = M2.tolist()
    Alphabet = np.array(['A', 'R'  ,'N' , 'D'  ,'C' , 'Q'  ,'E'  ,'G'  ,'H' , 'I'  ,'L'  ,'K'  ,'M' , 'F' , 'P' , 'S'  ,'T' , 'W' , 'Y' , 'V' , 'B' , 'Z'  ,'X',  '*'])
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'BLOSUM62 Model', ax = ax[0])
L_alphabet = len(Alphabet)

#-------------------------PWM----------------------------------
#antigen_seq = np.random.randint(0, len(Alphabet), L)
#antigen = Alphabet[antigen_seq]
antigen = 'FMLFMAVFVMTSWYC'
L=len(antigen)
antigen_list = [i for i in antigen]
contributions = np.zeros(shape = (1,20))
antigen_seq = np.array([], dtype = int)
for i, aa in enumerate(antigen_list):
    index = Alphabet_list.index(aa)
    antigen_seq = np.append(antigen_seq, int(index))
print('Antige:', antigen)
antigen_seq, antigen
PWM = M2[:,antigen_seq]
#for i in np.arange(L):
#    PWM[:,i]-=np.min(PWM[:,i], axis=0)
PWM_list = PWM.tolist()
plot_PWM(PWM=PWM_list, Alphabet=Alphabet, sequence = antigen, title=r'PWM', ax = ax[1])
min_E = np.sum([np.min(PWM[:,i]) for i in range(len(PWM[0,:]))])
avg_E = np.sum([np.mean(PWM[:,i]) for i in range(len(PWM[0,:]))])
var_E = np.sum([np.var(PWM[:,i]) for i in range(len(PWM[0,:]))])
max_E = np.sum([np.max(PWM[:,i]) for i in range(len(PWM[0,:]))])

d_eff = 4
mean_lowest_effect=0
for i, index in enumerate(antigen_seq):
    ax[2].scatter(np.ones(d_eff)*i, np.sort(PWM[:,i])[:d_eff])
    mean_lowest_effect+=np.mean(np.sort(PWM[:,i])[1:d_eff]-np.sort(PWM[:,i])[0])
mean_lowest_effect/=L

mean_lowest_effects=np.array(20)
for d_i in [5]:
    mean_lowest_effect_i = 0
    for i, index in enumerate(antigen_seq):
        mean_lowest_effect_i+=np.mean(np.sort(PWM[:,i])[1:d_i])
    mean_lowest_effect_i/=L
    mean_lowest_effects = np.append(mean_lowest_effects, mean_lowest_effect_i)
    
ax[2].hlines(mean_lowest_effect+np.min(PWM), ax[2].get_xlim()[0], ax[2].get_xlim()[1])
ax[2].set_xticks(range(L))
ax[2].set_xticklabels(antigen)
ax[2].tick_params(labelsize = 22)

linear_PWM = np.reshape(PWM, (L*20,1))
energies = np.linspace(0, np.max(linear_PWM), 20)
ds_eff = np.zeros_like(energies)
for i, energy in enumerate(energies[:-1]):
    ds_eff[i] = np.size(np.where((linear_PWM<energies[i+1]) & (linear_PWM>energies[i])))/L

ax[3].hist(linear_PWM, bins = 50, alpha = .5, align='left')
ax[3].vlines(mean_lowest_effect, ax[3].get_ylim()[0], ax[3].get_ylim()[1], color = 'darkblue', linestyle = 'dashed')
my_plot_layout(ax=ax[3])

fig.savefig('../../Figures/5_Geometric_exponent/Matrix_and_PWM.png')

#-------------------------RELATIVE CLONE SIZE _ N----------------------------------
T=.5
T0 = 0
Tf = 25
dT = 0.5
alpha = 1
gamma = 0
beta = 0.5
Ns = np.array([2e2, 2e3, 2e4, 2e5, 2e6])
markers = ['o', '^', '*', 's', 'X']
exponents_rcs = np.array([])
lambda_simulations = np.array([])
vars_rcs = np.array([])
energy_model='MJ'
linear = 0
fig, ax = plt.subplots(figsize = (10,8))
for n, N in enumerate(Ns):
    # ------- from relative clone size -------
    data_bcells_ensemble = np.loadtxt(Text_files_path + 'Dynamics/Ensemble/bcells_ensemble_L-%d_N-%d_Antigen-'%(L, N)+antigen+'_alpha-%.6f_beta-%.6f_gamma-%.6f_linear-%d_'%(alpha, beta, gamma, linear)+energy_model+'.txt')
    N_final_active = np.loadtxt(Text_files_path + "Dynamics/Ensemble/N_final_active_L-%d_N-%d_Antigen-"%(L, N)+antigen+"_alpha-%.6f_beta-%.6f_gamma-%.6f_Linear-%d_"%(alpha, beta, gamma, linear)+energy_model+".txt")
    N_final_active = np.concatenate(([0], N_final_active))
    N_final_active_cum = np.cumsum(N_final_active)
    N_clones = 20
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
    lambda_simulations = np.append(lambda_simulations, -beta/(alpha*(popt_csd[1])))
    vars_rcs = np.append(vars_rcs, pcov_csd[1,1])
    ax.plot(n_array, np.exp(my_linear_func(np.log(n_array), *popt_csd)), linestyle = '-', marker = '', ms = '10', linewidth = 2, color=plt.cm.Set1(n), alpha= .8)

    # ------- from density of sequences -------
    #data_dynamics_tail=np.loadtxt(Text_files_path + 'Dynamics/ensemble/energies_tail_ensemble_L-%d_N-%d_Antigen-'%(L, N_dynamics)+antigen+'.txt')
    #data_dynamics_tail_rho = np.histogram(data_dynamics_tail[:,0], bins = 8, density = False) # Creating histograms
    #e_dynamics_tail = data_dynamics_tail_rho[1][:-1] + abs(data_dynamics_tail_rho[1][1:]-data_dynamics_tail_rho[1][:-1])/2
    #rho_e_dynamics_tail = data_dynamics_tail_rho[0]*1e8
    #e_lineal = e_dynamics_tail[:4]
    #rho_e_lineal = rho_e_dynamics_tail[:4]
    #popt_ds, pcov_ds = curve_fit(my_linear_func, e_lineal, np.log(rho_e_lineal), p0 = (60,1))
    
    #exponents_ds = np.append(exponents_ds, (popt_ds[1]*alpha)/beta)
    #vars_ds = np.append(vars_ds, pcov_ds[1,1])
    
my_plot_layout(ax = ax, xscale='log', yscale= 'log', xlabel=r'Clone $n$', ylabel=r'Relative clone size $x_n$',
                  ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax.legend(fontsize = 22, loc = 0)
fig.savefig('../../Figures/1_Dynamics/relative_clone_size_n_N.png')


#-------------------------lambda _ N----------------------------------

es = np.linspace(min_E, max_E, 10000)
de = es[1]-es[0]

def P_e_gaussian(avg_E, var_E, es):
    return (2*np.pi*var_E)**(-0.5)*np.exp(-(es-avg_E)**2/(2*var_E))
def P_min_e(N, es):
    return (N*(1-np.cumsum(P_e_gaussian(avg_E, var_E, es)*de))**(N-1)*(P_e_gaussian(avg_E, var_E, es)))
fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18})
lambd_gaussian = np.array([])
Ns_array = np.logspace(2, 7, 50)
for N in Ns_array:
    p_min_e = -np.diff((1-np.cumsum(P_e_gaussian(avg_E, var_E, es)*de))**N)/np.diff(es)
    avg_min_E = np.sum(es[:-1]*p_min_e*de)
    lambd_gaussian = np.append(lambd_gaussian, (avg_E-avg_min_E)/var_E)

ax.plot(Ns_array, lambd_gaussian, label = 'Gaussian', linestyle = '--', linewidth= 2, color = 'indigo')

for n, N in enumerate(Ns):
    #ax.vlines(N/10, 0, lambd_gaussian[np.where(Ns_array<N)][-1], color = 'black', linestyle = '--', linewidth = 1, alpha = .4)  
    ax.scatter(N, lambda_simulations[n], color = 'indigo', s = 60)
    ax.errorbar(x = N, y = lambda_simulations[n], yerr = 1.96*np.sqrt(vars_rcs[n]), color = 'indigo', linestyle = '', capsize = 4)

my_plot_layout(ax=ax, xscale = 'log', yscale = 'linear', ylabel = '$\lambda(N)$', xlabel = '$N$')
ax.legend(fontsize=24)
fig.savefig('../../Figures/5_Geometric_exponent/lambda_simulations_N_prueba.png')

end = time.time()
print('\nFinished in', '%.2f'%((end-start)/60), 'minutes')





