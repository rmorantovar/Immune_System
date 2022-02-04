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

#Matrix = 'BLOSUM62'
Matrix = 'MJ2'
#Matrix = 'MM'

fig_M, ax_M = plt.subplots(1, 4, figsize=(40,8))

#-------------------------MATRIX----------------------------------
#w, v = LA.eig(M2_list)
if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    M2_list = M2.tolist()
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'MJ Model', ax = ax_M[0])
if(Matrix == 'MM'):
    M2 = (np.loadtxt(Text_files_path + Matrix + '.txt')+1)*e0
    M2_list = M2.tolist()
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'MJ Model', ax = ax_M[0])
if(Matrix == 'BLOSUM62'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,25))
    M2_list = M2.tolist()
    Alphabet = np.array(['A', 'R'  ,'N' , 'D'  ,'C' , 'Q'  ,'E'  ,'G'  ,'H' , 'I'  ,'L'  ,'K'  ,'M' , 'F' , 'P' , 'S'  ,'T' , 'W' , 'Y' , 'V' , 'B' , 'Z'  ,'X',  '*'])
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'BLOSUM62 Model', ax = ax_M[0])
L_alphabet = len(Alphabet)

#-------------------------PWM----------------------------------
fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})

#antigen_seq = np.random.randint(0, len(Alphabet), L)
#antigen = Alphabet[antigen_seq]
antigen = 'FMLFMAVFVMTSWYC'
antigens = ['NTKTAATNLF', 'TACNSEYPNTTK', 'FMLFMAVFVMTSWYC', 'TACNSFVMTSATNLFSEYPN', 'TACNSFVMTSSEYPNPNEFYMAWYC']
colors_L = plt.cm.ocean(np.linspace(0,1,2*len(antigens)+4))

for l, antigen in enumerate(antigens):
    L=len(antigen)
    antigen_list = [i for i in antigen]
    contributions = np.zeros(shape = (1,20))
    antigen_seq = np.array([], dtype = int)
    for i, aa in enumerate(antigen_list):
        index = Alphabet_list.index(aa)
        antigen_seq = np.append(antigen_seq, int(index))
    print('Antigen:', antigen)
    antigen_seq, antigen
    PWM = M2[:,antigen_seq]
    for i in np.arange(L):
        PWM[:,i]-=np.min(PWM[:,i], axis=0)
    PWM_list = PWM.tolist()
    #plot_PWM(PWM=PWM_list, Alphabet=Alphabet, sequence = antigen, title=r'PWM', ax = ax_M[1])
    min_E = np.sum([np.min(PWM[:,i]) for i in range(len(PWM[0,:]))])
    avg_E = np.sum([np.mean(PWM[:,i]) for i in range(len(PWM[0,:]))])
    var_E = np.sum([np.var(PWM[:,i]) for i in range(len(PWM[0,:]))])
    max_E = np.sum([np.max(PWM[:,i]) for i in range(len(PWM[0,:]))])

    print(min_E, avg_E, max_E)
        
    ax_M[2].set_xticks(range(L))
    ax_M[2].set_xticklabels(antigen)
    ax_M[2].tick_params(labelsize = 22)

    linear_PWM = np.reshape(PWM, (L*20,1))
    energies = np.linspace(0, np.max(linear_PWM), 20)
    ds_eff = np.zeros_like(energies)
    for i, energy in enumerate(energies[:-1]):
        ds_eff[i] = np.size(np.where((linear_PWM<energies[i+1]) & (linear_PWM>energies[i])))/L

    ax_M[3].hist(linear_PWM, bins = 50, alpha = .5, align='left')
    my_plot_layout(ax=ax_M[3])

    fig_M.savefig('../../Figures/5_Geometric_exponent/Matrix_and_PWM.png')

    #-------------------------RELATIVE CLONE SIZE _ N----------------------------------
    # ------- from relative clone size -------

    alpha = 1
    gamma = 0
    beta = 0.5
    Ns = np.array([2e2, 2e3, 2e4, 2e5, 2e6])
    Ns = np.array([2e2, 2e3, 2e4, 2e5])
    markers = ['o', '^', '*', 's', 'X']
    exponents_rcs = np.array([])
    vars_rcs = np.array([])
    energy_model='MJ'
    linear = 0
    growth_models = ['exponential', 'linear']

    lambda_Ns_file = open(Text_files_path+'Dynamics/Ensemble/lambda_Ns_alpha-%.2f_beta-%.2f_L-%d_'%(alpha, beta, L)+growth_models[linear]+'_'+energy_model+'.pkl','rb')
    data =  pickle.load(lambda_Ns_file)
    lambdas_simulation = data[0]
    vars_rcs_exponent = data[1]
    lambda_Ns_file.close()


    #-------------------------lambda_N----------------------------------

    Tmin = .1238
    Tmax = 10

    #Ts = np.logspace(Log10Tmin, Log10Tmax, 2000)
    Ts = np.linspace(Tmin, Tmax, 10000)
    lambdas = 1/Ts[:-1]
    F_PWM = -Ts*np.log(Z_PWM(PWM, Ts))

    Us = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts))
    dU = np.diff(Us)
    s = np.cumsum((lambdas)[:-1]*dU)
    Ms = 1/np.cumsum(np.exp(s)/(2*np.cumsum(np.exp(s)*dU)[-1])*dU)
    Ms = (2*np.cumsum(np.exp(s)*dU)[-1])/np.cumsum(np.exp(s)*dU)

    es = np.linspace(min_E-20, max_E+20, 10000)
    de = es[1]-es[0]

    def P_e_gaussian(avg_E, var_E, es):
        return (2*np.pi*var_E)**(-0.5)*np.exp(-(es-avg_E)**2/(2*var_E))
    def P_min_e(N, es):
        return (N*(1-np.cumsum(P_e_gaussian(avg_E, var_E, es)*de))**(N-1)*(P_e_gaussian(avg_E, var_E, es)))

    #lambd_gaussian = np.array([])
    Ns_array = np.logspace(1, 8, 200)
    lambdas_N = np.array([])
    Us_N = np.array([])

    for N in Ns_array:
        #p_min_e = -np.diff((1-np.cumsum(P_e_gaussian(avg_E, var_E, es)*de))**N)/np.diff(es)
        #avg_min_E = np.sum(es[:-1]*p_min_e*de)
        #lambd_gaussian = np.append(lambd_gaussian, (avg_E-avg_min_E+2.95)/var_E)
        Us_N = np.append(Us_N, Us[np.where(Ms>N)[0][-1]])
        lambdas_N = np.append(lambdas_N, lambdas[np.where(Us<(Us_N[-1])+(2.5/15)*L)[0][-1]]) #+(2.5/15)*L


    #ax.plot(Ns_array, lambd_gaussian, label = 'Gaussian', linestyle = '--', linewidth= 4, color = 'darkolivegreen', alpha = .6)
    #ax.plot(Ns_array, lambdas_N, color = 'darkolivegreen', linewidth = 2, label = r'$\int_{\epsilon_{MS}}^{\epsilon_m}\Lambda(\epsilon)d\epsilon\sim\frac{1}{M} $', alpha = .6, linestyle = '--')
    ax.plot(Ns_array/(20**L), lambdas_N, color = colors_L[2*l+4], linewidth = 3, label = r'%d'%(L), alpha = .6, linestyle = '--')


    for n, N in enumerate(Ns):
        ax.scatter(N/(20**L), lambdas_simulation[n], color = colors_L[2*l+4], s = 50, marker = 's', alpha = 1)
        ax.errorbar(x = N/(20**L), y = lambdas_simulation[n], yerr = 1.96*2*(alpha*(lambdas_simulation[n])**2)/(beta)*np.sqrt(vars_rcs_exponent[n]), color = colors_L[2*l+4], linestyle = '', capsize = 4)

my_plot_layout(ax=ax, xscale = 'log', yscale = 'linear', ylabel = '$\lambda(N)$', xlabel = '$M/20^L$')
ax.legend(fontsize=24, title = r'$L$', title_fontsize = 26)
fig.savefig('../../Figures/5_Geometric_exponent/lambda_simulations_N_L.pdf')

end = time.time()
print('\nFinished in', '%.2f'%((end-start)/60), 'minutes')





