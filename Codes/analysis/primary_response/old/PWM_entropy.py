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

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

RT = 0.593 #kcal/mol
N_A = 6.02214076e23

fig_lambda, ax_lambda = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15})
fig_Q_P0, ax_Q_P0 = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15})
#ax_densities = ax_Q_P0.inset_axes([0.7, 0.7, 0.3, 0.3])
fig_Q_P0_R, ax_Q_P0_R = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15})
#ax_densities_R = ax_Q_P0_R.inset_axes([0.7, 0.7, 0.3, 0.3])
fig_DH, ax_DH = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15})

#Matrix = 'BLOSUM62'
Matrix = 'MJ2'
#Matrix = 'MM'

if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
if(Matrix == 'MM'):
    M2 = (np.loadtxt(Text_files_path + Matrix + '.txt')+1)*e0
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
if(Matrix == 'BLOSUM62'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,25))
    Alphabet = np.array(['A', 'R'  ,'N' , 'D'  ,'C' , 'Q'  ,'E'  ,'G'  ,'H' , 'I'  ,'L'  ,'K'  ,'M' , 'F' , 'P' , 'S'  ,'T' , 'W' , 'Y' , 'V' , 'B' , 'Z'  ,'X',  '*'])
L_alphabet = len(Alphabet)
#----------------------------------------------------------------------


file_names = np.array(['PWM_Adams_etal_2016_1.pkl', 'PWM_Adams_etal_2016_2.pkl', 'Collesano_2022.csv'], dtype = object)
file_names = np.array(['PWM_Adams_etal_2016_1.pkl'], dtype = object)
colors_data = np.array(['royalblue', 'indigo', 'indianred'])
labels_data = np.array(['H1', 'H3', 'MHC'])
#labels_data = np.array(['H1', 'MHC'])
file_exts = ['pkl', 'pkl', 'csv']
#file_exts = ['pkl', 'csv']
seqs = ['TFSDYWMNWV', 'GSYYGMDYWG', Alphabet[np.random.randint(0, len(Alphabet), 9)]]
#seqs = ['TFSDYWMNWV', Alphabet[np.random.randint(0, len(Alphabet), 9)]]

alphabet_data = np.array(['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'S', 'T', 'N', 'Q', 'C', 'P', 'H', 'K', 'R', 'D', 'E'])
WT_seq = np.array(['T', 'F', 'S', 'D', 'Y', 'W', 'M', 'N', 'W', 'V'])

L = len(WT_seq)
gamma = 1/60 #s^-1
#gamma = gamma*(60*60*24) # days^-1
k_on = 1e6 #(M*s)^-1
#k_on = k_on*(60*60*24) #(M*days)^-1
alpha = 1 # days^-1
alpha = alpha/(60*60*24)
e_MS = -23

#print(np.exp(min_E_data), np.exp(avg_E_data))
#----------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------------------------

Tmin = .01
Tmax = 50

es = np.linspace(e_MS, e_MS+50, 100)
de = es[1]-es[0]
#e0 = mean_lowest_effect
#d = d_eff
#L = 10
Ts = np.linspace(Tmin, Tmax, 20000)
lambdas = 1/Ts[:-1]

Hs_shannon = np.array([])
Hs_shannon_random = np.array([])

N_MJ = 1000
e_R_gauge  = 0

print('Simulating MJ matrices ...')

for n in tqdm(np.arange(N_MJ)):
    antigen_seq = np.random.randint(0, len(Alphabet), L)
    PWM = M2[:,antigen_seq]
    for i in np.arange(L):
        PWM[:,i]-=np.min(PWM[:,i], axis=0)

    min_E = np.sum([np.min(PWM[:,i]) for i in np.arange(len(PWM[0,:]))])
    avg_E = np.sum([np.mean(PWM[:,i]) for i in np.arange(len(PWM[0,:]))])
    var_E = np.sum([np.var(PWM[:,i]) for i in np.arange(len(PWM[0,:]))])
    max_E = np.sum([np.max(PWM[:,i]) for i in np.arange(len(PWM[0,:]))])
    F_PWM = -Ts*np.log(Z_PWM(PWM, Ts))
    Us = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts))+ e_MS
    dU = np.diff(Us)
    #e_R = Us[np.where(lambdas<1)][0] + e_R_gauge
    e_R = 8
    #----------------------------------------------------------------------
    Omega = 20**L
    Z = Z_PWM(PWM, 1)
    S = np.cumsum(lambdas[:-1]*dU)
    S_gaussian = -(Us[:-1]-avg_E)**2/(2*var_E) + L*np.log(20) - 0.5*np.log(2*np.pi*var_E)
    
    P0 = 1/Omega
    Q0 = np.exp(S)*P0
    Q0_gaussian = np.exp(S_gaussian)*P0

    P = np.exp(-Us[:-1])/Z
    Q = np.exp(S)*P
    Q_guassian = np.exp(S_gaussian)*P

    R = (1/(1+np.exp(Us[:-1]-e_R)))
    QR = np.exp(S)*R

    KL_D = -np.sum(Q*np.log(Q)*dU)
    #----------------------------------------------------------------------
    #ax.plot(Us, H_random, linestyle = ':')
    
    #ax_densities.plot(Us[:-1], H_random, color = 'grey', linestyle = '--', alpha = .1)

    #if(n%int(N_MJ/10)==0):
        #ax_Q_P0.plot(Us[:-1][::10], np.log(Q/P0)[::10], linewidth = 1, alpha = .2,  linestyle = '-', marker = '', ms = 8, color = 'grey')
        #ax_Q_P0.plot(Us[:-1][::10], (H+np.log(R)-H_random)[::10], linewidth = 1, alpha = .4,  linestyle = '--', marker = '', ms = 8,  color = 'grey')
    
    KL_D = np.sum(np.log(Q/P0)*Q*dU)
    #KL_D2 = np.sum((H+np.log(R)-H_random)*np.exp(H+np.log(R))*dU)
    #----------------------------------------------------------------------
    PWM_exp = np.exp(-PWM)
    for i in np.arange(L):
        PWM_exp[:, i]=PWM_exp[:, i]/np.sum(PWM_exp[:, i])
    H_shannon = 0
    for i in np.arange(L):
        H_shannon+=(np.sum(-PWM_exp[:,i]*np.log(PWM_exp[:,i])))
    Hs_shannon = np.append(Hs_shannon, H_shannon)
    H_shannon_random = L*(0.5+0.5*np.log(2*np.pi*var_E))
    H_shannon_random = L*np.log(20)
    Hs_shannon_random = np.append(Hs_shannon_random, H_shannon_random)
    if(N_MJ==1):
        print('Entropy of PWM:', H_shannon/np.log(2), 'bits')
        print('Entropy of random sequences:', H_shannon_random/np.log(2), 'bits')
        print('Entropy of uniform sequences:', L*np.log(20)/np.log(2), 'bits')
        print('KL_D: ', (KL_D)/np.log(2), 'bits')
        #print(KL_D2/np.log(2), 'bits')
    #----------------------------------------------------------------------
    if(n%int(N_MJ/10)==0):
        ax_lambda.plot(Us[::10], lambdas[::10], linewidth = 1, alpha = .2, linestyle = '-', color = 'grey')
        #ax_densities.plot(Us[:-1], S, color = 'grey', linestyle = '-', alpha = .1)
    #----------------------------------------------------------------------

hist_H = ax_DH.hist((Hs_shannon_random - Hs_shannon)/(Hs_shannon_random) , bins = 20, density=True, label = 'MJ', color = 'grey', alpha = .8)
ax_DH.vlines([0, np.sum(hist_H[1][1:]*hist_H[0]*np.diff(hist_H[1]))], 0, ax_DH.get_ylim()[1], colors = ['grey'], linestyles = ['-', '--'], linewidth = 4)
#----------------------------------------------------------------------
for k, file_name in enumerate(file_names[:]):
    L = len(seqs[k])
    print('Analizing ', labels_data[k], ' ...')
    if(file_exts[k]=='pkl'):
        #loading PWM from Adams et.al.
        file_PWM = open(Text_files_path + 'Data/'+ file_name, 'rb')

        log10Kd_WT = -8.91169118

        #PWM_data = (pickle.load(file_PWM) - log10Kd_WT)
        PWM_data = (pickle.load(file_PWM) - log10Kd_WT)*np.log(10)#*RT
        file_PWM.close()
    if(file_exts[k]=='csv'):
        PWM_data = pd.read_csv(Text_files_path + 'Data/' + file_name).to_numpy()
        PWM_data = PWM_data*np.log(10)

    #logKd_OPT = logKd_WT - ((PWM_H1_data[0, 2]) + (PWM_H1_data[15, 3]) + (PWM_H3_data[1, 1]) + (PWM_H3_data[9, 2]) + (PWM_H3_data[19, 6]) + (PWM_H3_data[4, 8]))

    #Change values by the minimum
    for i in np.arange(L):
        PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)
    
    #PWM_data_list = PWM_data.tolist()

    min_E_data = np.sum([np.min(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + e_MS
    avg_E_data = np.sum([np.mean(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + e_MS
    var_E_data = np.sum([np.var(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))])
    max_E_data = np.sum([np.max(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + e_MS

    ax_lambda.plot(es[:], -(es[:]-avg_E_data)/(var_E_data), alpha = .6, linewidth = 3, color = colors_data[k], linestyle = '--') #printi gaussian exponent
    #ax_lambda.plot(es[:]/Us[np.where(lambdas<1)][0], -(es[:]-avg_E_data)/(var_E_data), alpha = .6, linewidth = 3, color = colors_data[k], linestyle = '--') #printi gaussian exponent
    #---------------------------------------------------------------------- #Calculate Free energy from Z(PWM)
    Tmin = .01
    Tmax = 50
    Ts = np.linspace(Tmin, Tmax, 20000)
    lambdas = 1/Ts[:-1]
    F_PWM = -Ts*np.log(Z_PWM(PWM_data, Ts))
    Us = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts)) + e_MS
    print(Us[0])
    dU = np.diff(Us)
    ax_lambda.plot(Us[::1], lambdas[::1], linewidth = 3, linestyle = '-', marker = '', ms = 8, color = colors_data[k], label = labels_data[k])
    #ax_lambda.plot(Us[::10]/Us[np.where(lambdas<1)][0], lambdas[::10], linewidth = 3, linestyle = '-', marker = '', ms = 8, color = colors_data[k], label = labels_data[k])
    #---------------------------------------------------------------------- #Calculate entropy functions and distributions
    Omega = 20**L
    Z = Z_PWM(PWM_data, 1)
    S = np.cumsum(lambdas[:-1]*dU)
    S_gaussian = -(Us[:-1]-avg_E_data)**2/(2*var_E_data) + L*np.log(20) - 0.5*np.log(2*np.pi*var_E_data)
    
    N_r = 1e9
    P0 = 1/Omega
    Q0 = np.exp(S)*P0
    Q0_gaussian = np.exp(S_gaussian)*P0

    P_eq = np.exp(-(Us[:-1]-e_MS))/Z
    Q_eq = np.exp(S)*P_eq
    Q_guassian = np.exp(S_gaussian)*P_eq

    p_a = (gamma/(gamma + (k_on*np.exp(Us[:-1]))**2 ) )

    KL_D = -np.sum(Q_eq*np.log(Q_eq)*dU)
    #----------------------------------------------------------------------
    #ax_densities.plot(Us[:-1], S - np.log(Omega) + np.log(1e9), color = colors_data[k], linestyle = '-')

    ax_Q_P0.plot(np.log10(np.exp(1))*Us[:-1][:], np.log10(Q0)[:], linewidth = 2, alpha = .4, linestyle = ':', marker = '', ms = 6, color = colors_data[k], label = labels_data[k])
    ax_Q_P0.plot(np.log10(np.exp(1))*Us[:-1][:], np.log10(Q_eq)[:], linewidth = 2, alpha = 1, linestyle = '-', marker = '', ms = 6, color = colors_data[k])

    ax_Q_P0.vlines(np.log10(np.exp(1))*Us[np.where((np.log10(Q_eq))==np.max((np.log10(Q_eq))))][0], ax_Q_P0.get_ylim()[0], 9, color = colors_data[k], linestyle = '--')
    ax_Q_P0.vlines(np.log10(np.exp(1))*Us[np.where(lambdas<1)][0], ax_Q_P0.get_ylim()[0], 9, color = colors_data[k], linestyle = ':')
    #----------------------------------------------------------------------
    #ax_densities_R.plot(np.log10(np.exp(1))*Us[:-1], S - np.log(Omega) + np.log(1e9),  color = colors_data[k], linestyle = '-')
    ax_Q_P0_R.plot(np.log10(np.exp(1))*Us[:-1][:], np.log10(Q0)[:], linewidth = 2, alpha = .4, linestyle = ':', marker = '', ms = 6, color = colors_data[k], label = labels_data[k])
    for r, rho_A in enumerate((1e3*np.logspace(0, 8, 20)/N_A)):

        p_e = rho_A*k_on/alpha
        Q_R = Q0*p_e*p_a*N_r
        #ax_Q_P0_R.hlines(np.log10(p_e*p_a[0]), np.log10(np.exp(1))*Us[:-1][0], np.log10(np.exp(1))*Us[:-1][-1])
        ax_Q_P0_R.plot(np.log10(np.exp(1))*Us[:-1][:], np.log10(p_e*p_a), linewidth = 2, alpha = .05*r, linestyle = '--', marker = '', ms = 6, color = colors_data[k])
        ax_Q_P0_R.plot(np.log10(np.exp(1))*Us[:-1][:], np.log10(Q_R), linewidth = 2, alpha = .1*r, linestyle = '-', marker = '', ms = 6, color = colors_data[k])
        ax_Q_P0_R.vlines(np.log10(np.exp(1))*Us[np.where((np.log10(Q_R))==np.max((np.log10(Q_R))))][0], ax_Q_P0_R.get_ylim()[0], 3, color = colors_data[k], linestyle = '-')
    ax_Q_P0_R.vlines(np.log10(np.exp(1))*Us[np.where(lambdas<1)][0], ax_Q_P0.get_ylim()[0], 3, color = colors_data[k], linestyle = ':')
    #----------------------------------------------------------------------
    PWM_data_exp = np.exp(-PWM_data)
    for i in np.arange(L):
        PWM_data_exp[:, i]=PWM_data_exp[:, i]/np.sum(PWM_data_exp[:, i])
    H_shannon_data = 0
    H_shannon_data_random = 0
    for i in np.arange(L):
        H_shannon_data+=(np.sum(-PWM_data_exp[:,i]*np.log(PWM_data_exp[:,i])))
        var_i = np.var(PWM_data[:,i])
        H_shannon_data_random+=(0.5+0.5*np.log(2*np.pi*var_i))
    H_shannon_data_random = L*np.log(20)
    ax_DH.vlines([((H_shannon_data_random- H_shannon_data))/(H_shannon_data_random)], 0, ax_DH.get_ylim()[1], colors = colors_data[k], linestyles = ['-'], label = labels_data[k], linewidth = 4)
       
    print('Entropy of PWM:', H_shannon_data/np.log(2), 'bits')
    print('Entropy of random sequences:', H_shannon_data_random/np.log(2), 'bits ;' , (H_shannon_data_random - H_shannon_data)/np.log(2), 'bits')
    #print('Entropy of uniform sequences:', L*np.log2(20), 'bits ;',  L*np.log2(20) - H_shannon_data/np.log(2) , 'bits')
    print('KL_D: ', (KL_D)/np.log(2), 'bits')
    #print(KL_D2/np.log(2), 'bits')


#----------------------------------------------------------------------
ax_lambda.hlines([0, 1], e_MS, e_MS + 50, color = 'grey', linestyle = '--')
ax_lambda.legend(loc = 0, fontsize = 18, title = r'Model', title_fontsize = 20)
my_plot_layout(ax=ax_lambda, xscale = 'linear', yscale = 'linear', ylabel = '$\lambda\equiv 1/T$', xlabel = r'$\epsilon\equiv U$')
ax_lambda.set_ylim(-.8, 4);
#ax_lambda.set_xlim(0, 4);
#----------------------------------------------------------------------
#ax_Q_P0.plot(Us[:-1][::100], np.log(R)[::100], linewidth = 2, alpha = .6, linestyle = '-', marker = '', ms = 6, color = 'darkred')
#ax_Q_P0.hlines([0], 0, avg_E_data, color = 'grey', linestyle = '--')
my_plot_layout(ax=ax_Q_P0, xscale = 'linear', yscale = 'linear', ylabel = r'$\log_{10}{Q(\epsilon)}$', xlabel = r'$\log_{10}{K_d} [M]$')
ax_Q_P0.legend(loc = 1, fontsize = 14, title = r'Sequence', title_fontsize = 16)
ax_Q_P0.set_ylim(bottom = -6, top = 10)
ax_Q_P0.set_xlim(right = 0, left = -10.5);
#my_plot_layout(ax=ax_densities, xscale = 'linear', yscale = 'linear', ylabel = '$\log{\Lambda}$', xlabel = r'$\epsilon$', ticks_labelsize = 16)
#ax_densities.xaxis.tick_top()
#ax_densities.xaxis.set_label_position('top') 
#----------------------------------------------------------------------
#ax_Q_P0.plot(Us[:-1][::100], np.log(R)[::100], linewidth = 2, alpha = .6, linestyle = '-', marker = '', ms = 6, color = 'darkred')
ax_Q_P0_R.hlines([0], np.log10(np.exp(1))*e_MS, np.log10(np.exp(1))*(e_MS+50), color = 'grey', linestyle = '--')
my_plot_layout(ax=ax_Q_P0_R, xscale = 'linear', yscale = 'linear', ylabel = r'$\log_{10}{Q(\epsilon)}$', xlabel = r'$\log_{10}{K_d} [M]$')
ax_Q_P0_R.legend(loc = 2, fontsize = 14, title = r'Sequence', title_fontsize = 16)
ax_Q_P0_R.set_ylim(bottom = -12, top = 10)
ax_Q_P0_R.set_xlim(right = 0, left = -10.5);
#ax_densities_R.hlines(0, np.log10(np.exp(1))*e_MS, np.log10(np.exp(1))*(e_MS+50), color = 'grey', linestyle = '--' )
#my_plot_layout(ax=ax_densities_R, xscale = 'linear', yscale = 'linear', ylabel = '$\log{\Lambda}$', xlabel =r'$\log_{10}{K_d} [M]$', ticks_labelsize = 10, x_fontsize=14, y_fontsize=16)
#ax_densities_R.set_xlim(right=0)
#ax_densities_R.xaxis.tick_top()
#ax_densities_R.xaxis.set_label_position('top') 

#----------------------------------------------------------------------
my_plot_layout(ax=ax_DH, xscale = 'linear', yscale = 'linear', xlabel = '$\Delta S$')
ax_DH.legend(loc = 2, fontsize = 18, title = r'Model', title_fontsize = 20)

fig_lambda.savefig('../../Figures/6_Data/Data_lambda.pdf')
fig_Q_P0.savefig('../../Figures/6_Data/Data_2016_Q_P0.pdf')
fig_Q_P0_R.savefig('../../Figures/6_Data/Data_2016_Q_P0_R.pdf')
fig_DH.savefig('../../Figures/6_Data/Data_2016_DH.pdf')

