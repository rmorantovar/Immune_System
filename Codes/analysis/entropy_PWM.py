
import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
from Immuno_models import*
import scipy.special as sc
import pickle
from matplotlib import style
from scipy.optimize import curve_fit
from tqdm import tqdm

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

fig_lambda, ax_lambda = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15})
fig_Q_P0, ax_Q_P0 = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15})
ax_densities = ax_Q_P0.inset_axes([0.65, 0.3, 0.35, 0.35])
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


file_names = np.array([Text_files_path + '/Data/PWM_Adams_etal_2016_1.pkl', Text_files_path + '/Data/PWM_Adams_etal_2016_2.pkl'], dtype = object)
colors_data = np.array(['royalblue', 'indigo'])
labels_data = np.array(['H1', 'H3'])

alphabet_data = np.array(['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'S', 'T', 'N', 'Q', 'C', 'P', 'H', 'K', 'R', 'D', 'E'])
WT_seq = np.array(['T', 'F', 'S', 'D', 'Y', 'W', 'M', 'N', 'W', 'V'])
L = len(WT_seq)


#print(np.exp(min_E_data), np.exp(avg_E_data))
#----------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------------------------

Tmin = .00001
Tmax = 10

es = np.linspace(0, 30, 100)
de = es[1]-es[0]
#e0 = mean_lowest_effect
#d = d_eff
#L = 10
Ts = np.linspace(Tmin, Tmax, 200000)
lambdas = 1/Ts[:-1]

Hs_shannon = np.array([])
Hs_shannon_random = np.array([])

N_MJ = 10
e_R_gauge  = 0

print('Simulating MJ matrices ...')

for n in tqdm(np.arange(N_MJ)):
    antigen_seq = np.random.randint(0, len(Alphabet), L)
    PWM = M2[:,antigen_seq]
    for i in np.arange(L):
        PWM[:,i]-=np.min(PWM[:,i], axis=0)

    Us_bar = np.sum([np.mean(PWM[:,i]) for i in np.arange(L)])
    Us_var = np.sum([np.var(PWM[:,i]) for i in np.arange(L)])
    F_PWM = -Ts*np.log(Z_PWM(PWM, Ts))
    Us = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts))
    e_R = Us[np.where(lambdas<1)][0] + e_R_gauge
    #e_R = 5
    if(n%int(N_MJ/10)==0):
        ax_lambda.plot(Us[::10], lambdas[::10], linewidth = 1, alpha = .2, linestyle = '-', color = 'grey')
    #----------------------------------------------------------------------
    dU = np.diff(Us)
    H = np.cumsum(lambdas[:-1]*dU)
    R = (1/(1+np.exp(Us[:-1]-e_R)))
    #Q = np.exp(H)/np.sum(np.exp(H)*dU)
    Q = (np.exp(H)*R)/np.sum(np.exp(H)*R*dU)
    #ax.plot(Us[:-1], H, linestyle = '--')
    H_random = -(Us[:-1]-Us_bar)**2/(2*Us_var) + L*np.log(20) - 0.5*np.log(2*np.pi*Us_var)
    P0 = np.exp(H_random)/np.sum(np.exp(H_random)*dU)
    #ax.plot(Us, H_random, linestyle = ':')
    ax_densities.plot(Us[:-1], H + np.log(R), color = 'grey', linestyle = '-', alpha = .1)
    ax_densities.plot(Us[:-1], H_random, color = 'grey', linestyle = '--', alpha = .1)

    if(n%int(N_MJ/10)==0):
        ax_Q_P0.plot(Us[:-1][::10], np.log(Q/P0)[::10], linewidth = 1, alpha = .2,  linestyle = '-', marker = '', ms = 8, color = 'grey')
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
    H_shannon_random = L*(0.5+0.5*np.log(2*np.pi*Us_var))
    Hs_shannon_random = np.append(Hs_shannon_random, H_shannon_random)
    if(N_MJ==1):
        print('Entropy of PWM:', H_shannon/np.log(2), 'bits')
        print('Entropy of random sequences:', H_shannon_random/np.log(2), 'bits')
        print('Entropy of uniform sequences:', L*np.log2(20), 'bits')
        print('KL_D: ', (KL_D)/np.log(2), 'bits')
        #print(KL_D2/np.log(2), 'bits')
    #----------------------------------------------------------------------
hist_H = ax_DH.hist(Hs_shannon_random/np.log(2) - Hs_shannon/np.log(2) , bins = 20, density=True, label = 'MJ', color = 'grey', alpha = .8)
ax_DH.vlines([0, np.sum(hist_H[1][1:]*hist_H[0]*np.diff(hist_H[1]))], 0, ax_DH.get_ylim()[1], colors = ['grey'], linestyles = ['-', '--'], linewidth = 4)
#----------------------------------------------------------------------
for k, file_name in enumerate(file_names):
    #loading PWM from Adams et.al.
    file_PWM = open(file_name, 'rb')

    log10Kd_WT = -8.91169118

    PWM_data = (pickle.load(file_PWM) - log10Kd_WT)
    #PWM_data = np.log(10**(pickle.load(file_PWM) - log10Kd_WT))
    file_PWM.close()

    #logKd_OPT = logKd_WT - ((PWM_H1_data[0, 2]) + (PWM_H1_data[15, 3]) + (PWM_H3_data[1, 1]) + (PWM_H3_data[9, 2]) + (PWM_H3_data[19, 6]) + (PWM_H3_data[4, 8]))

    # #PWMs from data with the same order of aas as my MJ matrix
    # PWM_data_2 = np.ones_like(PWM_H1_data)

    # for i in np.arange(L_alphabet):
    #     PWM_data_2[i,:] = PWM_H1_data[np.where(Alphabet==alphabet_data[i])[0][0], :]

    #Change values by the minimum
    for i in np.arange(L):
        PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)
        #PWM_data_2[:,i]-=np.min(PWM_H1_data_2[:,i], axis=0)

    #PWM_data_list = PWM_data.tolist()
    #PWM_data_2_list = PWM_H1_data_2.tolist()


    min_E_data = np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    avg_E_data = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    var_E_data = np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    max_E_data = np.sum([np.max(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])


    print('Analizing ', labels_data[k], ' ...')
    Us_bar = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    Us_var = np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    ax_lambda.plot(es[:], -(es[:]-Us_bar)/(Us_var), alpha = .6, linewidth = 3, color = colors_data[k], linestyle = '--')
    F_PWM = -Ts*np.log(Z_PWM(PWM_data, Ts))
    Us = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts))
    e_R = Us[np.where(lambdas<1)][0] + e_R_gauge
    #e_R = 3
    ax_lambda.plot(Us[::10], lambdas[::10], linewidth = 3, linestyle = '-', marker = '', ms = 8, color = colors_data[k], label = labels_data[k])
    #----------------------------------------------------------------------
    dU = np.diff(Us)
    H = np.cumsum(lambdas[:-1]*dU)
    R = (1/(1+np.exp(Us[:-1]-e_R)))
    #Q = np.exp(H)/np.sum(np.exp(H)*dU)
    Q = (np.exp(H)*R)/np.sum(np.exp(H)*R*dU)
    #ax.plot(Us[:-1], H, linestyle = '--')
    H_random = -(Us[:-1]-Us_bar)**2/(2*Us_var) + L*np.log(20) - 0.5*np.log(2*np.pi*Us_var)
    P0 = np.exp(H_random)/np.sum(np.exp(H_random)*dU)
    #ax.plot(Us, H_random, linestyle = ':')

    ax_densities.plot(Us[:-1], H + np.log(R), color = colors_data[k], linestyle = '-')
    ax_densities.plot(Us[:-1], H_random, color = colors_data[k], linestyle = '--')
    ax_Q_P0.plot(Us[:-1][::100], np.log(Q/P0)[::100], linewidth = 2, alpha = .8, linestyle = '-', marker = '', ms = 6, color = colors_data[k], label = labels_data[k])
    #ax_Q_P0.plot(Us[:-1][::100], (H+np.log(R)-H_random)[::100], linewidth = 2, alpha = .8, linestyle = '--', marker = '', ms = 6, color = colors_data[k])
    
    ax_Q_P0.vlines(Us[np.where(lambdas<1)][0], -4, 12, color = colors_data[k], linestyle = '-')
    ax_Q_P0.vlines(Us[np.where(np.log(Q/P0)==np.max(np.log(Q/P0)))][0], -4, 12, color = colors_data[k], linestyle = ':')
    #ax_Q_P0.vlines(Us[np.where((H-H_random)==np.max((H-H_random)))][0], -4, 12, color = colors_data[k], linestyle = ':')

    KL_D = np.sum(np.log(Q/P0)*Q*dU)
    #KL_D2 = np.sum((H+np.log(R)-H_random)*np.exp(H+np.log(R))*dU)

    #----------------------------------------------------------------------
    PWM_data_exp = np.exp(-PWM_data)
    for i in np.arange(L):
        PWM_data_exp[:, i]=PWM_data_exp[:, i]/np.sum(PWM_data_exp[:, i])
    H_shannon_data = 0
    for i in np.arange(L):
        H_shannon_data+=(np.sum(-PWM_data_exp[:,i]*np.log(PWM_data_exp[:,i])))
    H_shannon_data_random = L*(0.5+0.5*np.log(2*np.pi*Us_var))
    #ax_DH.vlines([L*np.log2(20) - H_shannon/np.log(2)], 0, 1.1, colors = ['black'], linestyles = ['-'], label = r'Adams etal 2016 ($Z$)')
    ax_DH.vlines([H_shannon_data_random/np.log(2) - H_shannon_data/np.log(2)], 0, ax_DH.get_ylim()[1], colors = colors_data[k], linestyles = ['-'], label = labels_data[k], linewidth = 4)
    ax_DH.vlines([L*np.log2(20) - H_shannon_data/np.log(2)], 0, ax_DH.get_ylim()[1], colors = colors_data[k], linestyles = ['--'], linewidth = 4)
    ax_DH.vlines([KL_D/np.log(2)], 0, ax_DH.get_ylim()[1], colors = colors_data[k], linestyle= ':', linewidth = 4)
    #ax_DH.vlines([KL_D2/np.log(2)], 0, 1.1, colors = 'red', linestyle= ':', label = r'Adams etal 2016 $(D_{KL}2)$')
    print('Entropy of PWM:', H_shannon_data/np.log(2), 'bits')
    print('Entropy of random sequences:', H_shannon_data_random/np.log(2), 'bits ;' , H_shannon_data_random/np.log(2) - H_shannon_data/np.log(2), 'bits')
    print('Entropy of uniform sequences:', L*np.log2(20), 'bits ;',  L*np.log2(20) - H_shannon_data/np.log(2) , 'bits')
    print('KL_D: ', (KL_D)/np.log(2), 'bits')
    #print(KL_D2/np.log(2), 'bits')

#----------------------------------------------------------------------
ax_lambda.hlines([0, 1], 0, Us_bar, color = 'grey', linestyle = '--')
ax_lambda.legend(loc = 0, fontsize = 18, title = r'Model', title_fontsize = 20)
my_plot_layout(ax=ax_lambda, xscale = 'linear', yscale = 'linear', ylabel = '$\lambda\equiv 1/T$', xlabel = r'$\epsilon\equiv U$')
ax_lambda.set_ylim(-.8, 4);
#----------------------------------------------------------------------
ax_Q_P0.hlines([0], 0, Us_bar, color = 'grey', linestyle = '--')
my_plot_layout(ax=ax_Q_P0, xscale = 'linear', yscale = 'linear', ylabel = '$\log{P/P_0}$', xlabel = r'$\epsilon$')
ax_Q_P0.legend(loc = 1, fontsize = 14, title = r'CDR3 segment (Adams et.al.)', title_fontsize = 16)
ax_Q_P0.set_ylim(-4.5, 12.5);
my_plot_layout(ax=ax_densities, xscale = 'linear', yscale = 'linear', ylabel = '$\Lambda$', xlabel = r'$\epsilon$', ticks_labelsize = 16)
ax_densities.xaxis.tick_top()
ax_densities.xaxis.set_label_position('top') 

#----------------------------------------------------------------------
my_plot_layout(ax=ax_DH, xscale = 'linear', yscale = 'linear', xlabel = '$\Delta H_s$')
ax_DH.legend(loc = 2, fontsize = 18, title = r'Model', title_fontsize = 20)

fig_lambda.savefig('../../Figures/6_Data/Data_Adams_2016_lambda.pdf')
fig_Q_P0.savefig('../../Figures/6_Data/Data_Adams_2016_Q_P0.pdf')
fig_DH.savefig('../../Figures/6_Data/Data_Adams_2016_DH.pdf')

