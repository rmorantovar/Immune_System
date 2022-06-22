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

fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig1, ax1 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig22, ax22 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig3, ax3 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

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

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'FTSENAYCGR'
antigen = 'TACNSEYPNTTK'
#antigen = 'TANSEYPNTK'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'

transparency_q = [1, .6, .3]
colors_q = ['darkred', 'darkred', 'darkred']
colors_q2 = ['darkred', 'indigo', 'darkred']
colors_R = ['tab:olive', 'olive', 'olive']
energy_models = ['MJ']
models_name = ['exponential', 'linear', ]
colors = ['tab:blue', 'tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = [0]#, 1]

L=len(antigen)

N_r = 2e5
T0 = 3
Tf = 8
#Tf = 9
dT = 0.1
days = np.linspace(3, Tf, 3)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_act = .1 # hour^-1
k_act = k_act*24 #days^-1
qs = [1, 2]
beta = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e3
E_ms = -27

print('max_u = %.2e'%(k_on*np.exp(Tf*lambda_A)/N_A))

print('k_on/k_act = %.1e'%(k_on/k_act))


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
Ks = np.exp(Es[:-1])
beta_r = lambdas[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]

ax2.plot(Ks, Q0*N_r, alpha = transparency_q[0], color = 'grey', linewidth = 5, linestyle = '-')
ax22.plot(Ks, Q0*N_r, alpha = transparency_q[0], color = 'grey', linewidth = 5, linestyle = '-')    
for n_q, q in enumerate(qs):
    
    #--------------------------p_a(E, t)---------------------------
    #ax0.vlines(k_act/k_on, ax0.get_ylim()[0], 1, color = 'grey', linestyle = ':')
    for n_t, t in enumerate(days[[-3, -2, -1]]):
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_act, np.exp(lambda_A*t)/N_A, Es, q, lambda_A, N_c, dE)
        M_r = N_r*N_c*np.sum(Q0*p_a*dE)
        #--------------------------R(E, t)---------------------------
        if q==2:
            ax1.plot(Ks, R, alpha = transparency_q[0], color = colors_R[n_t], linewidth = 5, linestyle = '-')
        #ax1.hlines(r_a[0]*N_c/lambda_A, ax1.get_xlim()[0], ax1.get_xlim()[1], alpha = transparency_q[n_q], color = colors_q[n_q], linestyle = ':' )
        ax1.set_ylim(bottom = 1e-8, top = 2)
        #----------------------------------------------------------------
        #--------------------------QR(E, t)---------------------------
        #if q==1:
            #ax2.plot(Ks, Q0*N_r, alpha = transparency_q[1], color = 'grey', linewidth = 5, linestyle = '-')
            #ax2.plot(Ks, QR*N_r, alpha = transparency_q[1], color = colors_R[n_t], linewidth = 5, linestyle = '-')
            #ax2.vlines(np.exp(Es)[QR == np.max(QR)][0], ax2.get_ylim()[0], np.max(QR*N_r), color = colors_R[n_t], linestyle = ':')
            #ax2.hlines([1, N_r], ax2.get_xlim()[0], ax2.get_xlim()[1], alpha = 1, color = 'black', linestyle = ':')
            #ax2.vlines(np.exp(Es)[lambdas[:-1] < q][0], ax2.get_ylim()[0], np.max(QR), color = colors_R[-1], linestyle = ':', linewidth = 2)
        if q==2:
            ax2.plot(Ks, QR*N_r, alpha = transparency_q[0], color = colors_R[n_t], linewidth = 5, linestyle = '-')
            #ax2.hlines([N_r], ax2.get_xlim()[0], ax2.get_xlim()[1], alpha = 1, color = 'black', linestyle = ':') 
            #-------FOR Q0--------- 
            #ax2.vlines(np.exp(E_r), ax2.get_ylim()[0], N_r*Q0[Ks<np.exp(E_r)][-1], color = 'black', linestyle = ':')       
            #ax2.hlines(N_r*Q0[Ks<np.exp(E_r)][-1], ax2.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')
            #---------------------- 

        #ax2.vlines(Ks[lambdas[:-1] < 1][0], ax2.get_ylim()[0], np.max(QR), color = 'brown', linestyle = ':', linewidth = 4)
        ax2.set_ylim(bottom = 1e-8, top = 2*N_r)
    ax0.plot(Ks, p_a, color = colors_q[n_q], alpha = transparency_q[n_q], linewidth = 5, linestyle = '-', label = '%d'%(q))
    ax22.plot(Ks, QR*N_r, alpha = transparency_q[0], color = colors_q2[n_q], linewidth = 5, linestyle = '-')
    ax22.hlines([1, N_r], ax22.get_xlim()[0], ax22.get_xlim()[1], alpha = 1, color = 'black', linestyle = ':')
    #ax22.vlines(Ks[QR == np.max(QR)][0], ax22.get_ylim()[0], np.max(QR*N_r), color = colors_q2[n_q], linestyle = ':')
    #ax22.vlines(Ks[lambdas[:-1] < q][0], ax22.get_ylim()[0], np.max(QR), color = colors_q2[n_q], linestyle = '--', linewidth = 2)
        #----------------------------------------------------------------
    #--------------------------m(t)---------------------------
    m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])
    m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))
    ax3.plot(time, m_bar, linewidth = 4, linestyle = '-', color = colors_q[n_q])
    ax3.plot(time, m_bar_approx, linewidth = 3, linestyle = '--', color = 'black')
    ax3.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
    #----------------------------------------------------------------

my_plot_layout(ax=ax0, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
ax0.legend(title = '$q$', title_fontsize = 35, fontsize = 30)
#ax0.legend(fontsize = 30, title_fontsize=33)
fig0.savefig('../../Figures/Summary/p_a.pdf')

my_plot_layout(ax=ax1, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
fig1.savefig('../../Figures/Summary/R_clone.pdf')

my_plot_layout(ax=ax2, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
fig2.savefig('../../Figures/Summary/QR.pdf')
my_plot_layout(ax=ax22, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
fig22.savefig('../../Figures/Summary/Q2.pdf')

my_plot_layout(ax=ax3, yscale = 'log', ticks_labelsize = 30)
fig3.savefig('../../Figures/Summary/activation_rate.pdf')










