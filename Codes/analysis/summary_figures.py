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
from matplotlib.collections import LineCollection
from scipy.optimize import curve_fit

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

fig_beta, ax_beta = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

N_A = 6.02214076e23
k_BT = 1.380649e-23*293
style.use('seaborn-paper')

M1 = np.loadtxt(Text_files_path+'MJ.txt', skiprows= 1, usecols=range(1,21)).tolist()
M2 = (np.loadtxt(Text_files_path+'MJ2.txt', skiprows= 1, usecols=range(1,21))).tolist()
M3 = np.loadtxt(Text_files_path+'BLOSUM62.txt', skiprows= 1, max_rows = 23, usecols=range(1,24)).tolist()
Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'P_act', 's', 't']
Alphabet2 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'P_act', 's', 't', 'u', 'v', 'w']
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
antigen = 'TACNSEYPNTTKCGRWYC'
#antigen = 'TANSEYPNTK'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'

transparency_q = [1, 1, .3, 0]
colors_theta = ['darkred', 'olive', 'darkblue']
colors_theta2 = ['tab:red', 'tab:olive', 'tab:blue']

colors_R = [['tab:red', 'tab:red', 'tab:red', 'darkred'], ['tab:olive', 'tab:olive', 'tab:olive', 'olive'], ['tab:blue', 'tab:blue', 'tab:blue', 'darkblue']]

energy_models = ['MJ']
models_name = ['exponential', 'linear', ]
growth_models = [0] #, 1]

L=len(antigen)
print('L=%.d'%L)

N_r = 1e8
T0 = 3
Tf = 6
#Tf = 9
dT = 0.1
days = np.linspace(2, Tf, 5)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_pr = .1 # hour^-1
k_pr = k_pr*24 #days^-1
thetas = [1, 1.5, 2]
beta = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
E_ms = -27.63

print('K_d_ms=%.1e'%np.exp(E_ms))

print('max_u = %.2e'%(k_on*np.exp(Tf*lambda_A)/N_A))

print('k_pr/k_on = %.1e'%(k_on/k_pr)**(-1))


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
Ks = np.exp(Es[:-1])
beta_r = lambdas[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]

#----------------------------------------------------------------
points = np.array([Ks, lambdas[:-1]]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
norm = plt.Normalize(0, 5)
lc = LineCollection(segments, cmap='turbo_r', norm=norm)
# Set the values used for colormapping
lc.set_array(lambdas)
lc.set_linewidth(5)
line = ax_beta.add_collection(lc)
# fig_beta.colorbar(line, ax=ax_beta)

my_plot_layout(ax=ax_beta, yscale = 'linear', xscale = 'log', ticks_labelsize = 38)
ax_beta.set_xticks([])
ax_beta.set_yticks(np.arange(0, 6))
ax_beta.set_ylim(top = 5, bottom = -.2)
fig_beta.savefig('../../Figures/_Summary/beta.pdf')
#----------------------------------------------------------------


for n_thetas, theta in enumerate(thetas):
    E_theta = Es[lambdas>theta][-1]
    fig_P_act, ax_P_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_Q_act, ax_Q_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_Q_act2, ax_Q_act2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_m, ax_m = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_H, ax_H = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

    ax_Q_act.plot(Ks, Q0*N_r, alpha = transparency_q[0], color = 'grey', linewidth = 5, linestyle = '-')
    ax_Q_act2.plot(Ks, Q0*N_r, alpha = transparency_q[0], color = 'grey', linewidth = 5, linestyle = '-')    

    print('theta = %.1f'%theta)
    #--------------------------m(t)---------------------------
    u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, theta, lambda_A, N_c, dE)
    M_r = N_r*N_c*np.sum(Q0*p_a*dE)
    m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])
    m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))

    ax_m.plot(time, m_bar, linewidth = 4, linestyle = '-', color = colors_theta[n_thetas])
    ax_m.plot(time, m_bar_approx, linewidth = 3, linestyle = '--', color = 'black')
    ax_m.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
    #----------------------------------------------------------------

    t_act = time[m_bar>1][0]

    print(t_act)
    
    days = np.linspace(2, t_act, 4)
    for n_t, t in enumerate(days[[-4, -3, -2, -1]]):
        u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t)/N_A, Es, theta, lambda_A, N_c, dE)
        #ax_P_act.hlines(r_a[0]*N_c/lambda_A, ax_P_act.get_xlim()[0], ax_P_act.get_xlim()[1], alpha = transparency_q[n_thetas], color = colors_theta[n_thetas], linestyle = ':' )
        ax_P_act.set_ylim(bottom = 1e-11, top = 2)
        #----------------------------------------------------------------
        #--------------------------Q_act(E, t)---------------------------
        if theta!=0:
            #--------------------------P_act(E, t)---------------------------
            ax_P_act.plot(Ks, P_act, alpha = transparency_q[0], color = colors_R[n_thetas][n_t], linewidth = 5, linestyle = '-')
            ax_Q_act.plot(Ks, Q_act*N_r, alpha = transparency_q[0], color = colors_R[n_thetas][n_t], linewidth = 5, linestyle = '-')
            #-------FOR Q0--------- 
            ax_Q_act.vlines(np.exp(E_r), ax_Q_act.get_ylim()[0], N_r*Q0[Ks<np.exp(E_r)][-1], color = 'black', linestyle = ':')
            ax_Q_act.vlines(np.exp(E_theta), ax_Q_act.get_ylim()[0], N_r*Q0[Ks<np.exp(E_theta)][-1], color = 'black', linestyle = ':')       
            #ax_Q_act.hlines(N_r*Q0[Ks<np.exp(E_r)][-1], ax_Q_act.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')
            #---------------------- 
        # if theta==2:
        #     #--------------------------P_act(E, t)---------------------------
        #     ax_P_act.plot(Ks, P_act, alpha = transparency_q[0], color = colors_R2[n_t], linewidth = 5, linestyle = '-')

        #     ax_Q_act.plot(Ks, Q_act*N_r, alpha = transparency_q[0], color = colors_R2[n_t], linewidth = 5, linestyle = '-')
        #     #-------FOR Q0--------- 
        #     ax_Q_act.vlines(np.exp(E_r), ax_Q_act.get_ylim()[0], N_r*Q0[Ks<np.exp(E_r)][-1], color = 'black', linestyle = ':')
        #     ax_Q_act.vlines(np.exp(E_theta), ax_Q_act.get_ylim()[0], N_r*Q0[Ks<np.exp(E_theta)][-1], color = 'black', linestyle = ':')      
        #     #ax_Q_act.hlines(N_r*Q0[Ks<np.exp(E_r)][-1], ax_Q_act.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')
        #     #---------------------- 


    for n_t, t in enumerate(time[::4]):
        u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t)/N_A, Es, theta, lambda_A, N_c, dE)
        # ------------ H ---------------
        dE_temp = dE[Q_act!=0]
        Q0_temp = Q0[Q_act!=0]
        Ks_temp = Ks[Q_act!=0]
        Q_act_temp = Q_act[Q_act!=0]
        D_KL_t = np.sum(dE_temp*(Q_act_temp/np.sum(Q_act_temp*dE_temp))*(np.log((Q_act_temp/np.sum(Q_act_temp*dE_temp)))-np.log(Q0_temp)))
        ax_H.plot(t, D_KL_t, marker = 's', ms = 12, color = colors_theta2[n_thetas])
        
    u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t_act)/N_A, Es, theta, lambda_A, N_c, dE)
    dE_temp = dE[Q_act!=0]
    Q0_temp = Q0[Q_act!=0]
    Ks_temp = Ks[Q_act!=0]
    Q_act_temp = Q_act[Q_act!=0]
    D_KL_t = np.sum(dE_temp*(Q_act_temp/np.sum(Q_act_temp*dE_temp))*(np.log((Q_act_temp/np.sum(Q_act_temp*dE_temp)))-np.log(Q0_temp)))
    ax_H.plot(t_act, D_KL_t, marker = 's', ms = 16, color = colors_theta[n_thetas])
    # -----------------------------

    ax_Q_act.set_ylim(bottom = 1e-11, top = 2*N_r)
    
    ax0.plot(Ks, p_a, color = colors_theta[n_thetas], alpha = transparency_q[n_thetas], linewidth = 5, linestyle = '-', label = '%.1f'%(theta))

    ax_Q_act2.plot(Ks, Q_act*N_r, alpha = transparency_q[0], color = colors_theta[n_thetas], linewidth = 5, linestyle = '-')
    u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*days[-1])/N_A, Es, theta, lambda_A, N_c, dE)
    m = np.sum(dE*Q_act*N_r)

    print('# of activated lineages : %d'%m , m_bar[-1])

    ax_Q_act2.hlines([1, N_r], ax_Q_act2.get_xlim()[0], ax_Q_act2.get_xlim()[1], alpha = 1, color = 'black', linestyle = ':')
    #ax_Q_act2.vlines(Ks[Q_act == np.max(Q_act)][0], ax_Q_act2.get_ylim()[0], np.max(Q_act*N_r), color = colors_theta2[n_thetas], linestyle = ':')
    #ax_Q_act2.vlines(Ks[lambdas[:-1] < theta][0], ax_Q_act2.get_ylim()[0], np.max(Q_act), color = colors_theta2[n_thetas], linestyle = '--', linewidth = 2)
        #----------------------------------------------------------------

    my_plot_layout(ax=ax_P_act, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
    ax_P_act.set_xticks([])
    fig_P_act.savefig('../../Figures/_Summary/R_clone.pdf')

    my_plot_layout(ax=ax_Q_act, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
    #ax_Q_act.set_xticks([])
    fig_Q_act.savefig('../../Figures/_Summary/QR_q-%.1f.pdf'%theta)
    my_plot_layout(ax=ax_Q_act2, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
    fig_Q_act2.savefig('../../Figures/_Summary/QR2_q-%.1f.pdf'%theta)

    my_plot_layout(ax=ax_m, yscale = 'log', ticks_labelsize = 30)
    fig_m.savefig('../../Figures/_Summary/activation_rate_q-%.1f.pdf'%theta)

    my_plot_layout(ax=ax_H, yscale = 'linear', ticks_labelsize = 30)
    ax_H.set_ylim(bottom = 0, top = 17)
    fig_H.savefig('../../Figures/_Summary/entropy_q-%.1f.pdf'%theta)

my_plot_layout(ax=ax0, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
ax0.legend(title = '$theta$', title_fontsize = 35, fontsize = 30)
#ax0.legend(fontsize = 30, title_fontsize=33)
fig0.savefig('../../Figures/_Summary/p_a.pdf')










