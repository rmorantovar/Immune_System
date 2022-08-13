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
thetas = [2.2, 2.0, 1.8, 1.5, 1]
#thetas = [2.2, 2.0, 1.8]

colors_theta = ['lightblue','darkblue', 'olive', 'orange', 'darkred']
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
E_ms = -27.63
C = 5e3

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'FTSENAYCGR'
antigen = 'TACNSEYPNTTK'
antigen = 'TACNSEYPNTTKCGRWYC'


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

Es, dE, Q0, betas = calculate_Q0(0.01, 50, 200000, PWM_data, E_ms, L)
S = np.log(Q0)
Kds = np.exp(Es[:-1])

beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
print('beta_r = %.1f'%beta_r)
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
Kd_r = np.exp(E_r)

E_pr = Es[:-1][Kds<(k_pr/k_on)][-1]
Kd_pr = np.exp(E_pr)
t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
#----------------------------------------------------------------

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
colors = ['tab:blue', 'tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = [0]
linear = 0

fig_total_pop, ax_total_pop = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
fig_H, ax_H = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
fig_time, ax_time = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
fig_clone_size, ax_clone_size = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
fig_m, ax_m = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
fig_NC, ax_NC = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})

for i_theta, theta in enumerate(thetas):

    E_t = lambda t, theta:lambda_A*t/theta - np.log((lambda_A*N_A)/(k_on*N_c))/theta + np.log(k_pr/k_on) 
    
    beta_theta = betas[betas>theta][-1]
    E_theta = Es[betas>theta][-1]
    Kd_theta = np.exp(E_theta)

    delta_t_theta = (E_theta-E_pr)*theta/lambda_A
    t_theta = t_prime + delta_t_theta

    time1 = np.linspace(0, t_theta, 100)
    time2 = np.linspace(t_theta, Tf, 100)

    ax_H.plot(time1, -1*np.ones_like(time1)*S[Es[:-1]<E_theta][-1], color = colors_theta[i_theta])
    ax_H.plot(time2, [-S[Es[:-1]<E][-1] for E in E_t(time2, theta)] , label = '%d'%theta, color = colors_theta[i_theta])
        
    fig_a, ax_a = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
    #fig_act, ax_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})
    fig_clones, ax_clones = plt.subplots(1, 3, figsize=(30,10), gridspec_kw={'left':0.06, 'right':.98, 'bottom':.1, 'top': 0.9})
    fig_clones2, ax_clones2 = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

    colors_activation = []

    parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, theta, linear, N_ens)+energy_model
    data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies.txt', sep = '\t', header=None)

    min_e_data = np.min(data[0])
    max_e_data = np.max(data[0])

    data_active = data.loc[data[1]==1]

    activations_times = np.array(data_active[3])
    energies  = np.array(data_active[0])

    ar1, ar2 = np.histogram(activations_times, bins = time)

    data_N_active_linages = np.cumsum(ar1)

    print('Activated clones:',data_N_active_linages[-1], np.shape(data_active))
    
    t0 = np.where(data_N_active_linages==0)[0][-1]*dT

    for i in data[1]:
        if(i==0):
            colors_activation.append('silver')
        else:
            colors_activation.append(colors[0])

    #---- Antigen ----

    #---- B cell linages ----
    expfit = 1e-5*np.exp(lambda_B*time)
    clone_sizes = np.ones((int(data_N_active_linages[-1]), len(time)))
    for i in np.arange(0, int(data_N_active_linages[-1])):
        clone_sizes[i, int(activations_times[i]/dT):] = np.exp(lambda_B*(time[int(activations_times[i]/dT):] - activations_times[i] ))

    #---- Activation rate ----
    #--------------------------m(t)---------------------------
    u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, theta, lambda_A, N_c, dE)
    M_r = N_r*N_c*np.sum(Q0*p_a*dE)
    m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])
    m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))
    t_act = time[m_bar>1][0]
    #---------------------------------------------------------------- 
    #----------------------------------------------------------------
    u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf)/N_A, Es, theta, lambda_A, N_c, dE)
    psi = ((k_on*N_c)/(N_A))*p_a
    M_r = N_r*N_c*np.sum(Q0*p_a*dE)
    print(M_r)
    #----------------------------------------------------------------
    ax_a.plot(time[:-1], data_N_active_linages[:], linestyle = '', marker = 'o', ms = 5, linewidth = 2, label = 'simulation', color = colors_theta[i_theta], alpha = .6)
    ax_m.plot(time[:-1], data_N_active_linages[:], linestyle = '', marker = 'o', ms = 5, linewidth = 2, label = 'simulation', color = colors_theta[i_theta], alpha = .6)
    if(linear == 0):
        ax_a.plot(time, np.ones_like(time), color = 'black', linestyle = '-', linewidth = 1)
        ax_a.plot(time, m_bar, color = colors_theta[i_theta], label = models_name[0] + ' growth', linestyle = '--', linewidth = 3)
        ax_a.plot(time, m_bar_approx, color = 'grey', label = models_name[0] + ' growth', linestyle = ':', linewidth = 3)
        ax_m.plot(time, np.ones_like(time), color = 'black', linestyle = '-', linewidth = 1)
        ax_m.plot(time, m_bar, color = colors_theta[i_theta], label = models_name[0] + ' growth', linestyle = '--', linewidth = 3)
        ax_m.plot(time, m_bar_approx, color = 'grey', label = models_name[0] + ' growth', linestyle = ':', linewidth = 3)
    if(linear == 1):
        theory = ((k_on*M_r)/(N_A))*(1e8*time + 0.5*lambda_A*time**2)
    #ax_a.set_ylim(top=.5*N_r)

    #---- Entropy ----
    size_lim = 0
    if theta==1:
        size_lim=1000
    if theta==2:
        size_lim=1

    #---- Total Pop size ----
    total_pop = np.sum(clone_sizes, axis = 0)
    total_pop_active = total_pop - total_pop[0] + 1
    t_C = time[total_pop_active<C][-1] # Calculate time for reaching carrying capacity

    print('Delta_t = %.1f'%(t_C-t_act))

    filter_C = activations_times<t_C
    n_C = np.sum(filter_C)
    growth_time = time<t_C
    i_t_C = np.sum(growth_time)
    no_growth_time  = time>t_C
    i_t_C2 = np.sum(no_growth_time)
    clone_sizes_C = clone_sizes[filter_C, :]
    if(i_t_C<(len(time)-1)):
        print('freezing clone sizes...')
        final_sizes = clone_sizes_C[:, i_t_C].reshape((n_C, 1))
        final_sizes_list = final_sizes.tolist()
        clone_sizes_C[:,i_t_C+1:] = np.hstack([final_sizes_list]*i_t_C2)
    activations_times_C = activations_times[filter_C]
    energies_C = energies[filter_C]
    total_pop_active = np.sum(clone_sizes_C, axis = 0)
    ax_total_pop.plot(time, total_pop_active, color = colors_theta[i_theta])

    bcell_freqs = clone_sizes_C/total_pop_active
    entropy = -np.sum(bcell_freqs*np.log(bcell_freqs), axis = 0)

    if(t_C<t_theta):
        ax_H.plot(t_C, -S[Es[:-1]<E_t(t_theta, theta)][-1] , label = '%d'%theta, color = colors_theta[i_theta], linestyle = '', marker = 's', ms = 12)
    else:
        ax_H.plot(t_C, -S[Es[:-1]<E_t(t_C, theta)][-1]  , label = '%d'%theta, color = colors_theta[i_theta], linestyle = '', marker = 's', ms = 12)

    #---- Stackplots ----
    colors_muller = []
    min_bell_freq = np.min(bcell_freqs[:,-1])
    for c in range(int(len(clone_sizes_C[:,0]))):
        if bcell_freqs[c, -1]>(10*min_bell_freq/theta):
            colors_muller.append(colors_theta[i_theta])
        else:
            colors_muller.append('silver')
    ax_clones2.stackplot(time, bcell_freqs, colors = colors_muller);

    days_plot = np.linspace(Tf-1.5, Tf, 3)

    angles = np.random.random(len(energies_C))*2*np.pi
    energies_C = energies_C - (E_ms)
    radious = ((energies_C/(30)))*4

    for i_plot in range(len(days_plot)):
        for i_c in range(len(energies_C)):
            if (activations_times_C[i_c]<=days_plot[i_plot]):
                #print(i_c, positions[i_c], int(days_plot[i_plot]*len(time)/8)-1, time[int(days_plot[i_plot]*len(time)/8)-1], clone_sizes[i_c, int(days_plot[i_plot]*len(time)/8)-1], activations_times[i_c])
                #circle = plt.Circle(positions[i_c], np.sqrt(bcell_freqs[i_c, int(days_plot[i_plot]*len(time)/Tf)-1]/(np.pi*Tf)), color = colors_theta[i_theta], alpha = 1-(energies[i_c]/np.max(energies)))
                circle = plt.Circle((radious[i_c]*np.cos(angles[i_c]), radious[i_c]*np.sin(angles[i_c])), np.sqrt(clone_sizes_C[i_c, int(days_plot[i_plot]*len(time)/Tf)-1]/(C*np.pi)), color = colors_theta[i_theta], alpha = .8)
                ax_clones[i_plot].add_patch(circle)

        circle = plt.Circle((0, 0), 3, edgecolor = colors_theta[i_theta], facecolor = 'none')
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

    my_plot_layout(ax = ax_clones2, ticks_labelsize=34)
    ax_clones2.set_xlim(0, 6)
    ax_clones2.set_ylim(0, 1)
    ax_clones2.set_yticks([])
        
    my_plot_layout(ax = ax_a, yscale = 'log', xlabel = 'Time', ylabel = r'$\bar m$')
    ax_a.set_xlim(left = 3.5, right = Tf)
    ax_a.set_ylim(top=2*N_r, bottom = 1e-6)
    fig_a.savefig('../../Figures/1_Dynamics/Trajectories/activation_theta-%.1f.pdf'%theta)

    fig_clones.savefig('../../Figures/1_Dynamics/Trajectories/B_cell_clones_theta-%.1f.pdf'%(theta), dpi = 10)
    fig_clones2.savefig('../../Figures/1_Dynamics/Trajectories/B_cell_clones_2_theta-%.1f.pdf'%(theta), dpi = 10)

    ax_time.scatter(activations_times_C, np.exp(energies_C + (E_ms)), color = colors_theta[i_theta], alpha = .8)
    ax_time.hlines([Kd_pr, Kd_r, Kd_theta], 3.5, Tf, linestyle = ['-', '--', ':'], color = ['black', 'gray', colors_theta[i_theta]])

    ax_clone_size.scatter(clone_sizes_C[:,-1], np.exp(energies_C + (E_ms)), color = colors_theta[i_theta], alpha = .8)
    ax_clone_size.hlines([Kd_pr, Kd_r, Kd_theta], 1, 1e4, linestyle = ['-', '--', ':'], color = ['black', 'gray', colors_theta[i_theta]])

    Kds_C = np.exp(energies_C)
    NC = 1-np.array([np.product(1-1/(1+(Kds_C/((1e12*(clone_sizes_C[:,t]-1))/N_A)))) for t in np.arange(len(time))])
    ax_NC.plot(time, NC, color = colors_theta[i_theta])

my_plot_layout(ax = ax_m, yscale = 'log', xlabel = 'Time', ylabel = r'$\bar m$')
ax_m.set_xlim(left = 3., right = Tf)
ax_m.set_ylim(top=2*N_r, bottom = 1e-6)
fig_m.savefig('../../Figures/1_Dynamics/Trajectories/activation.pdf')

my_plot_layout(ax = ax_total_pop, yscale = 'linear', xlabel = 'Time', ylabel = r'$N_t$')
fig_total_pop.savefig('../../Figures/1_Dynamics/Trajectories/total_pop.pdf') 

my_plot_layout(ax = ax_H, yscale = 'linear', xlabel = 'Time', ylabel = r'$D_{KL}$')
ax_H.set_xlim(left = 3.5, right = Tf)
#ax_H.set_ylim(top=2*N_r, bottom = 1e-6)
fig_H.savefig('../../Figures/1_Dynamics/Trajectories/entropy.pdf') 

my_plot_layout(ax = ax_time, yscale = 'log', xscale = 'linear', ylabel = r'$K_d$', xlabel = r'times')
ax_time.set_xlim(left = 3.5, right = Tf+0.5)
#ax_time.set_ylim(top=2*N_r, bottom = 1e-6)
fig_time.savefig('../../Figures/1_Dynamics/Trajectories/times.pdf') 

my_plot_layout(ax = ax_clone_size, yscale = 'log', xscale = 'log', ylabel = r'$K_d$', xlabel = r'Clone size')
#ax_clone_size.set_xlim(left = 3.5, right = Tf)
#ax_clone_size.set_ylim(top=2*N_r, bottom = 1e-6)
fig_clone_size.savefig('../../Figures/1_Dynamics/Trajectories/clone_sizes.pdf') 

my_plot_layout(ax = ax_NC, yscale = 'log', xscale = 'linear', xlabel = r'times [days]', ylabel = r'Neutralization capacity')
ax_NC.set_xlim(left = 3.5, right = Tf+0.5)
#ax_NC.set_ylim(top=2*N_r, bottom = 1e-6)
fig_NC.savefig('../../Figures/1_Dynamics/Trajectories/Neutralization.pdf') 
    

