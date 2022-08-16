import sys
sys.path.append('../library/')
from Immuno_models import*
plt.rcParams['text.usetex'] = True


Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'


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
colors_N_r = ['darkred', 'olive', 'darkblue']
colors_N_r = ['tab:red', 'tab:olive', 'tab:blue']

colors_R = [['tab:red', 'tab:red', 'tab:red', 'darkred'], ['tab:olive', 'tab:olive', 'tab:olive', 'olive'], ['tab:blue', 'tab:blue', 'tab:blue', 'darkblue']]
colors_R = [['tab:red', 'tab:red', 'darkred', 'darkred'], ['tab:olive', 'tab:olive', 'olive', 'olive'], ['tab:blue', 'tab:blue', 'darkblue', 'darkblue']]

energy_models = ['MJ']
models_name = ['exponential', 'linear', ]
growth_models = [0] #, 1]

L=len(antigen)
print('L=%.d'%L)

N_rs = [1e6, 1e7, 1e8]
T0 = 3
Tf = 8
#Tf = 9
dT = 0.05
days = np.linspace(2, Tf, 5)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_pr = .5 # hour^-1
k_pr = k_pr*24 #days^-1

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

Es, dE, Q0, betas = calculate_Q0(0.01, 50, 200000, PWM_data, E_ms, L)
S = np.log(Q0)
Ks = np.exp(Es[:-1])

E_pr = Es[:-1][Ks<(k_pr/k_on)][-1]
Kd_pr = np.exp(E_pr)
t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))

#----------------------------------------------------------------

fig_H_t_C, ax_H_t_C = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_m_t_C, ax_m_t_C = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_t_C, ax_t_C = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

fig_H_t_act, ax_H_t_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_m_t_act, ax_m_t_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_t_act, ax_t_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

fig_E_act, ax_E_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_Delta_E_C, ax_Delta_E_C = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_best_E_C, ax_best_E_C = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

print('Initiating ...')
for i_N_r, N_r in enumerate(N_rs):

    print('Nr=%.0e'%N_r)

    beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
    E_r = Es[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
    Ks_r = np.exp(E_r)
    print('beta_r = %.2f'%beta_r)
    
    thetas = np.linspace(1, 2.8, 10)


    for i_theta, theta in enumerate(thetas):

        E_t = lambda t, theta:lambda_A*t/theta - np.log((lambda_A*N_A)/(k_on*N_c))/theta + np.log(k_pr/k_on) 
        
        beta_theta = betas[betas>theta][-1]
        E_theta = Es[betas>theta][-1]
        Kd_theta = np.exp(E_theta)

        delta_t_theta = (E_theta-E_pr)*theta/lambda_A
        t_theta = t_prime + delta_t_theta

        time1 = np.linspace(0, t_theta, 100)
        time2 = np.linspace(t_theta, Tf, 100)

        E_act = 0
        if(E_theta<E_r):
            E_act = E_r
        else:
            E_act = E_theta

        ax_E_act.plot(theta, E_act, marker = 's', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        #ax_H.plot(time1, -1*np.ones_like(time1)*S[Es[:-1]<E_theta][-1]-2.5, color = colors_N_r[i_N_r])
        #ax_H.plot(time2, np.array([-S[Es[:-1]<E][-1] for E in E_t(time2, theta)]) -2.5, color = colors_N_r[i_N_r])

        #--------------------------m(t)---------------------------
        u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, theta, lambda_A, N_c, dE)
        m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])

        t_act = time[m_bar>1][0]
        t_C = t_act + 1.23

        u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t_C)/N_A, Es, theta, lambda_A, N_c, dE)
        m_t_C = np.sum(dE*Q_act*N_r)

        u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t_act)/N_A, Es, theta, lambda_A, N_c, dE)
        m_t_act = np.sum(dE*Q_act*N_r)

        #print('# of activated lineages : %d'%m , m_bar[-1])
        ax_m_t_C.plot(theta, m_t_C, marker = 's', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        ax_t_C.plot(theta, t_C, marker = 's', ms = 12, color = colors_N_r[i_N_r], linestyle = '')

        ax_m_t_act.plot(theta, m_t_act, marker = 's', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        ax_t_act.plot(theta, t_act, marker = 's', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        
        #----------------------------------------------------------------
        
        if((i_N_r==1) and (i_theta%3==0)):

            fig_Q_act, ax_Q_act = plt.subplots(figsize=(20,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        
            ax_Q_act.plot(Ks, Q0*N_r, alpha = transparency_q[0], color = 'grey', linewidth = 5, linestyle = '-')

            days = np.linspace(2, t_C, 4)

            for n_t, t in enumerate(days[[-4, -3, -2, -1]]):
                u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t)/N_A, Es, theta, lambda_A, N_c, dE)
                #----------------------------------------------------------------
                #--------------------------Q_act(E, t)---------------------------
                ax_Q_act.plot(Ks, Q_act*N_r, alpha = transparency_q[0], color = colors_R[i_N_r][n_t], linewidth = 5, linestyle = '-')
                #-------FOR Q0--------- 
                ax_Q_act.vlines(np.exp(E_r), ax_Q_act.get_ylim()[0], N_r*Q0[Ks<np.exp(E_r)][-1], color = 'black', linestyle = 'dashed')
                ax_Q_act.vlines(np.exp(E_theta), ax_Q_act.get_ylim()[0], N_r*Q0[Ks<np.exp(E_theta)][-1], color = 'grey', linestyle = 'dotted', linewidth = 4)       
                #---------------------- 

        #------------------------------At activation-----------------------------------
        u_on, p_a, P_act, Q_act_t_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t_act)/N_A, Es, theta, lambda_A, N_c, dE)
        dE_temp = dE[Q_act_t_act!=0]
        Q0_temp = Q0[Q_act_t_act!=0]
        Ks_temp = Ks[Q_act_t_act!=0]
        Q_act_t_act_temp = Q_act_t_act[Q_act_t_act!=0]
        D_KL_t_act = np.sum(dE_temp*(Q_act_t_act_temp/np.sum(Q_act_t_act_temp*dE_temp))*(np.log((Q_act_t_act_temp/np.sum(Q_act_t_act_temp*dE_temp)))-np.log(Q0_temp)))
        ax_H_t_act.plot(theta, D_KL_t_act, marker = 's', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        #ax_H.hlines(L*np.log(20) - np.sum(betas[:-1][betas[:-1]>theta]*dE[betas[:-1]>theta]) - 5, 3, Tf, color = colors_theta[i_theta])
        #ax_H.hlines(-np.log(Q0[Es[:-1]<E_theta][-1])-2.5, 3, Tf, color = colors_theta[i_theta], linestyle = '--')
        Es_act = Es[:-1][(Q_act_t_act)==np.max(Q_act_t_act)]
        ax_E_act.plot(theta, Es_act, marker = '^', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        #------------------------------At carrying capacity----------------------------------
        u_on, p_a, P_act, Q_act_t_C = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t_C)/N_A, Es, theta, lambda_A, N_c, dE)
        dE_temp = dE[Q_act_t_C!=0]
        Q0_temp = Q0[Q_act_t_C!=0]
        Ks_temp = Ks[Q_act_t_C!=0]
        Q_act_t_C_temp = Q_act_t_C[Q_act_t_C!=0]
        D_KL_t_C = np.sum(dE_temp*(Q_act_t_C_temp/np.sum(Q_act_t_C_temp*dE_temp))*(np.log((Q_act_t_C_temp/np.sum(Q_act_t_C_temp*dE_temp)))-np.log(Q0_temp)))
        ax_H_t_C.plot(theta, D_KL_t_C, marker = 's', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        #ax_H.hlines(L*np.log(20) - np.sum(betas[:-1][betas[:-1]>theta]*dE[betas[:-1]>theta]) - 5, 3, Tf, color = colors_theta[i_theta])
        #ax_H.hlines(-np.log(Q0[Es[:-1]<E_theta][-1])-2.5, 3, Tf, color = colors_theta[i_theta], linestyle = '--')
        Es_C = Es[:-1][(N_r*Q_act_t_C)>1]
        Delta_E_C = Es_C[-1] - Es_C[0]
        ax_Delta_E_C.plot(theta, Delta_E_C, marker = 's', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        ax_best_E_C.plot(theta, Es_C[0], marker = '^', ms = 12, color = colors_N_r[i_N_r], linestyle = '')
        #----------------------------------------------------------------

        if((i_N_r==1) and (i_theta%3==0)):
            ax_Q_act.set_ylim(bottom = 1e-11, top = 2*N_r)
            my_plot_layout(ax=ax_Q_act, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
            #ax_Q_act.set_xticks([])
            ax_Q_act.set_xlim(right = 1e-3)
            fig_Q_act.savefig('../../Figures/9_Desing_principle/QR_theta-%.1f_Nr-%.0e.pdf'%(theta, N_r))
        
    k_act = np.min([Ks_r+k_on, Kd_theta*k_on])
    ax_t_C.plot(thetas, t_C+(2*(thetas-thetas[-1])/(k_pr+k_act)), color = colors_N_r[i_N_r])

my_plot_layout(ax=ax_m_t_C, yscale = 'log', ticks_labelsize = 30)
ax_m_t_C.set_ylim(bottom = .5)
fig_m_t_C.savefig('../../Figures/9_Desing_principle/m_t_C.pdf')

my_plot_layout(ax=ax_t_C, yscale = 'linear', ticks_labelsize = 30)
fig_t_C.savefig('../../Figures/9_Desing_principle/t_C.pdf')

my_plot_layout(ax=ax_H_t_C, yscale = 'linear', ticks_labelsize = 30)
#ax_H_t_C.set_ylim(bottom = -1, top = 18)
#ax_H_t_C.set_xlim(left = 2.5, right = Tf)
#ax_H_t_C.legend(fontsize = 28, title = r'$\theta$', title_fontsize = 32, loc = 3)
fig_H_t_C.savefig('../../Figures/9_Desing_principle/entropy_t_C.pdf')

my_plot_layout(ax=ax_m_t_act, yscale = 'log', ticks_labelsize = 30)
ax_m_t_act.set_ylim(top = ax_m_t_C.get_ylim()[1], bottom = .5)
fig_m_t_act.savefig('../../Figures/9_Desing_principle/m_t_act.pdf')

my_plot_layout(ax=ax_t_act, yscale = 'linear', ticks_labelsize = 30)
fig_t_act.savefig('../../Figures/9_Desing_principle/t_act.pdf')

my_plot_layout(ax=ax_H_t_act, yscale = 'linear', ticks_labelsize = 30)
#ax_H_t_act.set_ylim(bottom = -1, top = 18)
#ax_H_t_act.set_xlim(left = 2.5, right = Tf)
#ax_H_t_act.legend(fontsize = 28, title = r'$\theta$', title_fontsize = 32, loc = 3)
fig_H_t_act.savefig('../../Figures/9_Desing_principle/entropy_t_act.pdf')

my_plot_layout(ax=ax_E_act, yscale = 'linear', ticks_labelsize = 30)
fig_E_act.savefig('../../Figures/9_Desing_principle/E_act.pdf')

my_plot_layout(ax=ax_Delta_E_C, yscale = 'linear', ticks_labelsize = 30)
fig_Delta_E_C.savefig('../../Figures/9_Desing_principle/Delta_E_C.pdf')

my_plot_layout(ax=ax_best_E_C, yscale = 'linear', ticks_labelsize = 30)
fig_best_E_C.savefig('../../Figures/9_Desing_principle/Best_E_C.pdf')




