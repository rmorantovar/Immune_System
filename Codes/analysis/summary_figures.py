import sys
sys.path.append('../library/')
from Immuno_models import*
plt.rcParams['text.usetex'] = True


Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#fig_beta, ax_beta = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
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
antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'


transparency_q = [1, 1, .3, 0]
colors_theta = ['darkred', 'olive', 'darkblue']
colors_theta2 = ['tab:red', 'tab:olive', 'tab:blue']

colors_R = [['tab:red', 'tab:red', 'tab:red', 'darkred'], ['tab:olive', 'tab:olive', 'tab:olive', 'olive'], ['tab:blue', 'tab:blue', 'tab:blue', 'darkblue']]
colors_R = [['tab:red', 'tab:red', 'darkred', 'darkred'], ['tab:olive', 'tab:olive', 'olive', 'olive'], ['tab:blue', 'tab:blue', 'darkblue', 'darkblue']]

energy_models = ['MJ']
models_name = ['exponential', 'linear', ]
growth_models = [0] #, 1]

L=len(antigen)
print('L=%.d'%L)

N_r = 1e7
N_r = 1e8
N_r = 2e11

T0 = 0
Tf = 7
Tf = 3.6
dT = 0.1
days = np.linspace(2, Tf, 5)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
lambda_A = 8 #days^-1
k_pr = 2 # (M*hour)^-1
#k_pr = 180 # (M*hour)^-1
k_pr = k_pr*24 #(M*days)^-1
thetas = [1, 1.5, 2]
beta = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 2e5
E_ms = -27.63

print('K_d_ms=%.1e'%np.exp(E_ms))

print('max_u = %.2e'%(k_on*np.exp(Tf*lambda_A)/N_A))

print('k_pr/k_on = %.1e'%(k_on/k_pr)**(-1))

print('final antigen concetration = %.2e'%(np.exp(lambda_A*Tf)/1e3))
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
beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
Ks_r = np.exp(E_r)
print('beta_r = %.2f'%beta_r)

E_pr = Es[:-1][Ks<(k_pr/k_on)][-1]
Kd_pr = np.exp(E_pr)
beta_pr = betas[:-1][Ks<Kd_pr][-1]
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
#----------------------------------------------------------------
# points = np.array([Ks, betas[:-1]]).T.reshape(-1, 1, 2)
# segments = np.concatenate([points[:-1], points[1:]], axis=1)
# norm = plt.Normalize(0, 5)
# lc = LineCollection(segments, cmap='turbo_r', norm=norm)
# # Set the values used for colormapping
# lc.set_array(betas)
# lc.set_linewidth(5)
# line = ax_beta.add_collection(lc)
# # fig_beta.colorbar(line, ax=ax_beta)

# my_plot_layout(ax=ax_beta, yscale = 'linear', xscale = 'log', ticks_labelsize = 38)
# ax_beta.set_xticks([])
# ax_beta.set_yticks(np.arange(0, 6))
# ax_beta.set_ylim(top = 5, bottom = -.2)
# ax_beta.set_xlim(right = 1e-3)
# fig_beta.savefig('../../Figures/_Summary/beta.pdf')
#----------------------------------------------------------------

fig_H, ax_H = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

for i_theta, theta in enumerate(thetas):

    E_t = lambda t, theta:lambda_A*t/theta - np.log((lambda_A*N_A)/(k_on*N_c))/theta + np.log(k_pr/k_on) 
    
    beta_theta = betas[betas>theta][-1]
    E_theta = Es[betas>theta][-1]
    Kd_theta = np.exp(E_theta)

    delta_t_theta = (E_theta-E_pr)*theta/lambda_A
    t_theta = t_prime + delta_t_theta

    time1 = np.linspace(0, t_theta, 100)
    time2 = np.linspace(t_theta, Tf, 100)

    ax_H.plot(time1, -1*np.ones_like(time1)*S[Es[:-1]<E_theta][-1]-2.5, color = colors_theta[i_theta])
    ax_H.plot(time2, np.array([-S[Es[:-1]<E][-1] for E in E_t(time2, theta)]) -2.5, color = colors_theta[i_theta])

    fig_P_act, ax_P_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_Q_act, ax_Q_act = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_Q_act2, ax_Q_act2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_m, ax_m = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    

    ax_Q_act.plot(Ks, Q0*N_r, alpha = transparency_q[0], color = 'grey', linewidth = 5, linestyle = '-')
    ax_Q_act2.plot(Ks, Q0*N_r, alpha = transparency_q[0], color = 'grey', linewidth = 5, linestyle = '-')    

    print('theta = %.1f'%theta)
    #--------------------------m(t)---------------------------
    u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, theta, lambda_A, N_c, dE)
    M_r = N_r*N_c*np.sum(Q0*p_a*dE)
    m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])
    m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))

    ax_m.plot(time, m_bar, linewidth = 4, linestyle = '-', color = colors_theta[i_theta])
    ax_m.plot(time, m_bar_approx, linewidth = 3, linestyle = '--', color = 'black')
    ax_m.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
    t_act = time[m_bar>1][0]
    #----------------------------------------------------------------   
    
    days = np.linspace(0, t_act, 4)
    days = np.linspace(0, Tf, 4)
    for n_t, t in enumerate(days[[-4, -3, -2, -1]]):
        u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t)/N_A, Es, theta, lambda_A, N_c, dE)
        #ax_P_act.hlines(r_a[0]*N_c/lambda_A, ax_P_act.get_xlim()[0], ax_P_act.get_xlim()[1], alpha = transparency_q[i_theta], color = colors_theta[i_theta], linestyle = ':' )
        ax_P_act.set_ylim(bottom = 1e-11, top = 2)
        #----------------------------------------------------------------
        #--------------------------Q_act(E, t)---------------------------
        if theta!=0:
            #--------------------------P_act(E, t)---------------------------
            ax_P_act.plot(Ks, P_act, alpha = transparency_q[0], color = colors_R[i_theta][n_t], linewidth = 5, linestyle = '-')
            ax_Q_act.plot(Ks, Q_act*N_r, alpha = transparency_q[0], color = colors_R[i_theta][n_t], linewidth = 5, linestyle = '-')
            #-------FOR Q0--------- 
            #ax_Q_act.vlines(np.exp(E_r), ax_Q_act.get_ylim()[0], N_r*Q0[Ks<np.exp(E_r)][-1], color = 'black', linestyle = 'dashed')
            #ax_Q_act.vlines(np.exp(E_theta), ax_Q_act.get_ylim()[0], N_r*Q0[Ks<np.exp(E_theta)][-1], color = 'grey', linestyle = 'dotted', linewidth = 4)       
            #ax_Q_act.hlines(N_r*Q0[Ks<np.exp(E_r)][-1], ax_Q_act.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')
            #---------------------- 

    for n_t, t in enumerate(time[::4]):
        u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t)/N_A, Es, theta, lambda_A, N_c, dE)
        # ------------ H ---------------
        dE_temp = dE[Q_act!=0]
        Q0_temp = Q0[Q_act!=0]
        Ks_temp = Ks[Q_act!=0]
        Q_act_temp = Q_act[Q_act!=0]
        D_KL_t = np.sum(dE_temp*(Q_act_temp/np.sum(Q_act_temp*dE_temp))*(np.log((Q_act_temp/np.sum(Q_act_temp*dE_temp)))-np.log(Q0_temp)))
        ax_H.plot(t, D_KL_t, marker = 's', ms = 12, color = colors_theta2[i_theta])
        
    u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t_act)/N_A, Es, theta, lambda_A, N_c, dE)
    dE_temp = dE[Q_act!=0]
    Q0_temp = Q0[Q_act!=0]
    Ks_temp = Ks[Q_act!=0]
    Q_act_temp = Q_act[Q_act!=0]
    D_KL_t = np.sum(dE_temp*(Q_act_temp/np.sum(Q_act_temp*dE_temp))*(np.log((Q_act_temp/np.sum(Q_act_temp*dE_temp)))-np.log(Q0_temp)))
    ax_H.plot(t_act, D_KL_t, marker = 's', ms = 16, color = colors_theta[i_theta], linestyle = '', label = r'$%.1f$'%theta)
    #ax_H.hlines(L*np.log(20) - np.sum(betas[:-1][betas[:-1]>theta]*dE[betas[:-1]>theta]) - 5, 3, Tf, color = colors_theta[i_theta])
    #ax_H.hlines(-np.log(Q0[Es[:-1]<E_theta][-1])-2.5, 3, Tf, color = colors_theta[i_theta], linestyle = '--')
    # -----------------------------

    ax_Q_act.set_ylim(bottom = 1e-11, top = 2*N_r)
    
    ax0.plot(Ks, p_a, color = colors_theta[i_theta], alpha = 1, linewidth = 5, linestyle = '-', label = '%.1f'%(theta))

    ax_Q_act2.plot(Ks, Q_act*N_r, alpha = transparency_q[0], color = colors_theta[i_theta], linewidth = 5, linestyle = '-')
    u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*days[-1])/N_A, Es, theta, lambda_A, N_c, dE)
    m = np.sum(dE*Q_act*N_r)

    print('# of activated lineages : %d'%m , m_bar[-1])

    ax_Q_act2.hlines([1, N_r], ax_Q_act2.get_xlim()[0], ax_Q_act2.get_xlim()[1], alpha = 1, color = 'black', linestyle = ':')
    
        #----------------------------------------------------------------

    my_plot_layout(ax=ax_P_act, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
    ax_P_act.set_xticks([])
    ax_P_act.set_xlim(right = 1e-3)
    ax_P_act.set_ylim(top = 1.5)
    fig_P_act.savefig('../../Figures/_Summary/R_clone_theta-%.1f.pdf'%theta)

    my_plot_layout(ax=ax_Q_act, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
    #ax_Q_act.set_xticks([])
    ax_Q_act.set_xlim(right = 1e-3)
    fig_Q_act.savefig('../../Figures/_Summary/QR_theta-%.1f.pdf'%theta)
    my_plot_layout(ax=ax_Q_act2, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
    ax_Q_act2.set_xlim(right = 1e-3)
    fig_Q_act2.savefig('../../Figures/_Summary/QR2_theta-%.1f.pdf'%theta)

    my_plot_layout(ax=ax_m, yscale = 'log', ticks_labelsize = 30)
    fig_m.savefig('../../Figures/_Summary/activation_rate_theta-%.1f.pdf'%theta)

my_plot_layout(ax=ax_H, yscale = 'linear', ticks_labelsize = 30)
#ax_H.set_ylim(bottom = -1, top = 18)
ax_H.set_xlim(left = 0, right = Tf)
ax_H.legend(fontsize = 28, title = r'$\theta$', title_fontsize = 32, loc = 3)
fig_H.savefig('../../Figures/_Summary/entropy.pdf')

my_plot_layout(ax=ax0, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
ax0.legend(title = r'$\theta$', title_fontsize = 35, fontsize = 30)
#ax0.legend(fontsize = 30, title_fontsize=33)
fig0.savefig('../../Figures/_Summary/p_a.pdf')










