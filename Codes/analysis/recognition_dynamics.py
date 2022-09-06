import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 1
N_rs = [1e8]
T0 = 0
Tf = 10
Tf_sim = 6.5
#Tf = 10
dT = 0.01
lambda_A = 6
k_pr = 1 # hour^-1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

ns = [2.2, 2.0, 1.8, 1.5]#, 1]
ns = [1.4, 1.8, 2.2]
ns = [4, 2, 1]

transparency_n = [1]

colors_theta = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
colors_theta = ['tab:cyan','green', 'tab:orange', 'orange', 'darkred']
colors_R = [['tab:purple', 'tab:cyan', 'tab:cyan'], ['tab:blue', 'tab:green', 'tab:green'], ['tab:red', 'tab:orange', 'tab:orange']]

lambda_B = lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
#N_c = 1e5
E_ms = -27.63
C = 3e4

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

# antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
# antigen = 'FMLFMAVFVMTSWYC'
# antigen = 'FTSENAYCGR'
# antigen = 'TACNSEYPNTTK'
# antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'
antigen = 'TACNSEYPNTTKCGRWYC'
L=len(antigen)
#----------------------------------------------------------------
model = 'TCRen'
model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data = get_motif(antigen, model, Text_files_path)
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

avg_E = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]) + E_ms
var_E = np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])
beta_1, E_1, Kd_1 = get_n_properties(betas, Q0, Es, dE, 1)
#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))

#--------------------------Loops--------------------------

for N_r in N_rs:
    print('N_r = %.0e'%N_r)
    fig_Q0, ax_Q0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
    print('beta_r = %.1f'%beta_r)

    ax_Q0.plot(Kds, Q0, alpha = 1, color = 'grey', linewidth = 5, linestyle = '-')
    ax_Q0.plot(Kds, P_min_e(N_r, avg_E, var_E, Es[:-1], dE), linestyle = '--', marker = '',  color = 'black', ms = 2, linewidth = 4)
    #ax_Q0.plot(Kds, Kds**(beta_r)/(Ks[P_min_e(N_r, avg_E, var_E, Es[:-1], dE)==np.max(P_min_e(N_r, avg_E, var_E, Es[:-1], dE))]**(beta_r))*np.max(P_min_e(N_r, avg_E, var_E, Es[:-1], dE)), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 4 )

    for i_n, n in enumerate(ns):

        print('n = %.2f...'%n)

        fig_R, ax_R = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_QR, ax_QR = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_QR2, ax_QR2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_m_bar, ax_m_bar = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

        beta_n, E_n, Kd_n = get_n_properties(betas, Q0, Es, dE, n)
        delta_t_n = (E_n-E_pr)*n/lambda_A
        t_n = t_prime + delta_t_n

        time1 = np.linspace(0, t_n, 100)
        time2 = np.linspace(t_n, Tf, 100)

        ax_QR.plot(Kds, Q0*N_r, alpha = transparency_n[0], color = 'grey', linewidth = 5, linestyle = '-')

        ax_QR2.plot(Kds, Q0, alpha = .8, color = 'grey', linewidth = 5, linestyle = '-')
        ax_QR2.vlines([Kd_r], 0, 0.355, alpha = 1, color = 'black', linestyle = '--', linewidth = 5)

        #--------------------------m_bar(t)---------------------------
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, n, lambda_A, N_c, dE)
        M_r = N_r*N_c*np.sum(Q0*p_a*dE)
        
        #m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])# Review this
        m_bar = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, n, lambda_A, N_c, dE)[3]*dE) for t in time]) 
        m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))

        ax_m_bar.plot(time, m_bar, linewidth = 4, linestyle = '-', color = colors_theta[i_n])
        ax_m_bar.plot(time, m_bar_approx, linewidth = 3, linestyle = '--', color = 'black')
        ax_m_bar.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
        t_act = time[m_bar>1][0]
        print('t_act: %.2f'%t_act)
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act+1.3))/N_A, Es, n, lambda_A, N_c, dE)
        m_f_expected = np.sum(N_r*QR*dE)
        print('Activated clones expected:%.d'%m_f_expected)

        days = [t_act-1.3, t_act, t_act+1.3]
        #days = np.linspace(1, Tf-0.5, 3)
        for n_t, t in enumerate(days):
            u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t)/N_A, Es, n, lambda_A, N_c, dE)
            #ax_P_act.hlines(r_a[0]*N_c/lambda_A, ax_P_act.get_xlim()[0], ax_P_act.get_xlim()[1], alpha = transparency_q[i_n], color = colors_theta[i_n], linestyle = ':' )
            ax_R.set_ylim(bottom = 1e-11, top = 2)
            #----------------------------------------------------------------
            #--------------------------Q_act(E, t)---------------------------
            if n!=0:
                #--------------------------P_act(E, t)---------------------------
                ax_R.plot(Kds, R, alpha = transparency_n[0], color = colors_R[i_n][n_t], linewidth = 5, linestyle = '-')
                ax_QR.plot(Kds, QR*N_r, alpha = transparency_n[0], color = colors_R[i_n][n_t], linewidth = 5, linestyle = '-')
                #-------FOR Q0--------- 
                ax_QR.vlines(Kd_r, 0, .5, color = 'black', linestyle = 'dashed')
                ax_QR.vlines([Kd_n, Kd_1], ax_QR.get_ylim()[0], N_r*Q0[Kds<np.exp(E_n)][-1], color = 'grey', linestyle = 'dotted', linewidth = 4)       
                #ax_Q_act.hlines(N_r*Q0[Ks<np.exp(E_r)][-1], ax_Q_act.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')

        ax_QR2.plot(Kds, QR/np.sum(QR*dE), alpha = transparency_n[0], color = colors_R[i_n][n_t], linewidth = 5, linestyle = '-')

        my_plot_layout(ax=ax_R, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
        ax_R.set_xticks([])
        ax_R.set_xlim(right = 1e-2)
        ax_R.set_ylim(top = 1.5, bottom = 1e-8)
        fig_R.savefig('../../Figures/_Summary/single/R_n-%.1f_Nr-%.0e_'%(n, N_r)+model+'.pdf')
        plt.close(fig_R)

        my_plot_layout(ax=ax_QR, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
        #ax_QR.set_xticks([])
        ax_QR.set_xlim(right = 1e-2) #use 1e-3 for other plots
        ax_QR.set_ylim(bottom = 1e-9, top = 1.5*N_r)
        fig_QR.savefig('../../Figures/_Summary/single/QR_n-%.1f_Nr-%.0e_'%(n, N_r)+model+'.pdf')
        plt.close(fig_QR)

        my_plot_layout(ax=ax_QR2, yscale = 'linear', xscale = 'log', ticks_labelsize = 30)
        #ax_QR2.set_xlim(right = 1e-3)
        fig_QR2.savefig('../../Figures/_Summary/single/QR2_n-%.1f_Nr-%.0e_'%(n, N_r)+model+'.pdf')
        plt.close(fig_QR2)

        my_plot_layout(ax=ax_m_bar, yscale = 'log', ticks_labelsize = 30)
        fig_m_bar.savefig('../../Figures/_Summary/single/m_bar_n-%.1f_Nr-%.0e_'%(n, N_r)+model+'.pdf')
        plt.close(fig_m_bar)
        
    my_plot_layout(ax=ax_Q0, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
    #ax_Q_act.set_xticks([])
    ax_Q0.set_xlim(right = 1e1) #use 1e-3 for other plots
    ax_Q0.set_ylim(bottom = 1e-9, top = 1.5)
    fig_Q0.savefig('../../Figures/_Summary/single/Q0_Nr-%.0e_'%(N_r)+model+'.pdf')
    plt.close(fig_Q0)




