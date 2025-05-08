import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/primary_immune_response/in/'

#--------------- PARAMETERS ---------------------
N_ens = 1
N_rs = [2e8]
T0 = 3
Tf = 10
Tf_sim = 6.5
#Tf = 10
dT = 0.01
lambda_A = 6
k_pr = 1 # hour^-1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2, 3]

my_red = np.array((228,75,41))/256.
my_purple = np.array((125,64,119))/256.
my_purple2 = np.array((116,97,164))/256.
my_green = np.array((125,165,38))/256.
my_blue = np.array((76,109,166))/256.
my_gold = np.array((215,139,45))/256.
my_brown = np.array((182,90,36))/256.
my_blue2 = np.array((80,141,188))/256.
my_yellow = np.array((246,181,56))/256.
my_green2 = np.array((158,248,72))/256.
my_cyan = 'tab:cyan'

antigen_color = my_yellow/256.

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_red, my_green, my_blue2, my_gold])

#colors_kappa = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_kappa = np.flip(['tab:blue','tab:green','tab:red'])
colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

lambda_B = lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5
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
antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'
#antigen = 'TACNSEYPNTTKCGRWYC'
L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
model = 'TCRen'
#model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, model, Text_files_path)
print('min_e_PWM=%.4f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

avg_E = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]) + E_ms
var_E = np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])
beta_1, E_1, Kd_1 = get_kappa_properties(betas, Q0, Es, dE, 1)
#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
#--------------------------Loops--------------------------
print('Loops...')
for N_r in N_rs:
    print('________')
    print('N_r = %.0e'%N_r)
    fig_QRT_all, ax_QRT_all = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
    print('beta_r = %.1f'%beta_r)

    ax_QRT_all.plot(Kds, Q0, alpha = 1, color = 'grey', linewidth = 5, linestyle = '--')
    #ax_QRT_all.vlines(Kd_r, Q0[Kds==Kd_r]*N_r, 1, alpha = 1, color = 'grey', linewidth = 2, linestyle = ':')
    for i_kappa, kappa in enumerate(kappas):
        print('--------')
        print('kappa = %.2f...'%kappa)

        fig_QRT, ax_QRT = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_QRT2, ax_QRT2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

        beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
        delta_t_n = (E_kappa-E_pr)*kappa/lambda_A
        t_n = t_prime + delta_t_n

        time1 = np.linspace(0, t_n, 100)
        time2 = np.linspace(t_n, Tf, 100)

        ax_QRT.plot(Kds, Q0, alpha = transparency_n[0], color = 'grey', linewidth = 5, linestyle = '-')

        ax_QRT2.plot(Kds, Q0, alpha = .8, color = 'grey', linewidth = 5, linestyle = '-')
        ax_QRT2.vlines([Kd_r], 0, 0.355, alpha = 1, color = 'black', linestyle = '--', linewidth = 5)

        #--------------------------m_bar(t)---------------------------
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, kappa, lambda_A, N_c, dE)
        M_r = N_r*N_c*np.sum(Q0*p_a*dE)
        
        #m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])# Review this
        m_bar = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time]) 
        m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))

        t_act = time[m_bar>1][0]
        print('t_act: %.2f'%t_act)
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act+1.3))/N_A, Es, kappa, lambda_A, N_c, dE)
        m_f_expected = np.sum(N_r*QR*dE)
        print('Activated clones expected:%.d'%m_f_expected)

        days = [t_act-2*1.5, t_act-1.5, t_act, t_act+1.5]
        #days = np.linspace(1, Tf-0.5, 3)
        for i_t, t in enumerate(days):

            u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t)/N_A, Es, kappa, lambda_A, N_c, dE)
            #----------------------------------------------------------------
            #--------------------------QR_all(E, t)---------------------------
            if i_t==3:
                QRT = QR*N_r*np.exp(-Es[:-1]*lambda_B*kappa/lambda_A)/np.exp(Es[:-1])
                ax_QRT.plot(Kds, QRT/np.sum(QRT*dE), alpha = transparency_n[0], color = colors_R[i_kappa][i_t], linewidth = 5, linestyle = '-')

                ax_QRT_all.plot(Kds, QRT/np.sum(QRT*dE), alpha = transparency_n[0], color = colors_kappa[i_kappa], linewidth = 5, linestyle = '-', label = r'$%d$'%(kappa))
                #ax_QRT_all.plot(Kds[QR==np.max(QR)], (Q0*N_r)[QR==np.max(QR)], alpha = transparency_n[0], color = colors_kappa[i_kappa], marker = 'o', ms = 12)
                #ax_QRT_all.vlines(Kds[Kds==Kd_kappa], 1e-9, (QR*N_r)[Kds==Kd_kappa], alpha = transparency_n[0], color = colors_kappa[i_kappa], linewidth = 1.5, linestyle = '--')

            else:
                #--------------------------R(E, t) and QR(E, t)---------------------------
                QRT = QR*N_r*np.exp(-Es[:-1]*lambda_B*kappa/lambda_A)/np.exp(Es[:-1])
                ax_QRT.plot(Kds, QRT/np.sum(QRT*dE), alpha = transparency_n[0], color = colors_R[i_kappa][i_t], linewidth = 4, linestyle = '--')
                #-------FOR Q0--------- 
                #ax_QRT.vlines(Kd_r, 0, .5, color = 'black', linestyle = 'dashed')
                #ax_QRT.vlines([Kd_kappa, Kd_1], ax_QRT.get_ylim()[0], N_r*Q0[Kds<np.exp(E_n)][-1], color = 'grey', linestyle = 'dotted', linewidth = 4)       
                #ax_Q_act.hlines(N_r*Q0[Ks<np.exp(E_r)][-1], ax_Q_act.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')

        QRT = QR*N_r*np.exp(-Es[:-1]*lambda_B*kappa/lambda_A)/np.exp(Es[:-1])
        ax_QRT2.plot(Kds, QRT/np.sum(QRT*dE), alpha = transparency_n[0], color = colors_R[i_kappa][i_t], linewidth = 5, linestyle = '-')

        my_plot_layout(ax=ax_QRT, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
        #ax_QRT.set_xticks([])
        ax_QRT.set_xlim(right = 1e-2) #use 1e-3 for other plots
        ax_QRT.set_ylim(bottom = 1e-9)
        fig_QRT.savefig('../../Figures/_Summary/size/QRT_kappa-%.1f_Nr-%.0e_'%(kappa, N_r)+model+'.pdf')
        plt.close(fig_QRT)

        my_plot_layout(ax=ax_QRT2, yscale = 'linear', xscale = 'log', ticks_labelsize = 30)
        #ax_QRT2.set_xlim(right = 1e-3)
        fig_QRT2.savefig('../../Figures/_Summary/size/QRT2_kappa-%.1f_Nr-%.0e_'%(kappa, N_r)+model+'.pdf')
        plt.close(fig_QRT2)
    
    my_plot_layout(ax=ax_QRT_all, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
    #ax_QRT_all.set_xticks([])
    ax_QRT_all.set_xlim(right = 1e-0, left = 1e-11) #use 1e-3 for other plots
    ax_QRT_all.set_ylim(bottom = 1e-9)
    ax_QRT_all.legend(title = r'$p$', title_fontsize = 32, fontsize = 30)
    fig_QRT_all.savefig('../../Figures/_Summary/size/QRT_all_Nr-%.0e_'%(N_r)+model+'.pdf')
    plt.close(fig_QRT_all)

    


