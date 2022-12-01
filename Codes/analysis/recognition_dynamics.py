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
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
k_pr = k_pr*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*1000
#N_c = 1e5
#E_ms = -27.63
E_ms = -25
C = 1e4
AA = 1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = np.flip([1, 2, 2.5, 3, 4])
t_fs = np.flip([2.48, 3.52, 4, 4.15, 4.61])

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

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])
color_list = np.array([my_blue, my_blue, my_blue, my_blue, my_blue])

#colors_kappa = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_kappa = np.flip(['tab:blue','tab:green','tab:red'])
colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

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
#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'CMFILVWYAGTSQNEDHRKPFMRTPTRMCW' #L=30
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
Es, dE, Q0, betas = calculate_Q0(0.05, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])
beta_1, E_1, Kd_1 = get_kappa_properties(betas, Q0, Es, dE, 1)
#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_a = %.2f'%beta_pr, 'K_a = %.2e'%Kd_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
#--------------------------Loops--------------------------
print('Loops...')
for N_r in N_rs:
    print('________')
    print('N_r = %.0e'%N_r)
    fig_Q0, ax_Q0 = plt.subplots(figsize=(5,4), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
    print('beta_r = %.1f'%beta_r, 'K_d_r = %.2e'%Kd_r)

    ax_Q0.plot(Kds, Q0, alpha = 1, color = 'grey', linewidth = 5, linestyle = '--')
    ax_Q0.plot(Kds, P_min_e(N_r, avg_E, var_E, Es[:-1], dE), linestyle = '--', marker = '',  color = 'black', ms = 2, linewidth = 4)
    #ax_Q0.plot(Kds, Kds**(beta_r)/(Ks[P_min_e(N_r, avg_E, var_E, Es[:-1], dE)==np.max(P_min_e(N_r, avg_E, var_E, Es[:-1], dE))]**(beta_r))*np.max(P_min_e(N_r, avg_E, var_E, Es[:-1], dE)), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 4 )

    
    for i_kappa, kappa in enumerate(kappas):
        print('--------')

        #fig_QR_all_f, ax_QR_all_f = plt.subplots(figsize=(16*0.33,4), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.12, 'top': 0.96})
        fig_QR_all_f, ax_QR_all_f = plt.subplots(figsize=(4.5*2,3*0.8*2), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.12, 'top': 0.96})

        #ax_QR_all.plot(Kds, Q0*N_r, alpha = 1, color = 'grey', linewidth = 5, linestyle = '--')
        #ax_QR_all.vlines(Kd_r, Q0[Kds==Kd_r]*N_r, N_r, alpha = 1, color = 'grey', linewidth = 3, linestyle = ':')
        ax_QR_all_f.plot(Kds, Q0*N_r, alpha = 1, color = 'grey', linewidth = 5, linestyle = '--')
        #ax_QR_all_f.vlines(Kd_r, Q0[Kds==Kd_r]*N_r, N_r, alpha = 1, color = 'grey', linewidth = 3, linestyle = ':')
        #ax_QR_all_f.vlines(Kd_pr, Q0[Kds==Kd_pr]*N_r, N_r, alpha = 1, color = 'grey', linewidth = 1, linestyle = '-.')

        fig_R, ax_R = plt.subplots(figsize=(8,6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_QR, ax_QR = plt.subplots(figsize=(8,6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_QR2, ax_QR2 = plt.subplots(figsize=(8,6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_m_bar, ax_m_bar = plt.subplots(figsize=(8,6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

        beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
        print('beta_kappa = %.2f...'%kappa, 'K_kappa = %.2e...'%Kd_kappa)
        delta_t_n = (E_kappa-E_pr)*kappa/lambda_A
        t_n = t_prime + delta_t_n

        time1 = np.linspace(0, t_n, 100)
        time2 = np.linspace(t_n, Tf, 100)

        ax_QR.plot(Kds, Q0*N_r, alpha = transparency_n[0], color = 'grey', linewidth = 5, linestyle = '-')

        ax_QR2.plot(Kds, Q0, alpha = .8, color = 'grey', linewidth = 5, linestyle = '-')
        ax_QR2.vlines([Kd_r], 0, 0.355, alpha = 1, color = 'black', linestyle = '--', linewidth = 5)

        #--------------------------m_bar(t)---------------------------
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, kappa, lambda_A, N_c, dE)
        M_r = N_r*N_c*np.sum(Q0*p_a*dE)
        
        #m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])# Review this
        m_bar = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time]) 
        m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))

        ax_m_bar.plot(time, m_bar, linewidth = 4, linestyle = '-', color = colors_kappa[4-i_kappa])
        ax_m_bar.plot(time, m_bar_approx, linewidth = 3, linestyle = '--', color = 'black')
        ax_m_bar.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
        t_act = time[m_bar>1][0]
        print('t_act: %.2f'%t_act)
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act+1.3))/N_A, Es, kappa, lambda_A, N_c, dE)
        m_f_expected = np.sum(N_r*QR*dE)
        print('Activated clones expected:%.d'%m_f_expected)

        days = [t_act-2*1.2, t_act-1.2, t_act, t_fs[i_kappa]]
        #days = np.linspace(1, Tf-0.5, 3)


        for i_t, t in enumerate(days):

            u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*t)/N_A, Es, kappa, lambda_A, N_c, dE)
            #----------------------------------------------------------------
            #--------------------------QR_all(E, t)---------------------------
            if i_t==2:
                ax_R.plot(Kds, R, color = colors_R[4-i_kappa][i_t], linewidth = 5, linestyle = '-')
                
                ax_QR.plot(Kds, QR*N_r, color = colors_R[4-i_kappa][i_t], linewidth = 5, linestyle = '-')
                
                #ax_QR_all.plot(Kds, QR*N_r, color = colors_kappa[4-i_kappa], linewidth = 5, linestyle = '-', label = r'$%d$'%(kappa), alpha = .8)
                
                ax_QR_all_f.plot(Kds, QR*N_r, color = colors_kappa[4-i_kappa], linewidth = 5, linestyle = '--', alpha = .6)

                #ax_QR_all.vlines(Kds[QR==np.max(QR)], 1e-9, (Q0*N_r)[QR==np.max(QR)], color = colors_kappa[4-i_kappa], linewidth = 3, linestyle = '--')
                ax_QR_all_f.vlines(Kds[QR==np.max(QR)], 1e-9, (Q0*N_r)[QR==np.max(QR)], color = 'black', linewidth = 4, linestyle = '--', alpha = 1)
                
                #ax_QR_all.plot(Kds[Kds==Kd_kappa], (Q0*N_r)[Kds==Kd_kappa], markerfacecolor = colors_kappa[4-i_kappa], marker = '*', ms = 22, markeredgecolor='black', alpha = .6)
                #ax_QR_all.plot(Kds[QR==np.max(QR)], (Q0*N_r)[QR==np.max(QR)], markerfacecolor = colors_kappa[4-i_kappa], marker = 'D', ms = 16, markeredgecolor='black', alpha = .6)

                ax_QR_all_f.plot(Kds[QR==np.max(QR)], (Q0*N_r)[QR==np.max(QR)], markerfacecolor = colors_kappa[4-i_kappa], marker = 'D', ms = 16, markeredgecolor='black', alpha = .6)
                ax_QR_all_f.plot(Kds[Kds==Kd_kappa], (Q0*N_r)[Kds==Kd_kappa], markerfacecolor = colors_kappa[4-i_kappa], marker = '*', ms = 22, markeredgecolor='black', alpha = .6)

            if i_t==3:

                ax_R.plot(Kds, R, color = colors_R[4-i_kappa][i_t], linewidth = 5, linestyle = '-', alpha = .6)
                ax_QR.plot(Kds, QR*N_r, color = colors_R[4-i_kappa][i_t], linewidth = 5, linestyle = '-', alpha = .6)

                ax_QR_all_f.plot(Kds, QR*N_r, color = colors_kappa[4-i_kappa], linewidth = 5, linestyle = '-', alpha = .8)
                ax_QR_all_f.plot(Kds[QR==np.max(QR)], (Q0*N_r)[QR==np.max(QR)], markerfacecolor = colors_kappa[4-i_kappa], marker = 'o', ms = 16, alpha = .6, markeredgecolor='black')
                #if kappa!=1:
                ax_QR_all_f.vlines(Kds[QR==np.max(QR)], 1e-9, (Q0*N_r)[QR==np.max(QR)], color = 'black', linewidth = 4, linestyle = '--', alpha = 1)
            

                
            else:
                #--------------------------R(E, t) and QR(E, t)---------------------------
                ax_R.plot(Kds, R, alpha = transparency_n[0], color = colors_R[4-i_kappa][i_t], linewidth = 4, linestyle = '--')
                ax_QR.plot(Kds, QR*N_r, alpha = transparency_n[0], color = colors_R[4-i_kappa][i_t], linewidth = 4, linestyle = '--')
                #-------FOR Q0--------- 
                #ax_QR.vlines(Kd_r, 0, .5, color = 'black', linestyle = 'dashed')
                #ax_QR.vlines([Kd_kappa, Kd_1], ax_QR.get_ylim()[0], N_r*Q0[Kds<np.exp(E_n)][-1], color = 'grey', linestyle = 'dotted', linewidth = 4)       
                #ax_Q_act.hlines(N_r*Q0[Ks<np.exp(E_r)][-1], ax_Q_act.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')

        
        ax_QR2.plot(Kds, QR/np.sum(QR*dE), alpha = transparency_n[0], color = colors_R[4-i_kappa][i_t], linewidth = 5, linestyle = '-')

        my_plot_layout(ax=ax_R, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
        ax_R.set_xticks([])
        ax_R.set_xlim(right = 1e-2)
        ax_R.set_ylim(top = 1.5, bottom = 1e-8)
        fig_R.savefig('../../Figures/_Summary/affinity/L%d/R_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+model+'.pdf')
        plt.close(fig_R)

        my_plot_layout(ax=ax_QR, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
        #ax_QR.set_xticks([])
        ax_QR.set_xlim(right = 1e-2) #use 1e-3 for other plots
        ax_QR.set_ylim(bottom = 1e-9, top = 1.5*N_r)
        fig_QR.savefig('../../Figures/_Summary/affinity/L%d/QR_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+model+'.pdf')
        plt.close(fig_QR)

        my_plot_layout(ax=ax_QR2, yscale = 'linear', xscale = 'log', ticks_labelsize = 30)
        #ax_QR2.set_xlim(right = 1e-3)
        fig_QR2.savefig('../../Figures/_Summary/affinity/L%d/QR2_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+model+'.pdf')
        plt.close(fig_QR2)

        my_plot_layout(ax=ax_m_bar, yscale = 'log', ticks_labelsize = 30)
        fig_m_bar.savefig('../../Figures/_Summary/affinity/L%d/m_bar_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+model+'.pdf')
        plt.close(fig_m_bar)
    
        # my_plot_layout(ax=ax_QR_all, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
        # #ax_QR_all.set_xticks([])
        # ax_QR_all.set_xlim(right = 1e-0, left = 1e-11) #use 1e-3 for other plots
        # ax_QR_all.set_ylim(bottom = 1e-9, top = 1.5*N_r)
        # ax_QR_all.legend(title = r'$p$', title_fontsize = 32, fontsize = 30, loc = 4)
        # fig_QR_all.savefig('../../Figures/_Summary/affinity/L%d/QR_all_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+model+'.pdf')

        my_plot_layout(ax=ax_QR_all_f, yscale = 'log', xscale = 'log', ticks_labelsize = 34)
        #ax_QR_all_f.set_yticks([])
        ax_QR_all_f.set_xlim(right = 1e-0, left = 1e-11) #use 1e-3 for other plots
        ax_QR_all_f.set_ylim(bottom = 1e-9, top = 1.5*N_r)
        #ax_QR_all_f.legend(title = r'$p$', title_fontsize = 32, fontsize = 30)
        fig_QR_all_f.savefig('../../Figures/_Summary/affinity/L%d/QR_all_f_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+model+'.pdf')
        plt.close(fig_QR_all_f)

    my_plot_layout(ax=ax_Q0, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
    #ax_Q_act.set_xticks([])
    ax_Q0.set_xlim(right = 1e1) #use 1e-3 for other plots
    ax_Q0.set_ylim(bottom = 1e-13, top = 1.5)
    fig_Q0.savefig('../../Figures/_Summary/affinity/L%d/Q0_Nr-%.0e_'%(L, N_r)+model+'.pdf')
    plt.close(fig_Q0)




