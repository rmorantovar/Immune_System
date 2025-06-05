import sys
sys.path.append('../my_lib/')
sys.path.append('../../my_lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/primary_response/in/'
output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures'
#--------------- PARAMETERS ---------------------
N_ens = 1
L_0s = [1e9]
T0 = 0
Tf = 8
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
b0 = 1e5*10
#b0 = 1e5
#E_ms = -27.63
E_ms = -24
C = 1e4
AA = 1

ps = [1, 2, 2.5, 3, 4]
ps = [1, 4]
#t_ten = [3.1, 4.37, 4, 4.5, 4.9]
t_ten = [3.1, 4.5]
alpha_p = [.6, 1]
transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])
color_list = np.array([my_blue, my_red, my_red, my_red, my_red])

#colors_p = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_p = np.flip(['tab:blue','tab:green','tab:red'])
colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

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
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])
beta_1, E_1, Kd_1 = get_p_properties(betas, Q0, Es, dE, 1)
#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_a = %.2f'%beta_pr, 'K_a = %.2e'%Kd_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*b0))
print('--------')
#--------------------------Loops--------------------------
print('Loops...')
for L_0 in L_0s:
    print('________')
    print('L_0 = %.0e'%L_0)
    #fig_Q0, ax_Q0 = plt.subplots(figsize=(5,4), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_Q0, ax_Q0 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
    print('beta_r = %.1f'%beta_r, 'K_d_r = %.2e'%Kd_r)

    ax_Q0.plot(Kds, Q0*L_0, alpha = 1, color = 'grey', linewidth = 5, linestyle = '--')
    # ax_Q0.plot(Kds, P_min_e(L_0, avg_E, var_E, Es[:-1], dE), linestyle = '--', marker = '',  color = 'black', ms = 2, linewidth = 3, alpha = .6)
    ax_Q0.plot(Kds, P_min_e_Q0(L_0, Q0, dE)*L_0, linestyle = ':', marker = '',  color = 'grey', ms = 2, linewidth = 3, alpha = 1)
    ax_Q0.vlines(Kd_r, 1e-6, 1e9, linestyle = ':', color = 'grey', linewidth = 2, alpha = .6)
    u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*time[0])/N_A, Es, 3, lambda_A, b0, dE)
    K_array_t = np.exp(time*(lambda_A/3))*(k_step/k_on)
    p_a = (1+Kd_r/(k_step/k_on))**-3
    p_t_E = u_on*p_a*np.exp(lambda_A*np.log(K_array_t/(k_step/k_on))/(lambda_A/3))*np.exp(-(1/lambda_A)*(np.exp(lambda_A*np.log(K_array_t/(k_step/k_on))/(lambda_A/3)) - 1))/(lambda_A/3)
    ax_Q0.plot(K_array_t, p_t_E*L_0, linestyle = ':', marker = '',  color = my_blue, ms = 2, linewidth = 3, alpha = 1)
    print(np.exp(time*(lambda_A/3)))
    #ax_Q0.plot(Kds, Kds**(beta_r)/(Ks[P_min_e(L_0, avg_E, var_E, Es[:-1], dE)==np.max(P_min_e(L_0, avg_E, var_E, Es[:-1], dE))]**(beta_r))*np.max(P_min_e(L_0, avg_E, var_E, Es[:-1], dE)), linestyle = '-', marker = '',  color = 'black', ms = 2, linewidth = 4 )

    my_plot_layout(ax=ax_Q0, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
    #ax_Q_act.set_xticks([])
    ax_Q0.set_xlim(right = 1.2) #use 1e-3 for other plots
    ax_Q0.set_ylim(bottom = 1e-15)#, top = 1.5*L_0)
    fig_Q0.savefig(output_plot + '/primary_response/_Summary/affinity/L%d/Q0_Nr-%.0e_'%(L, L_0)+model+'1.pdf')

    my_plot_layout(ax=ax_Q0, yscale = 'linear', xscale = 'log', ticks_labelsize = 30)
    #ax_Q_act.set_xticks([])
    ax_Q0.set_xlim(right = 1.2) #use 1e-3 for other plots
    ax_Q0.set_ylim(bottom = 1e-15)#, top = 1.5*L_0)
    fig_Q0.savefig(output_plot + '/primary_response/_Summary/affinity/L%d/Q0_Nr-%.0e_'%(L, L_0)+model+'2.pdf')

    plt.close(fig_Q0)

    fig_QR_all_f, ax_QR_all_f = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    ax_QR_all_f.plot(Kds, Q0*L_0, alpha = 1, color = 'grey', linewidth = 5, linestyle = '--')

    fig_R_all, ax_R_all = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

    for i_p, p in enumerate(ps):
        print('--------')    
        
        fig_R, ax_R = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_QR, ax_QR = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_QR2, ax_QR2 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_m_bar, ax_m_bar = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

        beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
        print('beta_p = %.2f...'%p, 'K_p = %.2e...'%Kd_p)
        delta_t_n = (E_p-E_pr)*p/lambda_A
        t_n = t_prime + delta_t_n

        time1 = np.linspace(0, t_n, 100)
        time2 = np.linspace(t_n, Tf, 100)

        ax_QR.plot(Kds, Q0*L_0, alpha = transparency_n[0], color = 'grey', linewidth = 5, linestyle = '--')

        #ax_QR2.plot(Kds, Q0, alpha = .8, color = 'grey', linewidth = 5, linestyle = '-')
        ax_QR2.vlines([Kd_r], 0, 0.355, alpha = 1, color = 'black', linestyle = '--', linewidth = 3)

        #--------------------------m_bar(t)---------------------------
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*time[0])/N_A, Es, p, lambda_A, b0, dE)
        M_r = L_0*b0*np.sum(Q0*p_a*dE)
        
        #m_bar = np.array([L_0*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*b0))/lambda_A))*Q0*dE)) for t in time])# Review this
        m_bar = np.array([np.sum(L_0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, b0, dE)[3]*dE) for t in time]) 
        m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))

        ax_m_bar.plot(time, m_bar, linewidth = 4, linestyle = '-', color = colors_p[i_p])
        ax_m_bar.plot(time, m_bar_approx, linewidth = 3, linestyle = '--', color = 'black')
        ax_m_bar.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
        t_act = time[m_bar>1][0]
        print('t_act: %.2f'%t_act)
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t_act+1.3))/N_A, Es, p, lambda_A, b0, dE)
        m_f_expected = np.sum(L_0*QR*dE)
        print('Activated clones expected:%.d'%m_f_expected)

        days = [t_act-2*1.2, t_act-1.2, t_act, 1*(t_ten[i_p]-t_act)/2 + t_act, 2*(t_ten[i_p]-t_act)/2 + t_act,  4*(t_ten[i_p]-t_act)/2 + t_act]
        #days = np.linspace(1, Tf-0.5, 3)

        for i_t, t in enumerate(days):

            u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*t)/N_A, Es, p, lambda_A, b0, dE)
            #----------------------------------------------------------------
            #--------------------------QR_all(E, t)---------------------------
            if i_t==2:

                print('K_p = %.2e'%Kd_p)
                ax_R.plot(Kds, R, color = colors_R[i_p][i_t], linewidth = 5, linestyle = '-', alpha = alpha_p[i_p])
                
                ax_QR.plot(Kds, QR*L_0, color = colors_R[i_p][i_t], linewidth = 5, linestyle = '-', alpha = alpha_p[i_p])
                ax_QR.plot(Kds[Kds==Kd_r], (Q0*L_0)[Kds==Kd_r], markerfacecolor = 'grey', marker = '*', ms = 22, markeredgecolor='black', alpha = 0.9, zorder = 20)
                ax_QR.plot(Kds[Kds==Kd_p], (Q0*L_0)[Kds==Kd_p], markerfacecolor = colors_p[i_p], marker = 'o', ms = 18, markeredgecolor='black', alpha = alpha_p[i_p], zorder = 20)

                #ax_QR_all.plot(Kds, QR*L_0, color = colors_p[i_p], linewidth = 5, linestyle = '-', label = r'$%d$'%(p), alpha = .8)
                
                ax_QR_all_f.plot(Kds, QR*L_0, color = colors_p[i_p], linewidth = 5, linestyle = '-', alpha = alpha_p[i_p])
               #ax_QR_all_f.plot(Kds[QR==np.max(QR)], (Q0*L_0)[QR==np.max(QR)], markerfacecolor = 'white', marker = 'D', ms = 22, markeredgecolor='black', alpha = 1)
                ax_QR_all_f.plot(Kds[Kds==Kd_r], (Q0*L_0)[Kds==Kd_r], markerfacecolor = colors_p[i_p], marker = '*', ms = 22, markeredgecolor='black', alpha = .7, zorder = 20)
                ax_QR_all_f.plot(Kds[Kds==Kd_p], (Q0*L_0)[Kds==Kd_p], markerfacecolor = colors_p[i_p], marker = 'o', ms = 18, markeredgecolor='black', alpha = alpha_p[i_p], zorder = 20)
                #ax_QR_all_f.plot(Kds[Kds==Kds[0]], (Q0*L_0)[Kds==Kds[0]], markerfacecolor = colors_p[i_p], marker = 'o', ms = 20, markeredgecolor='black', alpha = .8)
               # ---- uncomment for manuscript ---- 
                # ax_QR.vlines(Kds[QR==np.max(QR)], 1e-9, (Q0*L_0)[QR==np.max(QR)], color = 'black', linewidth = 2, linestyle = '--', alpha = 1)

                ax_QR_all_f.vlines(Kds[QR==np.max(QR)], 1e-9, (Q0*L_0)[QR==np.max(QR)], color = 'black', linewidth = 2, linestyle = '--', alpha = 1)
                
                

            # ---- uncomment for manuscript ----
            if i_t==3 or i_t==4 : 

                ax_QR.plot(Kds, QR*L_0, color = colors_R[i_p][i_t], linewidth = 5, linestyle = '-', alpha = .6)

                ax_QR_all_f.plot(Kds, QR*L_0, color = colors_p[i_p], linewidth = 4, linestyle = '-', alpha = .6)
            
            # ---- uncomment for presentation ----        

            # if i_t == 1:
            #     ax_R.plot(Kds, R, color = colors_R[i_p][i_t], linewidth = 5, linestyle = '-', alpha = alpha_p[i_p])
            #     ax_QR.plot(Kds, QR*L_0, color = colors_R[i_p][i_t], linewidth = 5, linestyle = '-', alpha = alpha_p[i_p])
            #     ax_QR_all_f.plot(Kds, QR*L_0, color = colors_R[i_p][i_t], linewidth = 4, linestyle = '-', alpha = alpha_p[i_p])
            # if i_t == 4:
            #     ax_R.plot(Kds, R, color = colors_R[i_p][i_t], linewidth = 5, linestyle = '-', alpha = alpha_p[i_p])
            #     ax_QR.plot(Kds, QR*L_0, color = colors_R[i_p][i_t], linewidth = 5, linestyle = '-', alpha = alpha_p[i_p])
            #     ax_QR_all_f.plot(Kds, QR*L_0, color = colors_p[i_p], linewidth = 4, linestyle = '-', alpha = alpha_p[i_p])
            # else:
            #     #--------------------------R(E, t) and QR(E, t)---------------------------
            #     #ax_R.plot(Kds, R, alpha = transparency_n[0], color = colors_R[i_p][i_t], linewidth = 4, linestyle = '--')
            #     ax_QR2.plot(Kds, QR*L_0, alpha = transparency_n[0], color = colors_R[i_p][i_t], linewidth = 4, linestyle = '--')
            #     #-------FOR Q0--------- 
            #     #ax_QR.vlines(Kd_r, 0, .5, color = 'black', linestyle = 'dashed')
            #     #ax_QR.vlines([Kd_p, Kd_1], ax_QR.get_ylim()[0], L_0*Q0[Kds<np.exp(E_n)][-1], color = 'grey', linestyle = 'dotted', linewidth = 4)       
            #     #ax_Q_act.hlines(L_0*Q0[Ks<np.exp(E_r)][-1], ax_Q_act.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')

        
        ax_QR2.plot(Kds, QR*L_0, alpha = transparency_n[0], color = colors_R[i_p][i_t], linewidth = 5, linestyle = '-')
        ax_QR2.text(10**-1, ax_QR2.get_ylim()[1]/2, '%.1f'%(np.sum(QR*dE)*L_0), color = colors_R[i_p][i_t])
        ax_QR2.text(10**1, ax_QR2.get_ylim()[1]/2, '%.1e'%(np.sum(QR*dE)), color = colors_R[i_p][i_t])
        ax_QR2.text(10**-1, ax_QR2.get_ylim()[1]/3, '%.1e'%(np.sum(Q0[Kds<10**-5]*dE[Kds<10**-5])), color = colors_R[i_p][i_t])

        my_plot_layout(ax=ax_R, yscale = 'log', xscale = 'log', ticks_labelsize = 30, x_fontsize=30, y_fontsize=30 )
        #ax_R.set_xticks([])
        ax_R.set_xlim(right = 1e-2)
        ax_R.set_ylim(top = 1.5, bottom = 1e-8)
        fig_R.savefig(output_plot + '/primary_response/_Summary/affinity/L%d/R_p-%.1f_Nr-%.0e_'%(L, p, L_0)+model+'.pdf')
        plt.close(fig_R)

        my_plot_layout(ax=ax_QR, yscale = 'log', xscale = 'log', ticks_labelsize = 30, x_fontsize=30, y_fontsize=30 )
        #ax_QR.set_xticks([])
        ax_QR.set_xlim(right = 9e-2, left = np.exp(E_ms)) #use 1e-3 for other plots
        ax_QR.set_ylim(bottom = 1e-8, top = L_0)
        fig_QR.savefig(output_plot + '/primary_response/_Summary/affinity/L%d/QR_p-%.1f_Nr-%.0e_'%(L, p, L_0)+model+'.pdf')
        plt.close(fig_QR)

        my_plot_layout(ax=ax_QR2, yscale = 'linear', xscale = 'log', ticks_labelsize = 30, x_fontsize=30, y_fontsize=30 )
        #ax_QR2.set_xlim(right = 1e-3)
        fig_QR2.savefig(output_plot + '/primary_response/_Summary/affinity/L%d/QR2_p-%.1f_Nr-%.0e_'%(L, p, L_0)+model+'.pdf')
        plt.close(fig_QR2)

        my_plot_layout(ax=ax_m_bar, yscale = 'log', ticks_labelsize = 30, x_fontsize=30, y_fontsize=30 )
        fig_m_bar.savefig(output_plot + '/primary_response/_Summary/affinity/L%d/m_bar_p-%.1f_Nr-%.0e_'%(L, p, L_0)+model+'.pdf')
        plt.close(fig_m_bar)
    
        # my_plot_layout(ax=ax_QR_all, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
        # #ax_QR_all.set_xticks([])
        # ax_QR_all.set_xlim(right = 1e-0, left = 1e-11) #use 1e-3 for other plots
        # ax_QR_all.set_ylim(bottom = 1e-9, top = 1.5*L_0)
        # ax_QR_all.legend(title = r'$p$', title_fontsize = 32, fontsize = 30, loc = 4)
        # fig_QR_all.savefig('../../Figures/_Summary/affinity/L%d/QR_all_p-%.1f_Nr-%.0e_'%(L, p, L_0)+model+'.pdf')

    my_plot_layout(ax=ax_QR_all_f, yscale = 'log', xscale = 'log', ticks_labelsize = 30, x_fontsize=30, y_fontsize=30 )
    ax_QR_all_f.set_yticks([1e-8, 1e-4, 1e-0, 1e4, 1e8])
    ax_QR_all_f.set_xlim(right = 1e-1, left = np.exp(E_ms)) #use 1e-3 for other plots
    ax_QR_all_f.set_ylim(bottom = 1e-9, top = L_0)
    #ax_QR_all_f.legend(title = r'$p$', title_fontsize = 32, fontsize = 30)
    fig_QR_all_f.savefig(output_plot + '/primary_response/_Summary/affinity/L%d/QR_all_f_Nr-%.0e_'%(L, L_0)+model+'.pdf')
    plt.close(fig_QR_all_f)



