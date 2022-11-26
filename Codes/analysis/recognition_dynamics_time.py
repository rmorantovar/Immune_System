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
kappas = [1, 2, 2.5, 3, 4]

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

antigen_color = my_yellow

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])

#colors_kappa = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_kappa = np.flip(['tab:blue','tab:green','tab:red'])
colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['TRCen']
energy_model = 'TRCen'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0


#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18
L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
model = 'TCRen'
#energy_model = 'MJ2'
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
print('beta_a = %.2f'%beta_pr, 'K_a = %.2e'%Kd_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
#--------------------------Loops--------------------------

fig_antigen, ax_antigen = plt.subplots(figsize=(10,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
ax_antigen.plot(time, np.exp(lambda_A*time)/(5*1e3), linewidth = 5, color = antigen_color)

my_plot_layout(ax=ax_antigen, yscale = 'log', xscale = 'linear', ticks_labelsize = 38)
ax_antigen.set_xlim(right = Tf-2, left = T0)
ax_antigen.set_xticks([])
#ax_antigen.set_xlim(right = 1e-2, left = 1e-11) #use 1e-3 for other plots
ax_antigen.set_ylim(bottom = 1e3, top = 1e7)
#ax_antigen.legend(title = r'$\kappa$', title_fontsize = 34, fontsize = 32)
fig_antigen.savefig('../../Figures/_Summary/time/L%d/antigen.pdf'%(L))
plt.close(fig_antigen)

print('Loops...')
for N_r in N_rs:
    print('________')
    print('N_r = %.0e'%N_r)
    fig_QR_all, ax_QR_all = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_m_bar, ax_m_bar = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
    print('beta_r = %.1f'%beta_r, 'K_d_r = %.2e'%Kd_r)

    for i_kappa, kappa in enumerate(kappas):
        print('--------')
        print('kappa = %.2f...'%kappa)
        beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
        print('beta_kappa = %.2f...'%kappa, 'K_kappa = %.2e...'%Kd_kappa)

        fig_R, ax_R = plt.subplots(figsize=(10,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_K, ax_K = plt.subplots(figsize=(10,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_L, ax_L = plt.subplots(figsize=(10,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_QR, ax_QR = plt.subplots(figsize=(10,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        #fig_QR2, ax_QR2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        fig_N_b, ax_N_b = plt.subplots(figsize=(10,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        

        #-----------------Loading data----------------------------
        parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, kappa, N_c, linear, N_ens)+model
        # #data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies%d.txt'%rep, sep = '\t', header=None)
        data = get_data(folder_path = Text_files_path + 'Dynamics/Trajectories/L%d/'%L+parameters_path, rep = 0)
        # #-----------------Filtering data----------------------------
        # min_e_data = np.min(data[0])
        # max_e_data = np.max(data[0])

        data_active = data.loc[data[1]==1]
        # print('Activated clones before filter:%d'%len(data_active[0]))
        t_act_data = np.min(data_active[3])
        print('t_act_data: %.2f'%t_act_data)
        data_active = data_active.loc[data_active[3]<(t_act_data+1.0+0.1*(kappa-1))]
        activation_times = np.array(data_active[3])
        # print('Activated clones after filter:%d'%len(activation_times))
        energies  = np.array(data_active[0])
        # ar1, ar2 = np.histogram(activation_times, bins = time)
        # m_data = np.cumsum(ar1)

        ax_K.plot(time[time>t_prime], (k_pr/k_on)*np.exp(lambda_A*(time[time>t_prime] - t_prime)/kappa), alpha = transparency_n[0], color = colors_kappa[i_kappa], linewidth = 5, linestyle = '-')
        ax_K.hlines(Kd_pr, T0, t_prime, color = colors_kappa[i_kappa], linewidth = 5, linestyle = '-')
        ax_K.hlines(Kd_r, T0, Tf, color = 'grey', linestyle = ':')
        ax_K.hlines(Kd_kappa, T0, Tf, color = 'grey', linestyle = ':')
        ax_K.vlines(t_prime, 1e-14, k_pr/k_on, color = 'grey', linestyle = ':')
        ax_K.vlines(t_prime + kappa/lambda_A*(E_kappa - E_pr), 1e-14, Kd_kappa, color = 'grey', linestyle = ':')
        #--------------------------m_bar(t)---------------------------
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, kappa, lambda_A, N_c, dE)
        M_r = N_r*N_c*np.sum(Q0*p_a*dE)
        
        #m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])# Review this
        m_bar = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time]) 
        m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time))
        t_act = time[m_bar<1][-1]

        ax_m_bar.plot(time, m_bar, linewidth = 5, linestyle = '-', color = colors_kappa[i_kappa])
        ax_m_bar.plot(time, m_bar_approx, linewidth = 1, linestyle = '--', color = 'black')
        ax_m_bar.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
        ax_m_bar.vlines([t_prime, t_prime + kappa/lambda_A*(E_kappa - E_pr)], 1e-4, 1e7, color = 'grey', linestyle = ':')

        ax_L.plot(time, m_bar, linewidth = 5, linestyle = '-', color = colors_kappa[i_kappa])
        ax_L.plot(time, m_bar_approx, linewidth = 1, linestyle = '--', color = 'black')
        ax_L.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
        ax_L.vlines([t_prime + kappa/lambda_A*(E_kappa - E_pr)], 1e-5, 1e6, color = 'grey', linestyle = ':')

        #---------------------------- B cell linages ----------------------
        #clone_sizes = get_clones_sizes_C(int(m_data[-1]), time, activation_times, lambda_B, C, dT)
        #energies = np.array([E_r, E_r+np.log(5), E_r+np.log(10) , E_r+np.log(100), E_r+np.log(200), E_r+np.log(300), E_r+np.log(500), E_r+np.log(1000)])
        #activation_times = t_act + (energies - E_r)*kappa/lambda_A
        #print('Activation times:', activation_times)
        clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)
        lim_size = 2
        clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)

        sort_inds = clone_sizes_C[:, -1].argsort()
        clone_sizes_C_sorted = clone_sizes_C[sort_inds, :][-int(40*(4-3)):, :]
        activation_times_C_sorted = activation_times_C[sort_inds][-int(40*(4-3)):]
        energies_C_sorted = energies_C[sort_inds][-int(40*(4-3)):]

        
        delta_t_n = (E_kappa-E_pr)*kappa/lambda_A
        t_n = t_prime + delta_t_n
        time1 = np.linspace(0, t_n, 100)
        time2 = np.linspace(t_n, Tf, 100)

        print('t_act_theory: %.2f'%t_act)

        ax_N_b.plot(time, clone_sizes_C_sorted[-1, :], linewidth = 5, color = colors_kappa[i_kappa])
        ax_N_b.hlines(C, T0, Tf, color = 'grey', linestyle = ':')
        for k in range(2, len(clone_sizes_C_sorted[:, -1])):
            ax_N_b.plot(time, clone_sizes_C_sorted[-k, :], linewidth = 1.5, color = colors_kappa[i_kappa], linestyle= '-', alpha = .8)
        #ax_N_b.vlines([t_act], 0, C, linestyle = '--', linewidth = 1, color = 'grey')
        
        #ax_N_b.plot(time, np.sum(clone_sizes_C_sorted[:, :], axis = 0) - len(clone_sizes_C_sorted[:, -1]) + 1, linewidth = 5, color = colors_kappa[i_kappa], ls = ':')

        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act+1.3))/N_A, Es, kappa, lambda_A, N_c, dE)
        m_f_expected = np.sum(N_r*QR*dE)
        print('Activated clones expected:%.d'%m_f_expected)

        Kds_plot = [Kd_r*1000, Kd_r*100, Kd_r*10, Kd_r]
        Kds_plot = np.exp(np.array([E_r, E_r+np.log(5), E_r+np.log(10) , E_r+np.log(100), E_r+np.log(200), E_r+np.log(300), E_r+np.log(500), E_r+np.log(1000)]))
        Kds_plot = np.exp(energies_C)[::10]
        #days = np.linspace(1, Tf-0.5, 3)
        for i_Kd, Kd in enumerate(Kds_plot):
            rho_A_t = np.exp(lambda_A*time)/N_A
            u_on, p_a, R_t, QR_t = calculate_QR_t(Q0, k_on, k_pr, np.log(Kd), rho_A_t, Es, kappa, lambda_A, N_c, dE)
            #----------------------------------------------------------------
            #--------------------------QR_all(E, t)---------------------------
            if Kd<Kd_r:
                ax_R.plot(time, R_t, alpha = transparency_n[0], color = colors_kappa[i_kappa], linewidth = 5, linestyle = '-')
                ax_QR.plot(time, QR_t*N_r, alpha = transparency_n[0], color = colors_kappa[i_kappa], linewidth = 5, linestyle = '-')

                ax_QR_all.plot(time, QR_t*N_r, alpha = transparency_n[0], color = colors_kappa[i_kappa], linewidth = 5, linestyle = '-', label = r'$%d$'%(kappa))
                #ax_QR_all.plot(Kds[QR==np.max(QR)], (Q0*N_r)[QR==np.max(QR)], alpha = transparency_n[0], color = colors_R[i_kappa][i_Kd], marker = 'o', ms = 10)

            else:
                #--------------------------R(E, t) and QR(E, t)---------------------------
                ax_R.plot(time, R_t, alpha = .8, color = colors_kappa[i_kappa], linewidth = 2, linestyle = '-')
                ax_QR.plot(time, QR_t*N_r, alpha = transparency_n[0], color = colors_kappa[i_kappa], linewidth = 2, linestyle = '--')
                #-------FOR Q0--------- 
                #ax_QR.vlines(Kd_r, 0, .5, color = 'black', linestyle = 'dashed')
                #ax_QR.vlines([Kd_kappa, Kd_1], ax_QR.get_ylim()[0], N_r*Q0[Kds<np.exp(E_n)][-1], color = 'grey', linestyle = 'dotted', linewidth = 4)       
                #ax_Q_act.hlines(N_r*Q0[Ks<np.exp(E_r)][-1], ax_Q_act.get_xlim()[0], np.exp(E_r), alpha = 1, color = 'black', linestyle = ':')

        my_plot_layout(ax=ax_R, yscale = 'log', xscale = 'linear', ticks_labelsize = 38)
        ax_R.set_xticks([])
        ax_R.set_xlim(right = Tf-2, left = T0)
        ax_R.set_ylim(bottom = 1e-5, top = 1.5)
        fig_R.savefig('../../Figures/_Summary/time/L%d/R_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+energy_model+'.pdf')
        plt.close(fig_R)

        my_plot_layout(ax=ax_K, yscale = 'log', xscale = 'linear', ticks_labelsize = 38)
        ax_K.set_xticks([])
        ax_K.set_xlim(right = Tf-2, left = T0)
        ax_K.set_ylim(bottom = 1e-9, top = 1e-5)
        fig_K.savefig('../../Figures/_Summary/time/L%d/K_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+energy_model+'.pdf')
        plt.close(fig_K)

        my_plot_layout(ax=ax_L, yscale = 'log', xscale = 'linear', ticks_labelsize = 38)
        ax_L.set_xticks([])
        ax_L.set_xlim(right = Tf-2, left = T0)
        ax_L.set_ylim(bottom = 1e-1, top = 1e3)
        fig_L.savefig('../../Figures/_Summary/time/L%d/L_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+energy_model+'.pdf')
        plt.close(fig_L)

        my_plot_layout(ax=ax_QR, yscale = 'log', xscale = 'linear', ticks_labelsize = 38)
        #ax_QR.set_xticks([])
        #ax_QR.set_xlim(right = 1e-2) #use 1e-3 for other plots
        ax_QR.set_ylim(bottom = 1e-8)
        fig_QR.savefig('../../Figures/_Summary/time/L%d/QR_kappa-%.1f_Nr-%.0e_'%(L, kappa, N_r)+energy_model+'.pdf')
        plt.close(fig_QR)

        # my_plot_layout(ax=ax_QR2, yscale = 'linear', xscale = 'log', ticks_labelsize = 30)
        # #ax_QR2.set_xlim(right = 1e-3)
        # fig_QR2.savefig('../../Figures/_Summary/time/QR2_kappa-%.1f_Nr-%.0e_'%(kappa, N_r)+energy_model+'.pdf')
        # plt.close(fig_QR2)

        my_plot_layout(ax=ax_N_b, yscale = 'log', ticks_labelsize = 38)
        ax_N_b.set_xlim(right = Tf-2, left = T0)
        ax_N_b.set_ylim(bottom = 1e0, top = 1.05e4)
        fig_N_b.savefig('../../Figures/_Summary/time/L%d/Bcell_clones-%.1f_Nr-%.0e_'%(L, kappa, N_r)+energy_model+'.pdf')
        plt.close(fig_N_b)

    my_plot_layout(ax=ax_m_bar, yscale = 'log', ticks_labelsize = 30)
    fig_m_bar.savefig('../../Figures/_Summary/time/L%d/m_bar_Nr-%.0e_'%(L, N_r)+energy_model+'.pdf')
    plt.close(fig_m_bar)
    
    my_plot_layout(ax=ax_QR_all, yscale = 'log', xscale = 'linear', ticks_labelsize = 38)
    ax_QR_all.set_xticks([4, 6, 8])
    #ax_QR_all.set_xlim(right = 1e-2, left = 1e-11) #use 1e-3 for other plots
    ax_QR_all.set_ylim(bottom = 1e-9)
    ax_QR_all.legend(title = r'$\kappa$', title_fontsize = 34, fontsize = 32)
    fig_QR_all.savefig('../../Figures/_Summary/time/L%d/QR_all_Nr-%.0e_'%(L, N_r)+energy_model+'.pdf')
    plt.close(fig_QR_all)

    # my_plot_layout(ax=ax_Q0, yscale = 'log', xscale = 'log', ticks_labelsize = 38)
    # #ax_Q_act.set_xticks([])
    # ax_Q0.set_xlim(right = 1e1) #use 1e-3 for other plots
    # ax_Q0.set_ylim(bottom = 1e-9, top = 1.5)
    # fig_Q0.savefig('../../Figures/_Summary/time/Q0_Nr-%.0e_'%(N_r)+energy_model+'.pdf')
    # plt.close(fig_Q0)




