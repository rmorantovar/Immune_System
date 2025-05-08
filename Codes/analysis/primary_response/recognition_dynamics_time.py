import sys
sys.path.append('../lib/')
sys.path.append('../../lib/')
# from functions_2 import*
from funcs import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/primary_response/'

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
N_c = 1e5*10
#N_c = 1e5
#E_ms = -27.63
E_ms = -24
C = 1e4
AA = 1

ps = [1, 2, 2.5, 3, 4]
ps = [1, 3]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])
color_list = np.array([my_blue2, my_red, my_green])

#colors_p = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_p = np.flip(['tab:blue','tab:green','tab:red'])
colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['TRCen']
energy_model = 'TCRen'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0
alpha_p = [.6, 1]

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
motif = get_motif(from_aa_to_i(antigen, energy_model, Text_files_path), model, Text_files_path)
print('min_e_PWM=%.4f'%(np.sum([np.min(motif[:,i]) for i in range(len(motif[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(motif[:,i]) for i in range(len(motif[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    motif[:,i]-=np.min(motif[:,i], axis=0)

avg_E = np.sum([np.mean(motif[:,i]) for i in range(len(motif[0,:]))]) + E_ms
var_E = np.sum([np.var(motif[:,i]) for i in range(len(motif[0,:]))])

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_ms, L)
Kds = np.exp(Es[:-1])
beta_1, E_1, Kd_1 = get_p_properties(betas, Q0, Es, dE, 1)
#--------------------------Proofreading properties--------------------------
beta_step, E_step, K_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_a = %.2f'%beta_step, 'K_a = %.2e'%K_step)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
#--------------------------Loops--------------------------

fig_antigen, ax_antigen = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
ax_antigen.plot(time_array, np.exp(lambda_A*time_array)/(1e0), linewidth = 5, color = antigen_color)

my_plot_layout(ax=ax_antigen, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
ax_antigen.set_xlim(right = Tf-1, left = T0+1)
ax_antigen.set_xticks([])
#ax_antigen.set_xlim(right = 1e-2, left = 1e-11) #use 1e-3 for other plots
ax_antigen.set_ylim(bottom = 2e3, top = 2e12)
#ax_antigen.set_ylim(bottom = 1, top = 1e7)
#ax_antigen.legend(title = r'$\p$', title_fontsize = 34, fontsize = 32)
fig_antigen.savefig('../../../Figures/primary_response/_Summary/time/L%d/antigen.pdf'%(L))
plt.close(fig_antigen)

print('Loops...')
for L_0 in L_0s:
    print('________')
    print('L_0 = %.0e'%L_0)
    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
    print('beta_r = %.1f'%beta_r, 'K_d_r = %.2e'%Kd_r)

    fig_R, ax_R = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_K, ax_K = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_L, ax_L = plt.subplots(figsize=(8*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_N_b, ax_N_b = plt.subplots(figsize=(8.0*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

    for i_p, p in enumerate(ps):

        fig_N_b_p, ax_N_b_p = plt.subplots(figsize=(8.0*1.62,8*0.5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.2, 'top': 0.96})
        # ax_N_b_p.plot(time_array, np.exp(3*time_array)/(1e0), linewidth = 5, color = antigen_color)
        print('--------')
        print('p = %.2f...'%p)
        beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
        print('beta_p = %.2f...'%p, 'K_p = %.2e...'%Kd_p)
        
        #-----------------Loading data----------------------------
        parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+model
        # #data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies%d.txt'%rep, sep = '\t', header=None)
        #data = get_data(folder_path = Text_files_path + 'Dynamics/Trajectories/L%d/'%L+parameters_path, rep = 0)
        # #-----------------Filtering data----------------------------
        # min_e_data = np.min(data[0])
        # max_e_data = np.max(data[0])
        #--------------------------m_bar(t)---------------------------
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*time_array[0])/N_A, Es, p, lambda_A, N_c, dE)
        M_r = L_0*N_c*np.sum(Q0*p_a*dE)
        
        #m_bar = np.array([L_0*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time_array])# Review this
        m_bar = np.array([np.sum(L_0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array]) 
        m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time_array))
        t_act = time_array[m_bar<1][-1]
        print('t_act:', t_act)

        ax_L.plot(time_array, m_bar, linewidth = 5, linestyle = '-', color = colors_p[i_p])
        #ax_L.plot(time_array, m_bar_approx, linewidth = 1, linestyle = '--', color = 'black')
        #ax_L.hlines(1, T0, Tf, color = 'grey', linestyle = ':')
        #ax_L.vlines([t_prime + p/lambda_A*(E_p - E_step)], 1e-5, 1e6, color = 'grey', linestyle = ':')

        #---------------------------- B cell linages ----------------------
        t_f = Tf
        n_t = int(m_bar[time_array<t_f][-1])
        print('L_act at Tf = ', n_t)
        activation_times = np.array([time_array[m_bar<=n][-1] for n in range(1, 5000, 5)])
        #start_time = time.time()
        clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
        #print("--- %s seconds ---" % (time.time() - start_time))
        #--------------------------t_C filter-------------------------
        lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
        print('lower N_b limit =', lim_size)
        clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, np.ones_like(activation_times), lim_size)
        print('final L_act =', np.size(activation_times_C))

        delta_t_n = (E_p-E_step)*p/lambda_A
        t_n = t_prime + delta_t_n
        time1 = np.linspace(0, t_n, 100)
        time2 = np.linspace(t_n, Tf, 100)

        print('t_act_theory: %.2f'%t_act)

        #ax_N_b.plot(time_array, clone_sizes_C_sorted[-1, :], linewidth = 5, color = colors_p[i_p])
        #ax_N_b.plot(time_array, clone_sizes_C[0, :]-np.heaviside(activation_times_C[0] - time_array , 1), linewidth = 5, color = colors_p[i_p])

        bar_N_1 = clone_sizes_C[0, -1]

        ax_N_b.hlines(C, T0, Tf, color = 'grey', linestyle = ':')
        # for k in range(2, len(clone_sizes_C_sorted[:, -1])):
        #     ax_N_b.plot(time_array, clone_sizes_C_sorted[-k, :], linewidth = 1.5, color = colors_p[i_p], linestyle= '-', alpha = .8)
        my_colors = [my_green_a, my_green_b, my_green_c]
        for i_k, k in enumerate([0, 1, 10]):
            ax_N_b.plot(time_array, clone_sizes_C[k*10, :]-np.heaviside(activation_times_C[k*10] - time_array , 1), linewidth = 1.5, color = colors_p[i_p], linestyle= '-', alpha = .8)
            ax_N_b_p.plot(time_array, clone_sizes_C[k*10, :]-np.heaviside(activation_times_C[k*10] - time_array , 1), linewidth = 3, color = my_colors[i_k], linestyle= '-', alpha = .8)

        for f in [.99, .1, .01]:
            t_f = activation_times_C[clone_sizes_C[:, -1]>(bar_N_1*f)][-1]
            print('f:', f, 't_f:',  t_f)
            #ax_N_b.hlines(bar_N_1*f, T0, Tf, color = 'grey', linestyle = '--')
            #ax_N_b.vlines(t_f, 1, bar_N_1*f, color = 'grey', linestyle = '--')
        
        #ax_N_b.vlines([t_act], 0, C, linestyle = '--', linewidth = 1, color = 'grey')
        #ax_N_b.plot(time_array, np.sum(clone_sizes_C_sorted[:, :], axis = 0) - len(clone_sizes_C_sorted[:, -1]) + 1, linewidth = 5, color = colors_p[i_p], ls = ':')

        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t_act+1.3))/N_A, Es, p, lambda_A, N_c, dE)
        m_f_expected = np.sum(L_0*QR*dE)
        print('Activated clones expected:%.d'%m_f_expected)

        
        K_front = (k_step/k_on)*np.exp(lambda_A*(time_array[time_array>t_prime] - t_prime)/p)
        t_p = time_array[time_array>t_prime][K_front>Kd_p][0]
        print('t_p:', t_p)
        # ax_K.hlines(K_step, T0, t_prime, color = colors_p[i_p], linewidth = 5, linestyle = '--', alpha = .5)
        ax_K.plot(time_array[time_array>t_prime], K_front, alpha = .5, color = colors_p[i_p], linewidth = 5, linestyle = '--')
        ax_K.plot(time_array[time_array>1.15][0], Kd_r, markerfacecolor = 'grey', marker = '*', ms = 22, markeredgecolor='black', alpha = 0.8, zorder = 20)
        ax_K.plot(time_array[time_array>1.15][0], Kd_p, markerfacecolor = colors_p[i_p], marker = 'o', ms = 18, markeredgecolor='black', alpha = alpha_p[i_p], zorder = 20)

        if(p>beta_r):
            #ax_K.hlines(Kd_r, T0, t_act, color = colors_p[i_p], linewidth = 5, linestyle = '-', alpha = 1)
            ax_K.plot(time_array[time_array>=(t_act+0.2)], (Kd_r)*np.exp(lambda_A*(time_array[time_array>=(t_act+0.2)] - (t_act+0.2))/p), alpha = 1, color = colors_p[i_p], linewidth = 5, linestyle = '-')
        else:
            ax_K.hlines(Kd_p, t_act+0.15, t_p, color = colors_p[i_p], linewidth = 5, linestyle = '-', alpha = 1)
            ax_K.plot(time_array[time_array>=t_p], (Kd_p)*np.exp(lambda_A*(time_array[time_array>=t_p] - t_p)/p), alpha = 1, color = colors_p[i_p], linewidth = 5, linestyle = '-')
        print('t`:', t_prime, 't_act:', t_act)
        #ax_K.hlines(Kd_r, T0, Tf, color = 'grey', linestyle = '--')
        #ax_K.hlines(Kd_p, T0, Tf, color = 'grey', linestyle = ':')
        #ax_K.vlines(t_prime, 1e-14, k_step/k_on, color = 'grey', linestyle = ':')
        #ax_K.vlines(t_prime + p/lambda_A*(E_p - E_step), 1e-14, Kd_p, color = 'grey', linestyle = ':')

        my_plot_layout(ax=ax_N_b_p, yscale = 'log', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
        if p==1:
            ax_N_b_p.set_xticks([])
        ax_N_b_p.set_xlim(right = Tf-1, left = T0+1)
        ax_N_b_p.set_ylim(bottom = 1e0, top = 2e2)
        #ax_N_b.set_ylim(bottom = 1e0, top = C*1.1)
        fig_N_b_p.savefig('../../../Figures/primary_response/_Summary/time/L%d/Bcell_p-%.1f_Nr-%.0e_'%(L, p, L_0)+energy_model+'.pdf')
        plt.close(fig_N_b_p)

    my_plot_layout(ax=ax_N_b, yscale = 'log', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    ax_N_b.set_xlim(right = Tf-1, left = T0+1)
    ax_N_b.set_ylim(bottom = 1e0, top = 5e2)
    #ax_N_b.set_ylim(bottom = 1e0, top = C*1.1)
    fig_N_b.savefig('../../../Figures/primary_response/_Summary/time/L%d/Bcell_Nr-%.0e_'%(L, L_0)+energy_model+'.pdf')
    plt.close(fig_N_b)

    my_plot_layout(ax=ax_R, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    #ax_R.set_xticks([])
    ax_R.set_xlim(right = Tf-1, left = T0+1)
    ax_R.set_ylim(bottom = 1.5e-3, top = 1.5)
    fig_R.savefig('../../../Figures/primary_response/_Summary/time/L%d/R_Nr-%.0e_'%(L, L_0)+energy_model+'.pdf')
    plt.close(fig_R)

    my_plot_layout(ax=ax_K, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    # ax_K.set_xticks([])
    ax_K.set_xlim(right = Tf-1, left = T0+1)
    ax_K.set_ylim(bottom = 9e-10, top = 2e-4)
    fig_K.savefig('../../../Figures/primary_response/_Summary/time/L%d/K_Nr-%.0e_'%(L, L_0)+energy_model+'.pdf')
    plt.close(fig_K)

    my_plot_layout(ax=ax_L, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    ax_L.set_xticks([])
    ax_L.set_xlim(right = Tf-1, left = T0+1)
    ax_L.set_ylim(bottom = 9e-1, top = 9e3)
    ax_L.set_yticks([1e0, 1e2])
    fig_L.savefig('../../../Figures/primary_response/_Summary/time/L%d/L_Nr-%.0e_'%(L, L_0)+energy_model+'.pdf')
    plt.close(fig_L)





