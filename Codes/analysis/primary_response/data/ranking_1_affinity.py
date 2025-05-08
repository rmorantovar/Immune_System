import sys
sys.path.append('../../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------
N_ens = 100
L_0 = 1e8
T0 = 0
Tf = 12
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_As = [5.4]#, 7.5, 9]
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*10
#N_c = 1e5
#E_m = -27.63
E_m = -24
C = 1e4
AA = 1

# ps = [1, 1.5, 2.0, 2.5, 4.0]
# E_lims = [-8.5, -8.5, -12, -12, -12]
# t_lims = [4.0, 4.5, 5.2, 5.7, 7.2]
# color_list = np.array([my_blue, my_green, my_purple2, my_brown, my_red])#

# ps = [1, 2, 3, 4, 5]
# E_lims = [-8.5, -12, -12, -12, -12]
# t_lims = [4.0, 5.2, 6.2, 7.2, 8.2]
# color_list = np.array([my_blue, my_purple2, my_blue2, my_red, my_purple])#

ps = [1.0, 4]
E_lims = [-8.5, -12]
t_lims = [4, 7.2]
color_list = np.array([my_blue, my_red])#

ps = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4, 4.5, 5]
E_lims = [-8.5, -8.5, -12, -12, -12, -12, -12, -12, -12]
t_lims = np.array([4.0, 4.5, 5.2, 5.7, 6.2, 6.7, 7.2, 7.7, 8.2])# + 1
color_list = np.array([my_blue, my_green, my_purple2, my_brown, my_blue2, my_yellow, my_red, my_yellow, my_purple, my_green2])#


transparency_n = [1]

#color_list = np.array([my_blue2, my_red])
#color_list = np.array([my_blue2, my_blue2, my_green, my_green, my_red, my_red, my_gold])

colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18'

L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path + "primary_immune_response/in/")
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_m, L)
Kds = np.exp(Es[:-1])

#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_pr = %.2f'%beta_pr)

# t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------

for lambda_A in lambda_As:

    fig_ranking, ax_ranking = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    fig_ranking_log, ax_ranking_log = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

    for i_p, p in enumerate((ps)):
        E_lim = E_lims[i_p]
        t_lim = t_lims[i_p]
        # fig_ranking_i, ax_ranking_i = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
        print('--------')
        print('p = %.2f...'%p)
        beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
        beta_act = np.min([beta_r, beta_p])

        #-----------------Loading data----------------------------
        #parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
        #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
        #data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
        return_data_type = 0
        data, return_data_type = get_data_b(folder_path = '../../../out/primary_immune_response', data_type = 'ranking_1_affinity_p-%.1f_lamA-%.1f_L0-%.0e'%(p, lambda_A, L_0))
        n_first_clones = 100

        if(return_data_type):
            final_E = data[0]
            final_E_log = data[1]
            counts_final_E = data[2]
            trajectories = data[3]
            trajectories_rank = data[4]
        else:
            data = pd.read_csv(Text_files_path + 'primary_immune_response/simulations_lambda_A/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_lamA_%.1f/activated_population_'%(N_ens, L_0, p, k_step, E_lim, t_lim, lambda_A) + antigen + '.csv')
            final_E = np.zeros(n_first_clones)
            final_E_log = np.zeros(n_first_clones)
            counts_final_E = np.zeros(n_first_clones)
            max_rank = 100

            trajectories = np.array([], dtype = object)
            trajectories_rank = np.array([])

            for i_ens in tqdm(np.arange(N_ens)):
                data_active = data.loc[data['ens_id']==i_ens]
                #data_active = data_i.loc[data_i['active']==1]
                t_act_data = np.min(data_active['time'])
                data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
                activation_times = np.array(data_active['time'])
                energies  = np.array(data_active['E'])

                #---------------------------- B cell linages ----------------------
                clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

                #--------------------------t_C filter-------------------------
                lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
                clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
                
                sort_inds = energies_C.argsort()
                #clone_sizes_C_sorted = clone_sizes_C[sort_inds, :][:int(n_first_clones*(4-3)), :]
                #activation_times_C_sorted = activation_times_C[sort_inds][:int(n_first_clones*(4-3))]
                energies_C_sorted = energies_C[sort_inds][:int(n_first_clones*(4-3))]

                best_clone_i = energies_C_sorted[0]
                #ax_ranking.scatter(np.exp(energies_C_sorted[-1]), biggest_clone_i, color = colors_p[i_p], alpha = .25, edgecolor='black', linewidth=1, facecolor = colors_p[i_p])
                #activation_times_total = np.append(activation_times_total, activation_times_C_sorted)
                sorted_clones = np.exp(-energies_C_sorted)/np.exp(-best_clone_i)
                max_rank_i = len(sorted_clones)

                for i in range(max_rank_i):
                    final_E[i]+= sorted_clones[i]
                    final_E_log[i]+= np.log(sorted_clones[i])
                    counts_final_E[i] += 1
                if(max_rank_i<max_rank):
                    max_rank = max_rank_i
                if((i_ens%10==0) and (p==4.0)):
                    trajectories = np.append(trajectories, sorted_clones)
                    trajectories_rank = np.append(trajectories_rank, max_rank_i)

            f = open('../../../out/primary_immune_response' + '/processed_data_ranking_1_affinity_p-%.1f_lamA-%.1f_L0-%.0e.pkl'%(p, lambda_A, L_0), 'wb')
            pickle.dump([final_E, final_E_log, counts_final_E, trajectories, trajectories_rank], f, pickle.HIGHEST_PROTOCOL)  

        final_E = (final_E/counts_final_E)
        final_E_log = np.exp(final_E_log/counts_final_E)

        counter = 0
        for j in range(len(trajectories_rank)):
            ranks_j = np.arange(1, trajectories_rank[j]+1)
            len_rank_j = len(ranks_j)
            ax_ranking.plot(ranks_j, trajectories[counter:counter+len_rank_j], color = colors_p[i_p], linewidth = 1, alpha = .2)
            ax_ranking_log.plot(ranks_j, trajectories[counter:counter+len_rank_j], color = colors_p[i_p], linewidth = 1, alpha = .2)
            # ax_ranking_i.plot(ranks_j, trajectories[counter:counter+len_rank_j], color = colors_p[i_p], linewidth = 1, alpha = .2)
            counter += len_rank_j


        ranking = np.arange(1, n_first_clones+1)
        fit = ranking**(-1/(beta_act))
        fit*=final_E[0]
        
        ax_ranking.plot(ranking, final_E[:n_first_clones], color = colors_p[i_p], linewidth = 0, marker = '.', alpha = 1, ms = 12, label = r'$%.1f$'%(p))
        #ax_ranking.plot(ranking, fit, color = colors_p[i_p], linewidth = 5, label = r'$%.1f$'%(p), alpha = .8)
        if p == 4:
            ax_ranking_log.plot(ranking, final_E_log[:n_first_clones], color = colors_p[i_p], linewidth = 0, marker = '.', alpha = .8, ms = 12, label = r'$%.1f$'%(p))
            ax_ranking_log.plot(ranking, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)
        else:
            ax_ranking_log.plot(ranking, final_E_log[:n_first_clones], color = colors_p[i_p], linewidth = 0, marker = '.', alpha = .9, ms = 12, label = r'$%.1f$'%(p))

        # ax_ranking_i.plot(ranking, final_E[:n_first_clones], color = colors_p[i_p], linewidth = 0, marker = '.', alpha = 1, ms = 12)
        #ax_ranking_i.plot(ranking, fit, color = colors_p[i_p], linewidth = 5, label = r'$%.1f$'%(p), alpha = .8)

        # my_plot_layout(ax = ax_ranking_i, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
        # ax_ranking_i.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
        #ax_ranking_i.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
        # ax_ranking_i.set_ylim(top = 1e9, bottom = 1e4)
        #ax_ranking_i.set_yticks([1, 0.1, 0.01, 0.001])
        #ax_ranking_i.set_yticklabels([1, 0.1, 0.01])
        #fig_ranking_i.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2020/Ranking_1_affinity_p-%.2f'%(p)+'_'+energy_model+'_%.0e.pdf'%(L_0))

    my_plot_layout(ax = ax_ranking_log, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_ranking_log.legend(fontsize = 15, title_fontsize = 16, title = r'$p$')
    #ax_ranking_log.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
    ax_ranking_log.set_ylim(bottom = 1e-2)
    #ax_ranking_log.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_ranking_log.set_yticklabels([1, 0.1, 0.01])
    fig_ranking_log.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2020/Ranking_1_affinity_log_'+energy_model+'_lamA-%.1f_L0-%.0e.pdf'%(lambda_A, L_0))

    my_plot_layout(ax = ax_ranking, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_ranking.legend(fontsize = 15, title_fontsize = 16, title = r'$p$')
    #ax_ranking.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
    #ax_ranking.set_ylim(bottom = 1e-2)
    #ax_ranking.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_ranking.set_yticklabels([1, 0.1, 0.01])
    fig_ranking.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2020/Ranking_1_affinity_'+energy_model+'_lamA-%.1f_L0-%.0e.pdf'%(lambda_A, L_0))
    plt.close()
print('----END-----')




