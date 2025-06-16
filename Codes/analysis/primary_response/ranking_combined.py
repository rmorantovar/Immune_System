import sys
sys.path.append('../library/')
sys.path.append('../../my_lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------
N_ens = 400
L_0 = 1e9
T0 = 0
Tf = 12
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
#E_m = -27.63
E_m = -24
E_lims = [-8.5, -12]
t_lims = [4, 7.2]
C = 1e4
AA = 1

ps = [1.0, 4]
E_lims = [-8.5, -12]
t_lims = [4, 7.2]
color_list = np.array([my_blue, my_red])
color_list = np.array([my_red, my_blue2])

# ps = [1, 1.5, 2.0, 2.5, 4.0]
# E_lims = [-8.5, -8.5, -12, -12, -12]
# t_lims = [4.0, 4.5, 5.2, 5.7, 7.2]
# color_list = np.array([my_blue, my_green, my_purple2, my_brown, my_red])#


transparency_n = [1]


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
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path + "primary_response/in/")
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
print('beta_step = %.2f'%beta_pr, 'K_step = %.2e'%Kd_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_ranking, ax_ranking = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
fig_ranking_log, ax_ranking_log = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

for i_p, p in enumerate((ps)):
    E_lim = E_lims[i_p]
    t_lim = t_lims[i_p]
    fig_ranking_i, ax_ranking_i = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    print('--------')
    print('p = %.2f...'%p)
    beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
    beta_act = np.min([beta_r, beta_p])

    #-----------------Loading data----------------------------
    #parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
    return_data_type = 0
    data, return_data_type = get_data(folder_path = Text_files_path + 'primary_response/out/', data_type = 'ranking_combined_try_p-%.1f'%p)
    n_first_clones = 100

    if(return_data_type):
        final_Nb_size = data[0]
        final_E_size = data[1]
        final_E_log_size = data[2]
        counts_final_Nb_size = data[3]
        trajectories_size = data[4]
        trajectories2_size = data[5]
        trajectories_rank_size = data[6]

        final_Nb_E = data[7]
        final_E_E = data[8]
        final_E_log_E = data[9]
        counts_final_Nb_E = data[10]
        trajectories_E = data[11]
        trajectories2_E = data[12]
        trajectories_rank_E = data[13]

    else:
        data = pd.read_csv(Text_files_path + 'primary_response/simulations/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
        final_Nb_size = np.zeros(n_first_clones)
        counts_final_Nb_size = np.zeros(n_first_clones)
        final_E_size = np.zeros(n_first_clones)
        final_E_log_size = np.zeros(n_first_clones)
        max_rank = 100

        trajectories_size = np.array([], dtype = object)
        trajectories2_size = np.array([], dtype = object)
        trajectories_rank_size = np.array([])

        final_Nb_E = np.zeros(n_first_clones)
        counts_final_Nb_E = np.zeros(n_first_clones)
        final_E_E = np.zeros(n_first_clones)
        final_E_log_E = np.zeros(n_first_clones)
        max_rank = 100

        trajectories_E = np.array([], dtype = object)
        trajectories2_E = np.array([], dtype = object)
        trajectories_rank_E = np.array([])

        fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

        for i_ens in tqdm(np.arange(N_ens)):
            data_active = data.loc[data['ens_id']==i_ens]
            #data_active = data_i.loc[data_i['active']==1]
            t_act_data = np.min(data_active['time'])
            data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
            activation_times = np.array(data_active['time'])
            energies  = np.array(data_active['E'])
            if i_ens%10==0:
                ax.scatter(activation_times, energies)
            #---------------------------- B cell linages ----------------------
            clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

            #--------------------------t_C filter-------------------------
            lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
            clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
            #print(energies_C, clone_sizes_C[:, -1])
            #print(np.min(energies_C), np.max(clone_sizes_C[:, -1]))

            sort_inds_size = clone_sizes_C[:, -1].argsort()
            clone_sizes_C_sorted_size = clone_sizes_C[sort_inds_size[-int(n_first_clones):],:]
            biggest_clone_i_size = np.max(clone_sizes_C[:, -1])
            sorted_clones_size = np.flip(clone_sizes_C_sorted_size[:, -1])/biggest_clone_i_size
            energies_C_sorted_size = energies_C[sort_inds_size[-int(n_first_clones):]]
            #best_clone_i_size = np.min(energies_C_sorted_size)
            best_clone_i_size = np.min(energies_C)
            #best_clone_i_size = np.min(energies)
            sorted_clones2_size = np.exp(-np.flip(energies_C_sorted_size))/np.exp(-best_clone_i_size)
            max_rank_i = len(sorted_clones_size)
            if max_rank_i>1:
                for i in range(max_rank_i):
                    final_Nb_size[i]+= (sorted_clones_size[i])                    
                    final_E_size[i]+= (sorted_clones2_size[i])
                    final_E_log_size[i]+= np.log(sorted_clones2_size[i])
                    counts_final_Nb_size[i] += 1
                if(max_rank_i<max_rank):
                    max_rank = max_rank_i
                if((i_ens%10==0) and (p==4.0)):
                    trajectories_size = np.append(trajectories_size, sorted_clones_size)
                    trajectories2_size = np.append(trajectories2_size, sorted_clones2_size)
                    trajectories_rank_size = np.append(trajectories_rank_size, max_rank_i)


            sort_inds_E = energies_C.argsort()
            energies_C_sorted_E = energies_C[sort_inds_E[:int(n_first_clones)]]
            best_clone_i_E = np.min(energies_C)
            sorted_clones2_E = np.exp(-(energies_C_sorted_E))/np.exp(-best_clone_i_E)
            clone_sizes_C_sorted_E = clone_sizes_C[sort_inds_E[:int(n_first_clones)], :]
            #biggest_clone_i_E = np.max(clone_sizes_C_sorted_E[:, -1])
            biggest_clone_i_E = np.max(clone_sizes_C[:, -1])
            sorted_clones_E = (clone_sizes_C_sorted_E[:, -1])/biggest_clone_i_E
            max_rank_i = len(sorted_clones_E)
            if max_rank_i>1:
                for i in range(max_rank_i):
                    final_Nb_E[i]+= (sorted_clones_E[i])                    
                    final_E_E[i]+= (sorted_clones2_E[i])
                    final_E_log_E[i]+= np.log(sorted_clones2_E[i])
                    counts_final_Nb_E[i] += 1
                if(max_rank_i<max_rank):
                    max_rank = max_rank_i
                if((i_ens%10==0) and (p==4.0)):
                    trajectories_E = np.append(trajectories_E, sorted_clones_E)
                    trajectories2_E = np.append(trajectories2_E, sorted_clones2_E)
                    trajectories_rank_E = np.append(trajectories_rank_E, max_rank_i)


            #print(energies_C_sorted_size, clone_sizes_C_sorted_size[:, -1])
            #print(biggest_clone_i_size, best_clone_i_size)
            #print(energies_C_sorted_E, clone_sizes_C_sorted_E[:, -1])
            #print(biggest_clone_i_E, best_clone_i_E)

        my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
        #ax.legend(fontsize = 20, title_fontsize = 22, title = r'$p$')
        #ax.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
        #ax.set_ylim(bottom = 2e-2)
        #ax.set_yticks([1, 0.1, 0.01, 0.001])
        #ax.set_yticklabels([1, 0.1, 0.01])
        fig.savefig('/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/primary_response/1_Dynamics/CSV/hist_p-%.1f'%p+energy_model+'.pdf')


        f = open(Text_files_path + 'primary_response/out' + '/processed_data_ranking_combined_try_p-%.1f.pkl'%p, 'wb')
        pickle.dump([final_Nb_size, final_E_size, final_E_log_size, counts_final_Nb_size, trajectories_size, trajectories2_size, trajectories_rank_size, final_Nb_E, final_E_E, final_E_log_E, counts_final_Nb_E, trajectories_E, trajectories2_E, trajectories_rank_E], f, pickle.HIGHEST_PROTOCOL) 

    counter = 0
    for j in range(len(trajectories_rank_size)):
        ranks_j = np.arange(1, trajectories_rank_size[j]+1)
        len_rank_j = len(ranks_j)
        #ax_ranking.plot(trajectories[counter:counter+len_rank_j], trajectories2[counter:counter+len_rank_j], color = colors_p[i_p], linewidth = 1, alpha = .2)
        counter += len_rank_j

    final_Nb_size = final_Nb_size /counts_final_Nb_size 
    final_E_size  = (final_E_size/counts_final_Nb_size )
    final_E_log_size  = np.exp(final_E_log_size /counts_final_Nb_size )

    ranking = np.arange(1, n_first_clones+1)

    fit_N = ranking**(-p*lambda_B/(lambda_A*beta_act))
    fit_K = ranking**(1/(beta_act))

    e_array = np.logspace(np.log10(np.min(final_E_log_size)), np.log10(np.max(final_E_log_size)))
    #e_array = np.logspace(0, 0.8, 10) 
    x_array = np.logspace(-1.4, 0)
    fit = e_array**(p*lambda_B/(lambda_A))
    #fit/=e_array[0]**(-p*lambda_B/(lambda_A))
    #fit*=final_Nb[-1]
    fit/=e_array[-1]**(p*lambda_B/(lambda_A))
    fit*=final_Nb_size[0]

    #fit = x_array**(-lambda_A/(p*lambda_B))

    if(p>beta_r):
        ax_ranking.plot(final_E_size, final_Nb_size, color = colors_p[i_p], linewidth = 0, marker = '*', alpha = 1, ms = 12, label = r'$%.1f$'%(p))
        ax_ranking.plot(e_array, fit, color = colors_p[i_p], linewidth = 5, alpha = .8)

        ax_ranking_log.plot(final_E_log_size, final_Nb_size, color = colors_p[i_p], linewidth = 0, marker = '*', alpha = .8, ms = 10, label = r'$%.1f$'%(p))
        ax_ranking_log.plot(e_array, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)
    else:
        ax_ranking.plot(final_E_size, final_Nb_size, color = colors_p[i_p], linewidth = 0, marker = '*', alpha = .8, ms = 12, label = r'$%.1f$'%(p))

        ax_ranking_log.plot(final_E_log_size, final_Nb_size, color = colors_p[i_p], linewidth = 0, marker = '*', alpha = .8, ms = 10, label = r'$%.1f$'%(p))
        #ax_ranking_log.plot(e_array, fit, color = colors_p[i_p], linewidth = 2, alpha = .8, ls = '--')


    counter = 0
    for j in range(len(trajectories_rank_E)):
        ranks_j = np.arange(1, trajectories_rank_E[j]+1)
        len_rank_j = len(ranks_j)
        #ax_ranking.plot(trajectories[counter:counter+len_rank_j], trajectories2[counter:counter+len_rank_j], color = colors_p[i_p], linewidth = 1, alpha = .2)
        counter += len_rank_j

    final_Nb_E = final_Nb_E/counts_final_Nb_E
    final_E_E = (final_E_E/counts_final_Nb_E)
    final_E_log_E = np.exp(final_E_log_E/counts_final_Nb_E)

    ranking = np.arange(1, n_first_clones+1)

    fit_N = ranking**(-p*lambda_B/(lambda_A*beta_act))
    fit_K = ranking**(1/(beta_act))

    e_array = np.logspace(np.log10(np.min(final_E_log_E)), np.log10(np.max(final_E_log_E)))
    #e_array = np.logspace(0, 0.8, 10) 
    x_array = np.logspace(-1.4, 0)
    fit = e_array**(p*lambda_B/(lambda_A))
    #fit/=e_array[0]**(-p*lambda_B/(lambda_A))
    #fit*=final_Nb[-1]
    fit/=e_array[-1]**(p*lambda_B/(lambda_A))
    fit*=final_Nb_E[0]

    #fit = x_array**(-lambda_A/(p*lambda_B))

    if(p>beta_r):
        ax_ranking.plot(final_E_E, final_Nb_E, color = colors_p[i_p], linewidth = 0, marker = '.', alpha = 1, ms = 12, label = r'$%.1f$'%(p))
        ax_ranking.plot(e_array, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)
        #ax_ranking_log.plot(final_E_log_E, final_Nb_E, color = colors_p[i_p], linewidth = 0, marker = '.', alpha = 1, ms = 12, label = r'$%.1f$'%(p))
        #ax_ranking_log.plot(e_array, fit, color = colors_p[i_p], linewidth = 5, alpha = .8)
    else:
        ax_ranking.plot(final_E_E, final_Nb_E, color = colors_p[i_p], linewidth = 0, marker = '.', alpha = .8, ms = 12, label = r'$%.1f$'%(p))

        ax_ranking_log.plot(final_E_log_E, final_Nb_E, color = colors_p[i_p], linewidth = 0, marker = '.', alpha = .8, ms = 12, label = r'$%.1f$'%(p))


my_plot_layout(ax = ax_ranking, xscale='log', yscale= 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
#ax_ranking.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
#ax_ranking.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
#ax_ranking.set_ylim(bottom = 2e-2)
#ax_ranking.set_yticks([1, 0.1, 0.01, 0.001])
#ax_ranking.set_yticklabels([1, 0.1, 0.01])
fig_ranking.savefig('/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/primary_response/1_Dynamics/CSV/Ranking_combined_'+energy_model+'.pdf')

my_plot_layout(ax = ax_ranking_log, xscale='log', yscale= 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
#ax_ranking_log.legend(fontsize = 20, title_fontsize = 22, title = r'$p$')
#ax_ranking_log.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
#ax_ranking_log.set_ylim(bottom = 2e-2)
#ax_ranking_log.set_yticks([1, 0.1, 0.01, 0.001])
#ax_ranking_log.set_yticklabels([1, 0.1, 0.01])
fig_ranking_log.savefig('/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/primary_response/1_Dynamics/CSV/Ranking_combined_log_'+energy_model+'.pdf')



print('----END-----')




