import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
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

C = 1e4
AA = 1


ps = [1.4, 1.8, 2.2]
ps = [1, 2, 3, 4]

#ps = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4]

# ps = [1.0, 4]
# E_lims = [-8.5, -12]
# t_lims = [4, 7.2]
# color_list = np.array([my_blue2, my_red])

ps = [1, 1.5, 2.0, 2.5, 4.0]
E_lims = [-8.5, -8.5, -12, -12, -12]
t_lims = [4.0, 4.5, 5.2, 5.7, 7.2]
color_list = np.array([my_blue, my_green, my_purple2, my_brown, my_red])#

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
print('beta_step = %.2f'%beta_pr, 'K_step = %.2e'%Kd_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_ranking, ax_ranking = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_ranking_log, ax_ranking_log = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

for i_p, p in enumerate((ps)):
    E_lim = E_lims[i_p]
    t_lim = t_lims[i_p]
    fig_ranking_i, ax_ranking_i = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    print('--------')
    print('p = %.2f...'%p)
    beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
    beta_act = np.min([beta_r, beta_p])

    #-----------------Loading data----------------------------
    #parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
    n_first_clones = 50
    return_data_type = 0
    data, return_data_type = get_data_b(folder_path = '../../out/primary_immune_response', data_type = 'ranking_combined_2_p-%.1f'%p)

    if(return_data_type):
        final_Nb = data[0]
        final_E = data[1]
        final_E_log = data[2]
        counts_final_Nb = data[3]
        trajectories = data[4]
        trajectories2 = data[5]
        trajectories_rank = data[6]
    else:
        data = pd.read_csv(Text_files_path + 'primary_immune_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
        final_Nb = np.zeros(n_first_clones)
        counts_final_Nb = np.zeros(n_first_clones)
        final_E = np.zeros(n_first_clones)
        final_E_log = np.zeros(n_first_clones)
        max_rank = 50

        trajectories = np.array([], dtype = object)
        trajectories2 = np.array([], dtype = object)
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
            energies_C_sorted = energies_C[sort_inds[:int(n_first_clones)]]
            best_clone_i = energies_C_sorted[0]
            sorted_clones2 = np.exp(-(energies_C_sorted))/np.exp(-best_clone_i)
            clone_sizes_C_sorted = clone_sizes_C[sort_inds[:int(n_first_clones)], :]
            #biggest_clone_i = clone_sizes_C_sorted[-1, -1]
            biggest_clone_i = np.max(clone_sizes[:, -1])
            sorted_clones = (clone_sizes_C_sorted[:, -1])/biggest_clone_i
            max_rank_i = len(sorted_clones)
            if max_rank_i>1:
                for i in range(max_rank_i):
                    final_Nb[i]+= (sorted_clones[i])                    
                    final_E[i]+= (sorted_clones2[i])
                    final_E_log[i]+= np.log(sorted_clones2[i])
                    counts_final_Nb[i] += 1
                if(max_rank_i<max_rank):
                    max_rank = max_rank_i
                if((i_ens%10==0) and (p==3.0)):
                    trajectories = np.append(trajectories, sorted_clones)
                    trajectories2 = np.append(trajectories2, sorted_clones2)
                    trajectories_rank = np.append(trajectories_rank, max_rank_i)

        f = open('../../out/primary_immune_response' + '/processed_data_ranking_combined_2_p-%.1f.pkl'%p, 'wb')
        pickle.dump([final_Nb, final_E, final_E_log, counts_final_Nb, trajectories, trajectories2, trajectories_rank], f, pickle.HIGHEST_PROTOCOL) 

    counter = 0
    for j in range(len(trajectories_rank)):
        ranks_j = np.arange(1, trajectories_rank[j]+1)
        len_rank_j = len(ranks_j)
        #ax_ranking.plot(trajectories[counter:counter+len_rank_j], trajectories2[counter:counter+len_rank_j], color = colors_p[i_p], linewidth = 1, alpha = .2)
        counter += len_rank_j

    final_Nb = final_Nb/counts_final_Nb
    final_E = (final_E/counts_final_Nb)
    final_E_log = np.exp(final_E_log/counts_final_Nb)

    ranking = np.arange(1, n_first_clones+1)

    fit_N = ranking**(-p*lambda_B/(lambda_A*beta_act))
    fit_K = ranking**(1/(beta_act))

    e_array = np.logspace(np.log10(np.min(final_E_log)), np.log10(np.max(final_E_log)))
    #e_array = np.logspace(0, 0.8, 10) 
    x_array = np.logspace(-1.4, 0)
    fit = e_array**(p*lambda_B/(lambda_A))
    #fit/=e_array[0]**(-p*lambda_B/(lambda_A))
    #fit*=final_Nb[-1]
    fit/=e_array[-1]**(p*lambda_B/(lambda_A))
    fit*=final_Nb[0]

    #fit = x_array**(-lambda_A/(p*lambda_B))

    if(p>beta_r):
        ax_ranking.plot(final_E, final_Nb, color = colors_p[i_p], linewidth = 0, marker = '*', alpha = 1, ms = 12, label = r'$%.1f$'%(p))
        ax_ranking.plot(e_array, fit, color = colors_p[i_p], linewidth = 5, alpha = .8)

        ax_ranking_log.plot(final_E_log, final_Nb, color = colors_p[i_p], linewidth = 0, marker = '*', alpha = 1, ms = 12, label = r'$%.1f$'%(p))
        ax_ranking_log.plot(e_array, fit, color = colors_p[i_p], linewidth = 5, alpha = .8)
    else:
        ax_ranking.plot(final_E, final_Nb, color = colors_p[i_p], linewidth = 0, marker = 'o', alpha = .8, ms = 12, label = r'$%.1f$'%(p))

        ax_ranking_log.plot(final_E_log, final_Nb, color = colors_p[i_p], linewidth = 0, marker = 'o', alpha = .8, ms = 12, label = r'$%.1f$'%(p))
        #ax_ranking_log.plot(e_array, fit, color = colors_p[i_p], linewidth = 2, alpha = .8, ls = '--')

my_plot_layout(ax = ax_ranking, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_ranking.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
#ax_ranking.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
#ax_ranking.set_ylim(bottom = 2e-2)
#ax_ranking.set_yticks([1, 0.1, 0.01, 0.001])
#ax_ranking.set_yticklabels([1, 0.1, 0.01])
fig_ranking.savefig('../../../Figures/primary_immune_response/1_Dynamics/CSV/Ranking_combined_2_'+energy_model+'.pdf')

my_plot_layout(ax = ax_ranking_log, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_ranking_log.legend(fontsize = 20, title_fontsize = 22, title = r'$p$')
#ax_ranking_log.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
#ax_ranking_log.set_ylim(bottom = 2e-2)
#ax_ranking_log.set_yticks([1, 0.1, 0.01, 0.001])
#ax_ranking_log.set_yticklabels([1, 0.1, 0.01])
fig_ranking_log.savefig('../../../Figures/primary_immune_response/1_Dynamics/CSV/Ranking_combined_2_log_'+energy_model+'.pdf')


print('----END-----')




