import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 500
N_r = 1e8
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
kappas = [1, 2, 3, 4]
kappas = [1, 2, 2.5, 3, 4]
kappas = [1, 3]
#kappas = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_red])
#color_list = np.array([my_blue2, my_blue2, my_green, my_green, my_red, my_red, my_gold])


colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

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
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path)
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_a = %.2f'%beta_pr, 'K_a = %.2e'%Kd_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_ranking, ax_ranking = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

for i_kappa, kappa in enumerate((kappas)):
    fig_ranking_i, ax_ranking_i = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    print('--------')
    print('kappa = %.2f...'%kappa)
    beta_kappa, E_kappa, Kd_kappa = get_p_properties(betas, Q0, Es, dE, kappa)
    beta_act = np.min([beta_r, beta_kappa])

    #-----------------Loading data----------------------------
    parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
    #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
    #data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
    data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'ranking_size')
    n_first_clones = 50

    if(return_data_type):
        final_Nb = data[0]
        counts_final_Nb = data[1]
        trajectories = data[2]
        trajectories_rank = data[3]
    else:
        #activation_times_total = np.array([])
        
        final_Nb = np.zeros(n_first_clones)
        counts_final_Nb = np.zeros(n_first_clones)
        max_rank = 50

        trajectories = np.array([], dtype = object)
        trajectories_rank = np.array([])

        for i_ens in tqdm(np.arange(N_ens)):
            data_active = data.loc[data['i_ens']==i_ens]
            #data_active = data_i.loc[data_i['active']==1]
            t_act_data = np.min(data_active['act_time'])
            data_active = data_active.loc[data_active['act_time']<(t_act_data+1.2+0.3*(kappa-1))] # it was 1.0 + 0.1*...
            activation_times = np.array(data_active['act_time'])
            energies  = np.array(data_active['energy'])

            #---------------------------- B cell linages ----------------------
            clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

            #--------------------------t_C filter-------------------------
            lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
            clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
            
            sort_inds = clone_sizes_C[:, -1].argsort()
            clone_sizes_C_sorted = clone_sizes_C[sort_inds, :][-int(n_first_clones*(4-3)):, :]
            #activation_times_C_sorted = activation_times_C[sort_inds][-int(n_first_clones*(4-3)):]
            #energies_C_sorted = energies_C[sort_inds][-int(n_first_clones*(4-3)):]

            biggest_clone_i = clone_sizes_C_sorted[-1, -1]
            #ax_ranking.scatter(np.exp(energies_C_sorted[-1]), biggest_clone_i, color = colors_kappa[i_kappa], alpha = .25, edgecolor='black', linewidth=1, facecolor = colors_kappa[i_kappa])
            #activation_times_total = np.append(activation_times_total, activation_times_C_sorted)
            sorted_clones = np.flip(clone_sizes_C_sorted[:, -1])/biggest_clone_i
            max_rank_i = len(sorted_clones)
            if max_rank_i>1:
                for i in range(max_rank_i):
                    #final_Nb[i]+= np.log(sorted_clones[i])
                    final_Nb[i]+= (sorted_clones[i])
                    counts_final_Nb[i] += 1
                if(max_rank_i<max_rank):
                    max_rank = max_rank_i
                if((i_ens%10==0) and (kappa==3.0)):
                    trajectories = np.append(trajectories, sorted_clones)
                    trajectories_rank = np.append(trajectories_rank, max_rank_i)

        f = open(Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path+'/processed_data_ranking_size.pkl', 'wb')
        pickle.dump([final_Nb, counts_final_Nb, trajectories, trajectories_rank], f, pickle.HIGHEST_PROTOCOL) 
    
    counter = 0
    for j in range(len(trajectories_rank)):
        ranks_j = np.arange(1, trajectories_rank[j]+1)
        len_rank_j = len(ranks_j)
        ax_ranking.plot(ranks_j, trajectories[counter:counter+len_rank_j], color = colors_kappa[i_kappa], linewidth = 1, alpha = .2)
        ax_ranking_i.plot(ranks_j, trajectories[counter:counter+len_rank_j], color = colors_kappa[i_kappa], linewidth = 1, alpha = .2)
        counter += len_rank_j


    #final_Nb = np.exp(final_Nb/counts_final_Nb)
    final_Nb = final_Nb/counts_final_Nb
    
    ranking = np.arange(1, n_first_clones+1)
    fit = ranking**(-kappa*lambda_B/(lambda_A*beta_act))
    ax_ranking.plot(ranking, final_Nb, color = colors_kappa[i_kappa], linewidth = 0, marker = '*', alpha = 1, ms = 12)
    ax_ranking.plot(ranking, fit, color = colors_kappa[i_kappa], linewidth = 5, label = r'$%.d$'%(kappa), alpha = .8)

    ax_ranking_i.plot(ranking, final_Nb, color = colors_kappa[i_kappa], linewidth = 0, marker = '*', alpha = 1, ms = 12)
    ax_ranking_i.plot(ranking, fit, color = colors_kappa[i_kappa], linewidth = 5, label = r'$%.d$'%(kappa), alpha = .8)

    my_plot_layout(ax = ax_ranking_i, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_ranking_i.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
    #ax_ranking_i.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
    ax_ranking_i.set_ylim(bottom = 2e-2)
    #ax_ranking_i.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_ranking_i.set_yticklabels([1, 0.1, 0.01])
    fig_ranking_i.savefig('../../Figures/1_Dynamics/Ensemble/L%d/Ranking_p-%.2f'%(L, kappa)+'_'+energy_model+'.pdf')

kappas_theory = np.linspace(1, 4.1, 40)
exponent = []
for kappa in kappas_theory:
    beta_kappa, E_kappa, Kd_kappa = get_p_properties(betas, Q0, Es, dE, kappa)
    beta_act = np.min([beta_r, beta_kappa])
    exponent.append(-kappa*lambda_B/(lambda_A*beta_act))



my_plot_layout(ax = ax_ranking, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_ranking.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
#ax_ranking.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
ax_ranking.set_ylim(bottom = 5e-2)
#ax_ranking.set_yticks([1, 0.1, 0.01, 0.001])
#ax_ranking.set_yticklabels([1, 0.1, 0.01])
fig_ranking.savefig('../../Figures/1_Dynamics/Ensemble/L%d/Ranking_'%L+energy_model+'.pdf')


print('----END-----')




