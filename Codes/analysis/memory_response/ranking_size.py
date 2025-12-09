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
#E_m = -27.63
E_m = -24
E_lims = [-8.5, -12]
t_lims = [4.0, 7.2]
E_lims = [-12]
t_lims = [7.2]
# E_lims = [-8.5, -8.5, -12, -12, -12, -12, -12, -12, -12]
# t_lims = [4.0, 4.5, 5.2, 5.7, 6.2, 6.7, 7.2, 7.7, 8.2]
C = 1e4
AA = 1


ps = [1, 2, 2.5, 3, 4]
ps = [1, 4]
ps = [4]
# ps = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4, 4.5, 5.0]

transparency_n = [1]

color_list = np.array([my_blue, my_red, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([my_blue2, my_red])
#color_list = np.array([my_blue2, my_blue2, my_green, my_green, my_red, my_red, my_gold])
color_list = np.array([my_red, my_blue2])

colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

times = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0


project = 'memory_response'
subproject = 'clonal_structure'
output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
os.makedirs(output_plot, exist_ok=True)


#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18'

L=len(antigen)

#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path + "in/")
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)
#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_m, L)
Kds = np.exp(Es[:-1])
#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------

for i_p, p in enumerate((ps)):
    E_lim = E_lims[i_p]
    t_lim = t_lims[i_p]
    print('p = %.2f...'%p)
    beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
    beta_act = np.min([beta_r, beta_p])

    #-----------------Loading data----------------------------
    t_cst = 3.1 # loop between 3.1 and 4.4
    for t_cst in [3.1, 3.5, 4.4]:
        fig_ranking_i, ax_ranking_i = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
        fig_CSD_i, ax_CSD_i = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
        fig_CSD_2_i, ax_CSD_2_i = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

        return_data_type = 0
        data, return_data_type = get_data(folder_path = Text_files_path + 'memory_response/out/', data_type = 'ranking_size_p-%.1f_t_cst-%.1f'%(p, t_cst))
        # return_data_type = 0
        n_first_clones = 100
        
        if(return_data_type):
            final_Nb = data[0]
            counts_final_Nb = data[1]
            trajectories = data[2]
            trajectories_rank = data[3]
            clone_size_total = data[4]
            final_Nb_cst = data[5]
            counts_final_Nb_cst = data[6]
            trajectories_cst = data[7]
            trajectories_rank_cst = data[8]
            clone_size_total_cst = data[9]
            growth_rates_total_cst = data[10]
        else:
            #activation_times_total = np.array([])
            data = pd.read_csv(Text_files_path + 'primary_response/simulations/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))    
            print(f'Processing data ... ')
            final_Nb = np.zeros(n_first_clones)
            counts_final_Nb = np.zeros(n_first_clones)
            trajectories = np.array([], dtype = object)
            trajectories_rank = np.array([])
            clone_size_total = np.array([])

            final_Nb_cst = np.zeros(n_first_clones)
            counts_final_Nb_cst = np.zeros(n_first_clones)
            trajectories_cst = np.array([], dtype = object)
            trajectories_rank_cst = np.array([])
            clone_size_total_cst = np.array([])
            growth_rates_total_cst = np.array([])

            max_rank = 100

            for i_ens in tqdm(np.arange(N_ens)):
                data_active = data.loc[data['ens_id']==i_ens]
                t_act_data = np.min(data_active['time'])
                data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
                activation_times = np.array(data_active['time'])
                energies  = np.array(data_active['E'])

                #---------------------------- B cell linages ----------------------
                clone_sizes = get_clones_sizes_C(len(activation_times), times, activation_times, lambda_B, C, dT)
                #--------------------------t_C filter-------------------------
                lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
                clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
                clone_size_total = np.concatenate((clone_size_total, clone_sizes_C[:,-1]))

                sort_inds = clone_sizes_C[:, -1].argsort()
                clone_sizes_C_sorted = clone_sizes_C[sort_inds, :][-int(n_first_clones*(4-3)):, :]            
                biggest_clone_i = clone_sizes_C_sorted[-1, -1]
                sorted_clones = np.flip(clone_sizes_C_sorted[:, -1])/biggest_clone_i
                max_rank_i = len(sorted_clones)
                if max_rank_i>1:
                    for i in range(max_rank_i):
                        final_Nb[i]+= (sorted_clones[i])
                        counts_final_Nb[i] += 1
                    if(max_rank_i<max_rank):
                        max_rank = max_rank_i
                    if((i_ens%10==0) and (p==4.0)):
                        trajectories = np.append(trajectories, sorted_clones)
                        trajectories_rank = np.append(trajectories_rank, max_rank_i)

                #---------------------------- B cell linages cst----------------------
                growth_rates_cst = lambda_B*(1+(np.exp(energies_C)/(np.exp(lambda_A*(t_prime+t_cst))/N_A)))**(-1)
                clone_sizes_cst = get_clones_sizes_C(len(activation_times_C), times, np.ones(len(activation_times_C)), growth_rates_cst, C, dT)
                #--------------------------t_C filter cst-------------------------
                lim_size = np.max([int(np.max(clone_sizes_cst[:, -1])*0.01), 2])
                clone_sizes_C_cst, activation_times_C_cst, energies_C_cst, filter_C_cst, n_C_cst = apply_filter_C(clone_sizes_cst, activation_times_C, energies_C, lim_size)
                clone_size_total_cst = np.concatenate((clone_size_total_cst, clone_sizes_C_cst[:,-1]))
                growth_rates_total_cst = np.concatenate((growth_rates_total_cst, growth_rates_cst[filter_C_cst]))

                sort_inds_cst = clone_sizes_C_cst[:, -1].argsort()
                clone_sizes_C_sorted_cst = clone_sizes_C_cst[sort_inds_cst, :][-int(n_first_clones*(4-3)):, :]
                biggest_clone_i_cst = clone_sizes_C_sorted_cst[-1, -1]
                sorted_clones_cst = np.flip(clone_sizes_C_sorted_cst[:, -1])/biggest_clone_i_cst
                max_rank_i_cst = len(sorted_clones_cst)
                if max_rank_i_cst>1:
                    for i in range(max_rank_i_cst):
                        final_Nb_cst[i]+= (sorted_clones_cst[i])
                        counts_final_Nb_cst[i] += 1
                    if(max_rank_i_cst<max_rank):
                        max_rank = max_rank_i_cst
                    if((i_ens%10==0) and (p==4.0)):
                        trajectories_cst = np.append(trajectories_cst, sorted_clones_cst)
                        trajectories_rank_cst = np.append(trajectories_rank_cst, max_rank_i_cst)

            f = open(Text_files_path + 'memory_response/out/processed_data_ranking_size_p-%.1f_t_cst-%.1f.pkl'%(p, t_cst), 'wb')
            pickle.dump([final_Nb, counts_final_Nb, trajectories, trajectories_rank, clone_size_total, final_Nb_cst, counts_final_Nb_cst, trajectories_cst, trajectories_rank_cst, clone_size_total_cst, growth_rates_total_cst], f, pickle.HIGHEST_PROTOCOL) 
        
        counter = 0
        for j in range(len(trajectories_rank)):
            ranks_j = np.arange(1, trajectories_rank[j]+1)
            len_rank_j = len(ranks_j)
            ax_ranking_i.plot(ranks_j, trajectories[counter:counter+len_rank_j], color = colors_p[i_p], linewidth = 1, alpha = .2)
            counter += len_rank_j
        final_Nb = final_Nb/counts_final_Nb
        ranking = np.arange(1, n_first_clones+1)
        fit = ranking**(-p*lambda_B/(lambda_A*beta_act))
        ax_ranking_i.plot(ranking, final_Nb, color = colors_p[i_p], linewidth = 0, marker = '*', alpha = 1, ms = 10, label = r'$\textrm{EPR}$')
        ax_ranking_i.plot(ranking, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)

        counter_cst = 0
        for j in range(len(trajectories_rank_cst)):
            ranks_j = np.arange(1, trajectories_rank_cst[j]+1)
            len_rank_j = len(ranks_j)
            ax_ranking_i.plot(ranks_j, trajectories_cst[counter_cst:counter_cst+len_rank_j], color = colors_p[i_p+1], linewidth = 1, alpha = .2)
            counter_cst += len_rank_j
        final_Nb_cst = final_Nb_cst/counts_final_Nb_cst
        ax_ranking_i.plot(ranking, final_Nb_cst, color = colors_p[i_p+1], linewidth = 0, marker = '*', alpha = 1, ms = 10, label = r'$\textrm{CST}$')

        my_plot_layout(ax = ax_ranking_i, xscale='log', yscale= 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
        ax_ranking_i.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
        #ax_ranking_i.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
        ax_ranking_i.set_ylim(bottom = 2e-2, top = 1.05)
        #ax_ranking_i.set_yticks([1, 0.1, 0.01, 0.001])
        #ax_ranking_i.set_yticklabels([1, 0.1, 0.01])
        fig_ranking_i.savefig(output_plot + '/ranking_clone-size_p-%.2f_t_cst-%.1f.pdf'%(p, t_cst))

        bins = np.logspace(np.log10(np.min(clone_size_total)*0.5), np.log10(np.max(clone_size_total)*10), 200)
        len_clone_sizes = len(clone_size_total)
        print(len_clone_sizes)
        clone_size_distribution = plt.hist(clone_size_total, bins = bins, density = False, cumulative = True, alpha = 0)
        clone_size = clone_size_distribution[1][:-1]
        clone_size_counts = clone_size_distribution[0]/len_clone_sizes#/np.sum(clone_size_distribution[0]*(np.diff(clone_size_distribution[1])))
        
        print(np.sum(clone_size_counts[:]*np.diff(clone_size_distribution[1])))
        
        Nb_array = np.logspace(np.log10(np.min(clone_size_total)), np.log10(np.max(clone_size_total)), 50)
        fit = Nb_array**(-beta_act*lambda_A/(lambda_B*p))
        if(p==1):
            fit = Nb_array**(-beta_act*lambda_A/(lambda_B*p))
            fit = fit/fit[-1]*1.2e-3#*np.sum(clone_size_counts[:])*0.8
        else:
            fit = Nb_array**(-beta_act*lambda_A/(lambda_B*p))
            fit = fit/fit[-1]*1.1e-5#*np.sum(clone_size_counts[:])*0.8
        normalization = len(clone_size_total)
        normalization = 1
        #ax_CSD_i.plot(clone_size/1, 1-np.cumsum(clone_size_counts[:]*np.diff(clone_size_distribution[1]))/np.sum(clone_size_counts[:]*np.diff(clone_size_distribution[1])), color = colors_p[i_p], linewidth = 0, marker = 's', alpha = 1, ms = 5, label = r'$%.d$'%(p)) 
        ax_CSD_i.plot(clone_size/1, 1-clone_size_counts[:], color = colors_p[i_p], linewidth = 2, marker = '', alpha = 1, ms = 5, label = r'$\textrm{EPR}$')
        ax_CSD_i.plot(Nb_array/1, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)

        bins = np.logspace(np.log10(np.min(clone_size_total_cst)*0.5), np.log10(np.max(clone_size_total_cst)*10), 200)
        len_clone_sizes_cst = len(clone_size_total_cst)
        print(len_clone_sizes_cst)
        clone_size_distribution_cst = plt.hist(clone_size_total_cst, bins = bins, density = False, cumulative = True, alpha = 0)
        clone_size_cst = clone_size_distribution_cst[1][:-1]
        clone_size_counts_cst = clone_size_distribution_cst[0]/len_clone_sizes_cst#/np.sum(clone_size_distribution[0]*(np.diff(clone_size_distribution[1])))
        
        #ax_CSD_i.plot(clone_size/1, 1-np.cumsum(clone_size_counts[:]*np.diff(clone_size_distribution[1]))/np.sum(clone_size_counts[:]*np.diff(clone_size_distribution[1])), color = colors_p[i_p], linewidth = 0, marker = 's', alpha = 1, ms = 5, label = r'$%.d$'%(p)) 
        ax_CSD_i.plot(clone_size_cst/1, 1-clone_size_counts_cst[:], color = colors_p[i_p+1], linewidth = 2, marker = '', alpha = 1, ms = 5, label = r'$\textrm{CST}$')

        my_plot_layout(ax = ax_CSD_i, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
        ax_CSD_i.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
        ax_CSD_i.set_xlim(right = 1e3, left = 1)
        ax_CSD_i.set_ylim(bottom = 9e-5, top = 1.2)
        #ax_CSD_i.set_yticks([1, 0.1, 0.01, 0.001])
        #ax_CSD_i.set_yticklabels([1, 0.1, 0.01])
        fig_CSD_i.savefig(output_plot + '/clone-size_distribution_p-%.2f_lamAB-%.2f_t_cst-%.1f.pdf'%(p, lambda_A/lambda_B, t_cst))

        #### GROWTH RATES ####
        bins = np.logspace(np.log10(5e-2*0.5),np.log10(lambda_B*2), 40)

        ax_CSD_2_i.hist(growth_rates_total_cst, bins = bins, density = True, alpha = 0.7, color = colors_p[i_p+1], label = r'$\textrm{CST}$')

        my_plot_layout(ax = ax_CSD_2_i, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
        # ax_CSD_2_i.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
        ax_CSD_2_i.set_xlim(right = lambda_B*1.2, left = 0.05)
        ax_CSD_2_i.set_ylim(bottom = 1e-1, top = 2e2)
        #ax_CSD_2_i.set_yticks([1, 0.1, 0.01, 0.001])
        #ax_CSD_2_i.set_yticklabels([1, 0.1, 0.01])
        fig_CSD_2_i.savefig(output_plot + '/growth_rates_p-%.2f_lamAB-%.2f_t_cst-%.1f.pdf'%(p, lambda_A/lambda_B, t_cst))
   

print('----END-----')




