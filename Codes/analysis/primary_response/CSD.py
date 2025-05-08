import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'


#--------------- PARAMETERS ---------------------
N_enss = [400, 400]
L_0 = 1e9
T0 = 0
Tf = 10
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

ps = [2.2, 2.0, 1.8, 1.5]#, 1]
ps = [1.4, 1.8, 2.2]
ps = [2, 3, 4]
ps = [1, 4]

E_lims = [-8.5, -12]
t_lims = [4.0, 7.2]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_blue2, my_red, my_brown])
#color_list = np.array([my_green])

#colors_p = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_p = np.flip(['tab:blue','tab:green','tab:red'])
colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
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

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_CSD, ax_CSD = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
for i_p, p in enumerate(ps):

    fig_CSD_i, ax_CSD_i = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    N_ens = N_enss[i_p]
    E_lim = E_lims[i_p]
    t_lim = t_lims[i_p]
    markers_b = ['^', 's', 'o', '*']
    print('--------')
    print('p = %.2f...'%p)
    beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
    beta_act = np.min([beta_r, beta_p])
    #-----------------Loading data----------------------------
    parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
    #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
    return_data_type = 0
    data, return_data_type = get_data_b(folder_path = '../../out/primary_immune_response', data_type = 'CSD_p-%.2f'%p)

    if(return_data_type==1):
        clone_size_total = data[0]
    else:
        data = pd.read_csv(Text_files_path + 'primary_immune_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
        clone_size_total = []
        for i_ens in tqdm(np.arange(N_ens)):
            data_active = data.loc[data['ens_id']==i_ens]
            t_act_data = np.min(data_active['time'])
            data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
            activation_times = np.array(data_active['time'])
            energies  = np.array(data_active['E'])

            #---------------------------- B cell linages ----------------------
            clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
            #--------------------------t_C filter-------------------------
            lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
            clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
            
            clone_size_total = np.concatenate((clone_size_total, clone_sizes_C[:,-1]/np.max(clone_sizes_C[:,-1])))

        f = open('../../out/primary_immune_response/processed_data_CSD_p-%.2f.pkl'%p, 'wb')
        pickle.dump([clone_size_total], f, pickle.HIGHEST_PROTOCOL) 

    bins = np.logspace(np.log10(np.min(clone_size_total)*0.5),np.log10(np.max(clone_size_total)*5), 200)
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
        fit = fit/fit[-1]*1.2e-3#*np.sum(clone_size_counts[:])*0.8
    normalization = len(clone_size_total)
    normalization = 1
    #ax_CSD.plot(clone_size/1, 1-np.cumsum(clone_size_counts[:]*np.diff(clone_size_distribution[1]))/np.sum(clone_size_counts[:]*np.diff(clone_size_distribution[1])), color = colors_p[i_p], linewidth = 0, marker = markers_b[i_lambda_B], alpha = 1, ms = 4, label = r'$%.d$'%(p)) 
    #ax_CSD.plot(clone_size[:], (clone_size_counts[:]*(clone_size_distribution[1][1:]-clone_size_distribution[1][:-1])), color = colors_p[i_p], linewidth = 0, marker = 's', alpha = .8, ms = 10)
    ax_CSD.plot(clone_size/1, 1-clone_size_counts[:], color = colors_p[i_p], linewidth = 2, marker = '', alpha = 1, ms = 4, label = r'$%.d$'%(p))
    #ax_CSD_i.plot(clone_size/1, 1-np.cumsum(clone_size_counts[:]*np.diff(clone_size_distribution[1]))/np.sum(clone_size_counts[:]*np.diff(clone_size_distribution[1])), color = colors_p[i_p], linewidth = 0, marker = 's', alpha = 1, ms = 5, label = r'$%.d$'%(p)) 
    #ax_CSD.plot(clone_size[:], (clone_size_counts[:]*(clone_size_distribution[1][1:]-clone_size_distribution[1][:-1])), color = colors_p[i_p], linewidth = 0, marker = 's', alpha = .8, ms = 10)
    ax_CSD_i.plot(clone_size/1, 1-clone_size_counts[:], color = colors_p[i_p], linewidth = 2, marker = '', alpha = 1, ms = 5, label = r'$%.d$'%(p))

    ax_CSD.plot(Nb_array/1, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)
    ax_CSD_i.plot(Nb_array/1, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)

    my_plot_layout(ax = ax_CSD_i, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_CSD_i.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
    ax_CSD_i.set_xlim(right = 1.5)
    ax_CSD_i.set_ylim(bottom = 9e-4, top = 1.2)
    #ax_CSD_i.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_CSD_i.set_yticklabels([1, 0.1, 0.01])
    fig_CSD_i.savefig('../../../Figures/primary_immune_response/1_Dynamics/CSV/CSD_p-%.2f_lamAB-%.2f'%(p, lambda_A/lambda_B)+'_'+energy_model+'.pdf')

#ax_CSD.hlines(1, 4e-4*C, 1e-1*C, linestyle = 'dashed', color = 'black', linewidth = 1)
my_plot_layout(ax = ax_CSD, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_CSD.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
ax_CSD.set_xlim(right = 1.5)
ax_CSD.set_ylim(bottom = 9e-4, top = 1.2)
#ax_CSD.set_yticks([1, 0.1, 0.01, 0.001])
#ax_CSD.set_yticklabels([1, 0.1, 0.01])
fig_CSD.savefig('../../../Figures/primary_immune_response/1_Dynamics/CSV/CSD_lamAB-%.2f_'%(lambda_A/lambda_B)+energy_model+'.pdf')
print('----END-----')


