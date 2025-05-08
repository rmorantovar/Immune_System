import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

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
markers_p = ['.', 'D']
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
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path)
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
fig_size, ax_size = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_affinity, ax_affinity = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_scatter, ax_scatter = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
for i_p, p in enumerate(ps):

    fig_size_i, ax_size_i = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_affinity_i, ax_affinity_i = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_scatter_i, ax_scatter_i = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
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
    data, return_data_type = get_data_b(folder_path = '../out', data_type = 'statistics_p-%.2f'%p)

    if(return_data_type==1):
        clone_size_total = data[0]
        affinities_total = data[1]
    else:
        data = pd.read_csv(Text_files_path + 'primary_immune_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
        clone_size_total = []
        affinities_total = []
        for i_ens in tqdm(np.arange(N_ens)):
            data_active = data.loc[data['ens_id']==i_ens]
            t_act_data = np.min(data_active['time'])
            #data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
            activation_times = np.array(data_active['time'])
            energies  = np.array(data_active['E'])

            #---------------------------- B cell linages ----------------------
            clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
            #--------------------------t_C filter-------------------------
            lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
            clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
            
            clone_size_total = np.concatenate((clone_size_total, clone_sizes_C[:,-1]/np.max(clone_sizes_C[:,-1])))
            affinities_total = np.concatenate((affinities_total, energies_C[:]))

        f = open('../out/processed_data_statistics_p-%.2f.pkl'%p, 'wb')
        pickle.dump([clone_size_total, affinities_total], f, pickle.HIGHEST_PROTOCOL) 


    bins = np.logspace(np.log10(np.min(clone_size_total)*0.5),np.log10(np.max(clone_size_total)*5), 200)
    len_clone_sizes = len(clone_size_total)
    print(len_clone_sizes)
    clone_size_distribution = plt.hist(clone_size_total, bins = bins, density = False, cumulative = True, alpha = 0)
    clone_size = clone_size_distribution[1][:-1]
    clone_size_counts = clone_size_distribution[0]/len_clone_sizes#/np.sum(clone_size_distribution[0]*(np.diff(clone_size_distribution[1])))
    
    bins = np.linspace(np.min(affinities_total), np.max(affinities_total), 200)
    len_affinities = len(affinities_total)
    print(len_affinities)
    affinity_distribution = plt.hist(affinities_total, bins = bins, density = False, cumulative = False, alpha = 0)
    affinity = affinity_distribution[1][:-1]
    affinity_counts = affinity_distribution[0]/len_affinities#/np.sum(clone_size_distribution[0]*(np.diff(clone_size_distribution[1])))
    

    #print(np.sum(clone_size_counts[:]*np.diff(clone_size_distribution[1])))
    
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
    
    ax_size.plot(clone_size/1, 1-clone_size_counts[:], color = colors_p[i_p], linewidth = 2, marker = '', alpha = 1, ms = 4, label = r'$%.d$'%(p))
    ax_size_i.plot(clone_size/1, 1-clone_size_counts[:], color = colors_p[i_p], linewidth = 2, marker = '', alpha = 1, ms = 5, label = r'$%.d$'%(p))

    ax_affinity.plot(np.exp(affinity), affinity_counts[:], color = colors_p[i_p], linewidth = 2, ls = '', marker = 'D', alpha = 1, ms = 4, label = r'$%.d$'%(p))
    ax_affinity_i.plot(np.exp(affinity), affinity_counts[:], color = colors_p[i_p], linewidth = 2, ls = '', marker = 'D', alpha = 1, ms = 5, label = r'$%.d$'%(p))
    ax_affinity.plot(np.exp(Es[:-1]), Q0*len_affinities, color = colors_p[i_p], linewidth = 2, marker = '', alpha = 1, ms = 4)
    ax_affinity_i.plot(np.exp(Es[:-1]), Q0*len_affinities, color = colors_p[i_p], linewidth = 2, marker = '', alpha = 1, ms = 5)

    ax_size.plot(Nb_array/1, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)
    ax_size_i.plot(Nb_array/1, fit, color = colors_p[i_p], linewidth = 4, alpha = .8)

    ax_scatter.scatter(np.exp(affinities_total), clone_size_total, marker = markers_p[i_p], color = colors_p[i_p], alpha = .2)
    ax_scatter_i.scatter(np.exp(affinities_total), clone_size_total, marker = markers_p[i_p], color = colors_p[i_p], alpha = .2)

    my_plot_layout(ax = ax_size_i, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_size_i.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
    ax_size_i.set_xlim(left = 6e-3, right = 1.5)
    ax_size_i.set_ylim(bottom = 9e-4, top = 1.2)
    #ax_size_i.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_size_i.set_yticklabels([1, 0.1, 0.01])
    fig_size_i.savefig('../../Figures/1_Dynamics/CSV/size_p-%.2f_lamAB-%.2f'%(p, lambda_A/lambda_B)+'_'+energy_model+'.pdf')

    my_plot_layout(ax = ax_affinity_i, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_affinity_i.legend(fontsize = 22, title_fontsize = 26, title = r'$p$')
    ax_affinity_i.set_xlim(right = 1e-3, left = 1e-10)
    ax_affinity_i.set_ylim(bottom = 1e-7, top = 10)
    #ax_affinity_i.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_affinity_i.set_yticklabels([1, 0.1, 0.01])
    fig_affinity_i.savefig('../../Figures/1_Dynamics/CSV/affinity_p-%.2f_lamAB-%.2f'%(p, lambda_A/lambda_B)+'_'+energy_model+'.pdf')

    my_plot_layout(ax = ax_scatter_i, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_scatter_i.legend(fontsize = 22, title_fontsize = 26, title = r'$p$')
    ax_scatter_i.set_xlim(right = 1e-3, left = 1e-10)
    ax_scatter_i.set_ylim(bottom = 6e-3, top = 1.5)
    #ax_scatter_i.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_scatter_i.set_yticklabels([1, 0.1, 0.01])
    fig_scatter_i.savefig('../../Figures/1_Dynamics/CSV/scatter_p-%.2f_lamAB-%.2f'%(p, lambda_A/lambda_B)+'_'+energy_model+'.pdf')

#ax_size.hlines(1, 4e-4*C, 1e-1*C, linestyle = 'dashed', color = 'black', linewidth = 1)
my_plot_layout(ax = ax_size, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_size.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
ax_size.set_xlim(left = 6e-3, right = 1.5)
ax_size.set_ylim(bottom = 9e-4, top = 1.2)
#ax_size.set_yticks([1, 0.1, 0.01, 0.001])
#ax_size.set_yticklabels([1, 0.1, 0.01])
fig_size.savefig('../../Figures/1_Dynamics/CSV/size_lamAB-%.2f_'%(lambda_A/lambda_B)+energy_model+'.pdf')

#ax_size.hlines(1, 4e-4*C, 1e-1*C, linestyle = 'dashed', color = 'black', linewidth = 1)
my_plot_layout(ax = ax_affinity, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_affinity.legend(fontsize = 22, title_fontsize = 26, title = r'$p$')
ax_affinity.set_xlim(right = 1e-3, left = 1e-10)
ax_affinity.set_ylim(bottom = 1e-7, top = 10)
#ax_affinity.set_yticks([1, 0.1, 0.01, 0.001])
#ax_affinity.set_yticklabels([1, 0.1, 0.01])
fig_affinity.savefig('../../Figures/1_Dynamics/CSV/affinity_lamAB-%.2f_'%(lambda_A/lambda_B)+energy_model+'.pdf')

#ax_size.hlines(1, 4e-4*C, 1e-1*C, linestyle = 'dashed', color = 'black', linewidth = 1)
my_plot_layout(ax = ax_scatter, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_scatter.legend(fontsize = 22, title_fontsize = 26, title = r'$p$')
ax_scatter.set_xlim(right = 1e-3, left = 1e-10)
ax_scatter.set_ylim(bottom = 6e-3, top = 1.5)
#ax_scatter.set_yticks([1, 0.1, 0.01, 0.001])
#ax_scatter.set_yticklabels([1, 0.1, 0.01])
fig_scatter.savefig('../../Figures/1_Dynamics/CSV/scatter_lamAB-%.2f_'%(lambda_A/lambda_B)+energy_model+'.pdf')
print('----END-----')


