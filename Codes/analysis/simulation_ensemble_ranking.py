import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 10
N_r = 1e8
T0 = 3
Tf = 8
Tf_sim = 7
#Tf = 10
dT = 0.01
lambda_A = 6
k_pr = 1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [3, 2]

transparency_n = [1]

colors_kappa = ['lightskyblue', 'tab:cyan','tab:green', 'tab:red']
colors_kappa = ['tab:cyan','tab:green', 'tab:red']
colors_R = [['deepskyblue', 'lightskyblue', 'lightskyblue'], ['tab:purple', 'tab:cyan', 'tab:cyan'], ['tab:blue', 'tab:green', 'tab:green'], ['tab:red', 'tab:red', 'tab:red']]
colors_R = [['tab:purple', 'tab:cyan', 'tab:cyan'], ['tab:blue', 'tab:green', 'tab:green'], ['tab:red', 'tab:red', 'tab:red']]

lambda_B = lambda_A/2
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
#N_c = 1e5
E_ms = -27.63
C = 3e4

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
antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'
#antigen = 'TACNSEYPNTTKCGRWYC'
L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data = get_motif(antigen, energy_model, Text_files_path)
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
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')
#--------------------------Loops--------------------------
fig_N_K, ax_N_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
for i_kappa, kappa in enumerate((kappas)):
    print('--------')
    print('kappa = %.2f...'%kappa)
    beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
    beta_act = np.min([beta_r, beta_kappa])

    #-----------------Loading data----------------------------
    parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, linear, N_ens)+energy_model
    #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
    data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)

    #activation_times_total = np.array([])
    final_Nb = np.zeros(40)
    max_rank = 0
    for i_ens in tqdm(np.arange(N_ens)):
        data_i = data.loc[data[4]==i_ens]
        data_active = data_i.loc[data_i[1]==1]
        t_act_data = np.min(data_active[3])
        data_active = data_active.loc[data_active[3]<(t_act_data+1.3)]
        activation_times = np.array(data_active[3])
        energies  = np.array(data_active[0])

        #---------------------------- B cell linages ----------------------
        clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

        #--------------------------t_C filter-------------------------
        lim_size = 1
        clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
        
        sort_inds = clone_sizes_C[:, -1].argsort()
        clone_sizes_C_sorted = clone_sizes_C[sort_inds, :][-int(40*(4-3)):, :]
        activation_times_C_sorted = activation_times_C[sort_inds][-int(40*(4-3)):]
        energies_C_sorted = energies_C[sort_inds][-int(40*(4-3)):]

        biggest_clone_i = clone_sizes_C_sorted[-1, -1]
        #ax_N_K.scatter(np.exp(energies_C_sorted[-1]), biggest_clone_i, color = colors_kappa[i_kappa], alpha = .25, edgecolor='black', linewidth=1, facecolor = colors_kappa[i_kappa])
        #activation_times_total = np.append(activation_times_total, activation_times_C_sorted)
        sorted_clones = np.flip(clone_sizes_C_sorted[:, -1])/clone_sizes_C_sorted[-1, -1]
        max_rank_i = len(sorted_clones)
        for i in range(max_rank_i):
            final_Nb[i]+= sorted_clones[i]
        if(max_rank_i>max_rank):
            max_rank = max_rank_i
        ax_N_K.step(np.arange(1, max_rank_i+1), sorted_clones, color = colors_kappa[i_kappa], linewidth = 1, alpha = .2)

    final_Nb/=N_ens
    #final_Nb = np.exp(lambda_B*(Tf-activation_times_total))
    #a, b = np.polyfit(energies_total, np.log(final_Nb), 1)
    #print('Slope from simulation=%.2f'%(a))
    #print('Expected slope=%.2f'%(-kappa*lambda_B/lambda_A))
    ranking = np.arange(1, 40+1)
    fit = ranking**(-kappa*lambda_B/(lambda_A*beta_act))
    ax_N_K.plot(ranking, final_Nb, color = colors_kappa[i_kappa], linewidth = 0, marker = '*', alpha = 1)
    ax_N_K.plot(ranking, fit, color = colors_kappa[i_kappa], linewidth = 5, label = r'$%.d$'%(kappa))

my_plot_layout(ax = ax_N_K, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_N_K.legend(fontsize = 32, title_fontsize = 34, title = r'$\kappa$')
#ax_N_K.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#ax_N_K.set_ylim(bottom = 1e-4)
#ax_N_K.set_yticks([1, 0.1, 0.01, 0.001])
#ax_N_K.set_yticklabels([1, 0.1, 0.01])
fig_N_K.savefig('../../Figures/1_Dynamics/Ensemble/Ranking_'+energy_model+'.pdf')
print('----END-----')




