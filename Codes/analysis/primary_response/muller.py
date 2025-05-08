import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
#N_enss = [401, 402, 403]# elite
N_enss = [4000]# aging
#L_0 = 1e8
L_0 = 1e9
T0 = 0
Tf = 9
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

E_lims = [-12]
t_lims = [7.0]
ps = [4.]

t_fs = [4.15] 
#t_fs = [4.37] # for 1e8
#t_fs = [4.93] #for aging

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_red, my_green, my_blue2, my_gold])
#color_list = np.array([my_green])

#colors_p = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_p = np.flip(['tab:blue','tab:green','tab:red'])
colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))
colors_p = ['limegreen']

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0


antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20

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
Es, dE, Q0, betas = calculate_Q0(0.001, 80, 1000000, PWM_data, E_m, L)
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
min_E = -19
max_E = -7.5

#min_E = np.log(1e-9)
max_E = np.log(1e-4)

for i_p, p in enumerate((ps)):
    E_lim = E_lims[i_p]
    t_lim = t_lims[i_p]
    fig_1_2, ax_1_2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    print('--------')
    print('p = %.2f...'%p)
    beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
    beta_act = np.min([beta_r, beta_p])
    #-----------------------------Activation time------------------------
    m_bar_theory = np.array([np.sum(L_0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
    t_act_theory = time_array[m_bar_theory>1][0]
    for N_ens in N_enss:
        #-----------------Loading data----------------------------
        parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
        #data = get_data_ensemble_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path)
        #data = get_data_ensemble_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d_elite/'%L+parameters_path)
        data = pd.read_csv(Text_files_path + 'primary_immune_response/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_E_m_%.1f/filtered_sequence_properties.csv'%(N_ens, L_0, p, k_step, E_lim, t_lim, E_m))
        for i_ens in tqdm(np.arange(N_ens)):
        #for i_ens in tqdm(n_ens): #38 aging, 332 elit ----- 689.0 at 1030

            data_active = data.loc[data['ens_id']==i_ens]
            #data_active = data_i.loc[data_i['active']==1]
            t_act_data = np.min(data_active['time'])
            data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.2*(p-1))] # it was 1.0 + 0.1*...
            activation_times = np.array(data_active['time'])
            energies  = np.array(data_active['E'])
            #---------------------------- B cell linages ----------------------
            clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
            #--------------------------t_C filter-------------------------
            lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
            clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
            Kds_C = np.exp(energies_C)
            Potency = np.sum(((clone_sizes_C[:, -1]-1))/Kds_C)
            #lim_size = 10
            #clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
            sort_inds = activation_times_C.argsort()
            if (Potency>(2.80863588e11/5)):
                #print('----------- FOUND ELITE:', i_ens, 'at', N_ens, '-----------')
                #print(Potency/2.80863588e11)
                #fig_muller, ax_muller = plt.subplots(figsize=(4.5,3), linewidth = 0, gridspec_kw={'left':0.005, 'right':.995, 'bottom':.02, 'top': 0.98}, dpi = 700, edgecolor = 'black')
                #fig_muller, ax_muller = plt.subplots(figsize=(4.5,3), linewidth = 0, gridspec_kw={'left':0.005, 'right':.995, 'bottom':.02, 'top': 0.98}, dpi = 700, edgecolor = 'black')
                fig_muller, ax_muller = plt.subplots(figsize=(5*1.62,5), linewidth = 0, gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96}, dpi = 700, edgecolor = 'black')
                ax_muller.spines["top"].set_linewidth(1)
                ax_muller.spines["left"].set_linewidth(1)
                ax_muller.spines["right"].set_linewidth(1)
                ax_muller.spines["bottom"].set_linewidth(1)
                #---------------------- Total Pop size ----------------------
                total_pop_active = np.sum(clone_sizes_C, axis = 0) #C
                #bcell_freqs = clone_sizes_C/total_pop_active
                bcell_freqs = clone_sizes_C/(np.sum(clone_sizes_C[:,-1]))
                entropy = -np.sum(bcell_freqs*np.log(bcell_freqs), axis = 0)

                #------------------------- Stackplots -------------------------
                greys = plt.cm.get_cmap('cividis', 50)                
                delta_E = max_E - min_E
                for c in np.flip(range(len(clone_sizes_C[:,0]))):
                    color_c = greys(int(50*(1-abs((energies_C[c]-min_E)/delta_E))))
                    ax_muller.stackplot(time_array, [(bcell_freqs[c, -1] - bcell_freqs[c, :])/2 + np.ones_like(bcell_freqs[0, :])*np.sum(bcell_freqs[:c, -1]), bcell_freqs[c, :], (bcell_freqs[c, -1] - bcell_freqs[c, :])/2], colors = ['white', color_c, 'white']);
                    if activation_times_C[c] in activation_times_C[sort_inds[:3]]:
                        ax_muller.scatter(activation_times_C[c], (bcell_freqs[c, -1] - bcell_freqs[c, 0])/2 + np.sum(bcell_freqs[:c, -1]), marker = 'D', edgecolor='black', linewidth=1, facecolor = 'black', s = 50, zorder = 20)

                #ax_muller.vlines([t_act_theory, t_fs[i_p]], 0, 1, color = 'black', linewidth = 2, alpha = 1, linestyles = ['--', ':'])
                ax_muller.vlines(t_act_theory, 0, 1, color = 'black', linewidth = 2, alpha = 1, linestyle = '--')

                my_plot_layout(ax = ax_muller, ticks_labelsize=30, yscale = 'linear')
                ax_muller.set_yticks([])
                #ax_muller.set_xticks(np.arange(Tf))
                ax_muller.set_xticks(range(0, Tf, 2))
                ax_muller.set_xlim(T0, Tf)
                ax_muller.set_ylim(0, 1)
                fig_muller.savefig('../../Figures/1_Dynamics/CSV/Mullers/Clones_p-%.1f_L_0%.e_N_ens-%d_n_ens-%d_Z-%.2e_'%(p, L_0, N_ens, i_ens, Potency/2.80863588e11)+energy_model+'.pdf', edgecolor=fig_muller.get_edgecolor())
                plt.close(fig_muller)

print('----END-----')




