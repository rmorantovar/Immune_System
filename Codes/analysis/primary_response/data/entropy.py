import sys
sys.path.append('../../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------
N_ens = 100
N_enss = [N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens, N_ens]

L_0 = 1e8
transparencies_p = [.8, 1, .8, .8, .8]

T0 = 0
Tf = 10
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

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0


E_lims = [-8.5, -8.5, -12, -12, -12, -12, -12, -12, -12]
t_lims = np.array([4.0, 4.5, 5.2, 5.7, 6.2, 6.7, 7.2, 7.7, 8.2])
ps = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4, 4.5, 5.0]

# E_lims = [-8.5, -8.5, -12, -12]#, -12, -12, -12, -12, -12]
# t_lims = np.array([4.0, 4.5, 5.2, 5.7]) + 1#, 6.2, 6.7, 7.2, 7.7, 8.2]) + 1
# ps = [1, 1.5, 2.0, 2.5]#, 3.0, 3.5, 4, 4.5, 5.0]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])
color_list = np.array([my_blue2, my_blue2, my_green, my_green, my_red, my_red, my_gold, my_gold, my_gold])


colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

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
print('L_0 = %.e'%L_0)
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_pr = %.2f'%beta_pr)

# t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')

#--------------------------Loops--------------------------
# ax_mf = ax_mf.twinx()
print('--------')
print('--- Processing all clones ---')
print('--------')
max_potency_simulations = dict()
max_potency_simulations_std = dict()
max_potency_theory = dict()

print('--- Entropy ---')
#for freq in np.array([0.0500, 0.0100, 0.0050, 0.0010, 0.0001])*10:
for lambda_A in lambda_As:
    print('--------')
    print('lamA=%.1f'%lambda_A)
    Delta_S = dict()
    N_acts = dict()
    fig_mf, ax_mf = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    fig_entropy, ax_entropy = plt.subplots(figsize=(4.5*2,3.0*1.805), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    fig_N_act, ax_N_act = plt.subplots(figsize=(4.5*2,3.0*1.805), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

    for freq in np.array([0.0100])*10:
        # ENTROPYYYYYYYY
        N_acts[freq] = dict()
        max_entropy_simulations = dict()
        max_entropy_simulations_std = dict()
        max_entropy_theory = dict()
        S1 = 1
        S4 = 1
        for i_p, p in enumerate(ps):
            N_ens = N_enss[i_p]
            E_lim = E_lims[i_p]
            t_lim = t_lims[i_p]
            print('p = %.2f...'%p)
            beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
            #-----------------Loading data----------------------------
            parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, N_c, linear, N_ens)+energy_model
            return_data_type = 0
            data, return_data_type = get_data_b(folder_path = '../../../out/primary_immune_response', data_type = 'S_mf_p-%.2f_%.3f_lamA-%.1f_L0-%.0e'%(p, freq, lambda_A, L_0))
        
            if(return_data_type):
                S_final = data[0]
                N_final = data[1]
                N_act = data[2]
                Counter = data[3]
            else:
                data = pd.read_csv(Text_files_path + 'primary_immune_response/simulations_lambda_A/output_N_ens_%d_L0_%d_p_%.1f_k_step_%.1f_E_lim_%.1f_t_lim_%.1f_lamA_%.1f/activated_population_'%(N_ens, L_0, p, k_step, E_lim, t_lim, lambda_A) + antigen + '.csv')
                S_final = []
                N_final = []
                N_act = []
                Counter = 0
                for i_ens in tqdm(np.arange(N_ens)):
                    data_active = data.loc[data['ens_id']==i_ens]
                    #data_active = data_i.loc[data_i['active']==1]
                    t_act_data = np.min(data_active['time'])
                    data_active = data_active.loc[data_active['time']<(t_act_data+1.2+0.3*(p-1))] # it was 1.0 + 0.1*...
                    activation_times = np.array(data_active['time'])
                    energies  = np.array(data_active['E'])
                    #---------------------------- B cell linages ----------------------
                    clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
                    #--------------------------t_C filter-------------------------
                    #lim_size = np.max([int(C*freq), 2])
                    lim_size = np.max([int(np.max(clone_sizes[:, -1])*freq), 2])
                    clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
                    N_act.append(len(activation_times_C))
                    #-------Simulations-------
                    #unactivated_clones = np.array([np.size(activation_times_C[activation_times_C>time_array[t]]) for t in range(len(time))])
                    #total_size = np.sum(clone_sizes_C, axis = 0)
                    #S_i = -np.array([np.sum(((clone_sizes_C[:, t])/total_size[t])*np.log((clone_sizes_C[:, t]-np.size(activation_times_C[activation_times_C>time_array[t]]))/total_size[t])) for t in range(len(time_array))])
                    #S_i = np.nan_to_num(S_i, nan = 0, posinf = 0, neginf = 0)
                    sort_inds = clone_sizes_C[:, -1].argsort()
                    clone_sizes_C_sorted = clone_sizes_C[sort_inds, :][:,:]
                    N_final_i = np.sum(clone_sizes_C_sorted[:, -1])
                    S_final_i = -np.sum((clone_sizes_C_sorted[:, -1]/N_final_i)*np.log(clone_sizes_C_sorted[:, -1]/N_final_i))
                    if(S_final_i!=0):
                        S_final.append(S_final_i)
                        N_final.append(N_final_i)#/np.log(N_final))
                        Counter+=1
                    #if(i_ens%1==0):
                    #   ax_mf.plot(time_array, entropy_i, color = colors_p[  i_p], alpha = .1, linewidth = 1)

                f = open('../../../out/primary_immune_response/processed_data_S_mf_p-%.2f_%.3f_lamA-%.1f_L0-%.0e.pkl'%(p,freq, lambda_A, L_0), 'wb')
                pickle.dump([S_final, N_final, N_act, Counter], f, pickle.HIGHEST_PROTOCOL) 
            
            if(p==1):
                S1 = np.mean(S_final)

            if(p==4):
                S4 = np.mean(S_final)

            S_final = np.array(S_final)
            N_acts[freq][p] = N_act
            max_entropy_simulations[p] = np.mean(S_final)
            max_entropy_simulations_std[p] = np.std(S_final)

        Delta_S[freq] = S1 - S4
        ps_theory = np.linspace(.5, 5, 400)
        y_interp2 = np.interp(ps_theory, ps, np.array(list(max_entropy_simulations.values())))

        ax_entropy.plot(ps, np.array(list(max_entropy_simulations.values())), linestyle = '-', marker = '^', linewidth = 3, label = r'$%.3f$'%freq, ms = 10, alpha = 1)
        ax_N_act.plot(ps, [np.mean(N_acts[freq][p]) for p in ps], linestyle = '-', marker = '^', linewidth = 3, label = r'$%.3f$'%freq, ms = 10, alpha = 1)
        #ax_entropy.plot(ps_theory, y_interp2, linestyle = '-', marker = '', linewidth = 3, ms = 10, alpha = 1)

        if(freq==0.1):
            ax_mf.plot(ps, np.array(list(max_entropy_simulations.values())), color = 'olive', linestyle = '-', marker = '^', linewidth = 3, label = 'Total', ms = 10, alpha = 1)
            #ax_mf.errorbar(x=ps, y=(np.array(list(max_entropy_simulations.values()))), yerr = (np.array(list(max_entropy_simulations_std.values()))), ls = 'none', color = 'olive', alpha = .6)

            if np.any(y_interp2>4.2):
                ax_mf.scatter(ps_theory[y_interp2>4.2][-1], 4.2, s = 180, facecolors='none', edgecolors='olive', marker = 'o', lw=2)
                ax_mf.errorbar(x = ps_theory[y_interp2>4.2][-1], y = 4.2, yerr = 2*0.2, xerr = np.array([[ps_theory[y_interp2>4.2][-1] - ps_theory[y_interp2>(4.2+2*0.2)][-1], ps_theory[y_interp2>(4.2-2*0.2)][-1] - ps_theory[y_interp2>4.2][-1]]]).T, color = 'olive', alpha = .6, fmt='o')
                print(ps_theory[y_interp2>4.2][-1], ps_theory[y_interp2>(4.2+2*0.2)][-1], ps_theory[y_interp2>(4.2-2*0.2)][-1])

    #ax_mf.vlines(beta_r, 2, 6, color = 'k', ls = '--')
    #ax_mf.spines['top'].set_color('olive')
    #ax_mf.xaxis.label.set_color('olive')
    #ax_mf.tick_params(axis='y', colors='olive')

    my_plot_layout(ax = ax_mf, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_mf.legend(fontsize = 28, title_fontsize = 30, loc = 8)
    ax_mf.set_xlim(left = 0.8, right = 5.2)
    #ax_mf.set_ylim(bottom = 0)
    #ax_mf.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_mf.set_yticklabels([1, 0.1, 0.01])

    fig_mf.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2020/entropy_'+energy_model+'_lamA-%.1f_L0-%.0e.pdf'%(lambda_A, L_0))
    plt.close()

