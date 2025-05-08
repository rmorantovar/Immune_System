import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 428# 200 aging
N_enss = [400+i for i in range(1, 32, 1)]
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

kappas = [3.0]
#t_fs = [4.15] 
t_fs = [4.37] # for 1e8
#t_fs = [4.93] #for aging

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_red, my_green, my_blue2, my_gold])
color_list = np.array([my_green])

#colors_kappa = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_kappa = np.flip(['tab:blue','tab:green','tab:red'])
colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))
colors_kappa = ['limegreen']

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
min_E = -19
max_E = -7.5

for i_kappa, kappa in enumerate((kappas)):
    fig_1_2, ax_1_2 = plt.subplots(figsize=(4.5*2,3.0*2), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    print('--------')
    print('kappa = %.2f...'%kappa)
    beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
    beta_act = np.min([beta_r, beta_kappa])
    
    first = []
    second = []
    first_elite = []
    second_elite = []
    for N_ens in N_enss:
        #-----------------Loading data----------------------------
        parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
        #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
        data = get_data_ensemble_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d_elite/'%L+parameters_path)
        for i_ens in tqdm(np.arange(N_ens)):
        #for i_ens in tqdm(n_ens): #38 aging, 332 elit ----- 689.0 at 1030

            data_active = data.loc[data['i_ens']==i_ens]
            #data_active = data_i.loc[data_i['active']==1]
            t_act_data = np.min(data_active['act_time'])
            data_active = data_active.loc[data_active['act_time']<(t_act_data+1.2+0.3*(kappa-1))] # it was 1.0 + 0.1*...
            activation_times = np.array(data_active['act_time'])
            energies  = np.array(data_active['energy'])

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
            if(len(clone_sizes_C[sort_inds, -1])==1):
                first.append(clone_sizes_C[sort_inds, -1][0])
                second.append(0)
            else:
                first.append(clone_sizes_C[sort_inds, -1][0])
                second.append(clone_sizes_C[sort_inds, -1][1])
            if (Potency>2.80863588e11):
                if(len(clone_sizes_C[sort_inds, -1])==1):
                    first_elite.append(clone_sizes_C[sort_inds, -1][0])
                    second_elite.append(0)
                else:
                    first_elite.append(clone_sizes_C[sort_inds, -1][0])
                    second_elite.append(clone_sizes_C[sort_inds, -1][1])

    ax_1_2.scatter(np.array(first_elite)/C, np.array(second_elite)/C, color = my_red, label = r'$%.2f$'%kappa, alpha = 1)
    ax_1_2.scatter(np.array(first)/C, np.array(second)/C, color = my_red, alpha = .05)
    ax_1_2.plot(np.linspace(10/C, 1, 100), np.linspace(10/C, 1, 100)*2**(-(lambda_B*kappa)/(lambda_A*beta_r)), color = 'k')
    ax_1_2.plot(np.linspace(10/C, 1, 100), np.linspace(10/C, 1, 100), color = 'grey')

    my_plot_layout(ax = ax_1_2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_1_2.legend(fontsize = 28, title_fontsize = 30, title = r'$p$', loc = 4)
    #ax_1_2.set_xlim(left = 2, right = Tf-1)
    #ax_1_2.set_ylim(bottom = 1e-3, top = 1e0)
    #ax_1_2.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_1_2.set_yticklabels([1, 0.1, 0.01])
    fig_1_2.savefig('../../Figures/1_Dynamics/Ensemble/L%d/Mullers/1_2_'%(L)+energy_model+'.pdf')
print('----END-----')




