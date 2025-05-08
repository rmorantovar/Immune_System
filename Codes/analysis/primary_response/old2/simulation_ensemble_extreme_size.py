import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 100
N_r = 2e8
T0 = 3
Tf = 10
Tf_sim = 7.0
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2, 3]
kappas = [2, 3, 4]


my_red = np.array((228,75,41))/256.
my_purple = np.array((125,64,119))/256.
my_purple2 = np.array((116,97,164))/256.
my_green = np.array((125,165,38))/256.
my_blue = np.array((76,109,166))/256.
my_gold = np.array((215,139,45))/256.
my_brown = np.array((182,90,36))/256.
my_blue2 = np.array((80,141,188))/256.
my_yellow = np.array((246,181,56))/256.
my_green2 = np.array((158,248,72))/256.
my_cyan = 'tab:cyan'

antigen_color = my_yellow/256.

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_red, my_green, my_blue2, my_gold])
color_list = np.array([my_green, my_blue2, my_gold])

#colors_kappa = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_kappa = np.flip(['tab:blue','tab:green','tab:red'])
colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

lambda_B = lambda_A/2
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5
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

avg_E = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]) + E_ms
var_E = np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])

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
fig_biggest_all, ax_biggest_all = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

biggest_clones = []

for i_kappa, kappa in enumerate((kappas)):
    fig_biggest, ax_biggest = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    print('--------')
    print('kappa = %.2f...'%kappa)
    beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
    beta_act = np.min([beta_r, beta_kappa])

    #-----------------Loading data----------------------------
    parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
    #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
    data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)

    #activation_times_total = np.array([])
    biggest_clones_p = []
    best_clones_p = []
    for i_ens in tqdm(np.arange(N_ens)):
        data_i = data.loc[data[4]==i_ens]
        data_active = data_i.loc[data_i[1]==1]
        t_act_data = np.min(data_active[3])
        data_active = data_active.loc[data_active[3]<(t_act_data+1.4)]
        activation_times = np.array(data_active[3])
        energies  = np.array(data_active[0])

        #---------------------------- B cell linages ----------------------
        clone_sizes = get_clones_sizes_C(len(activation_times), time, activation_times, lambda_B, C, dT)

        #--------------------------t_C filter-------------------------
        lim_size = 2
        clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
        
        biggest_clones_p.append(np.max(clone_sizes_C[:, -1]))

        biggest_clones.append(np.max(clone_sizes_C[:, -1]))

    
    data_biggest = np.histogram(np.array(biggest_clones_p), bins = np.logspace(np.log10(np.min(biggest_clones_p)), np.log10(np.max(biggest_clones_p)), 20), density = False)

    ax_biggest.plot(data_biggest[1][:-1], data_biggest[0]/len(biggest_clones_p), color = colors_kappa[i_kappa], label = r'$%d$'%kappa)
    ax_biggest_all.plot(data_biggest[1][:-1], data_biggest[0]/len(biggest_clones_p), color = colors_kappa[i_kappa], label = r'$%d$'%kappa, alpha = .8)

    my_plot_layout(ax = ax_biggest, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_biggest.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
    #ax_biggest.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
    #ax_biggest.set_ylim(bottom = 2e-2)
    #ax_biggest.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_biggest.set_yticklabels([1, 0.1, 0.01])
    fig_biggest.savefig('../../Figures/1_Dynamics/Ensemble/Biggest_p-%.2f'%(kappa)+'_'+energy_model+'.pdf')


my_plot_layout(ax = ax_biggest_all, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_biggest_all.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
#ax_biggest_all.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#ax_biggest_all.set_ylim(bottom = 2e-2)
#ax_biggest_all.set_yticks([1, 0.1, 0.01, 0.001])
#ax_biggest_all.set_yticklabels([1, 0.1, 0.01])
fig_biggest_all.savefig('../../Figures/1_Dynamics/Ensemble/Biggest_'+energy_model+'.pdf')

print('----END-----')




