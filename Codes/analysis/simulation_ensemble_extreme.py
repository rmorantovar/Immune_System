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
Tf = 8
Tf_sim = 6.5
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = [1, 2, 3]
kappas = [1, 2, 3]


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
fig_best_all, ax_best_all = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

best_clones = []

ax_best_all.plot(Kds,  Q0, alpha = 1, color = 'grey', linewidth = 5, linestyle = '--')

for i_kappa, kappa in enumerate((kappas)):
    fig_best, ax_best = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    print('--------')
    print('kappa = %.2f...'%kappa)
    beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
    beta_act = np.min([beta_r, beta_kappa])
    ax_best.plot(Kds, Q0, alpha = 1, color = 'grey', linewidth = 5, linestyle = '--')

    #-----------------Loading data----------------------------
    parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
    #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
    data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)

    #activation_times_total = np.array([])
    best_clones_p = []
    for i_ens in tqdm(np.arange(N_ens)):
        data_i = data.loc[data[4]==i_ens]
        data_active = data_i.loc[data_i[1]==1]
        t_act_data = np.min(data_active[3])
        data_active = data_active.loc[data_active[3]<(t_act_data+1.4)]
        activation_times = np.array(data_active[3])
        energies  = np.array(data_active[0])

        best_clones_p.append(np.min(energies))
        best_clones.append(np.min(energies))


    bata_best = np.histogram(np.exp(best_clones_p), bins = np.logspace(np.log10(3e-9), np.log10(7e-8), 10), density = False)
    ax_best.plot(bata_best[1][:-1], bata_best[0]/len(best_clones_p), linestyle = '', marker = 's', color = 'black', ms = 12)
    #ax_best.plot(Kds, P_min_e(N_r, avg_E, var_E, Es[:-1], dE), linestyle = '--', marker = '',  color = 'black', ms = 2, linewidth = 4)
    ax_best.plot(Kds, P_min_e_Q0(N_r, Q0, dE), linestyle = '--', marker = '',  color = 'black', ms = 2, linewidth = 4, alpha = .8)

    my_plot_layout(ax = ax_best, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_best.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
    #ax_best.set_xlim(right = 1e-2)
    ax_best.set_ylim(bottom = 1e-3, top = 1)
    #ax_best.set_yticks([1, 0.1, 0.01, 0.001])
    #ax_best.set_yticklabels([1, 0.1, 0.01])
    fig_best.savefig('../../Figures/1_Dynamics/Ensemble/Best_p-%.2f'%(kappa)+'_'+energy_model+'.pdf')

bata_best = np.histogram(np.exp(best_clones), bins = np.logspace(np.log10(3e-9), np.log10(7e-8), 10), density = False)
ax_best_all.plot(bata_best[1][:-1], bata_best[0]/len(best_clones), linestyle = '', marker = 's', color = 'black', ms = 12)
#ax_best_all.plot(Kds, P_min_e(N_r, avg_E, var_E, Es[:-1], dE), linestyle = '--', marker = '',  color = 'black', ms = 2, linewidth = 4)
ax_best_all.plot(Kds, P_min_e_Q0(N_r, Q0, dE), linestyle = '--', marker = '',  color = 'black', ms = 2, linewidth = 4, alpha = .8)

my_plot_layout(ax = ax_best_all, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_best_all.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
#ax_best_all.set_xlim(right = 1e-2)
ax_best_all.set_ylim(bottom = 1e-3, top = 1)
#ax_best_all.set_yticks([1, 0.1, 0.01, 0.001])
#ax_best_all.set_yticklabels([1, 0.1, 0.01])
fig_best_all.savefig('../../Figures/1_Dynamics/Ensemble/best'+energy_model+'.pdf')

print('----END-----')




