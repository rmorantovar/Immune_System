import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 500
N_rs = [1e7, 1e8]#, 1e9]
T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
#k_pr = 120
k_pr = k_pr*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_cs = [1e5*10000, 1e5*1000]
#N_c = 1e5
#E_ms = -27.63
E_ms = -25
C = 1e4
AA = 1

kappas = [1, 2, 3, 4]
p = 3
#kappas = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_green, my_red, my_gold])
color_list = np.array([my_blue2, my_red, my_green])
#color_list = np.array([my_blue2, my_blue2, my_green, my_green, my_red, my_red, my_gold])

colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

colors_R = []
for i in range(len(N_rs)):
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

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_a = %.2f'%beta_pr, 'K_a = %.2e'%Kd_pr)
print('--------')
print('Loops...')
#--------------------------Loops--------------------------

fig_exponents, ax_exponents = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

exponent_sim = []
exponent_sim2 = []
std_exponent_sim = []
std_exponent_sim2 = []

for i_N_r, N_r in enumerate(N_rs):
    N_c = N_cs[i_N_r]
    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
    print('beta_r = %.1f'%beta_r)
    t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
    print('--------')
    print('p = %.2f...'%p)
    beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
    beta_act = np.min([beta_r, beta_p])

    #-----------------Loading data----------------------------
    parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, p, N_c, linear, N_ens)+energy_model
    #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
    #data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
    data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'ranking_size_N0')
    n_first_clones = 50

    if(return_data_type):
        final_Nb = data[0]
        counts_final_Nb = data[1]
        trajectories = data[2]
        trajectories_rank = data[3]
    else:        
        print('Please run "simulation_ensemble_ranking_size.py" script first ...')
    
    final_Nb = final_Nb/counts_final_Nb
    ranking = np.arange(1, n_first_clones+1)
    popt, pcov = curve_fit(my_linear_func, np.log(ranking[~np.isnan(final_Nb)][:10]), np.log(final_Nb[~np.isnan(final_Nb)][:10]))

    exponent_sim.append(-popt[1])
    std_exponent_sim.append(np.sqrt(pcov[1,1]))

    print('--------')

    beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
    beta_act = np.min([beta_r, beta_p])

    #-----------------Loading data----------------------------
    parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, p, N_c, linear, N_ens)+energy_model
    #data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
    #data = get_data_ensemble(folder_path = Text_files_path + 'Dynamics/Ensemble/'+parameters_path)
    data, return_data_type = get_data_b(folder_path = Text_files_path + 'Dynamics/Ensemble/L%d/'%L+parameters_path, data_type = 'ranking_1_affinity_N0')
    
    n_first_clones = 50
    if(return_data_type):
        final_E = data[0]
        final_E_log = data[1]
        counts_final_E = data[2]
        trajectories = data[3]
        trajectories_rank = data[4]
    else:        
        print('Please run "simulation_ensemble_ranking_affinity.py" script first ...') 

    final_E_log = np.exp(final_E_log/counts_final_E)
    ranking = np.arange(1, n_first_clones+1)
    fit = ranking**(1/(beta_act))
    popt, pcov = curve_fit(my_linear_func, np.log(ranking[~np.isnan(final_E_log)][:50]), np.log(final_E_log[~np.isnan(final_E_log)][:50]))

    exponent_sim2.append(-popt[1])
    std_exponent_sim2.append(np.sqrt(pcov[1,1]))

    # ----- THEORY -------

    kappas_theory = np.linspace(1, 4.1, 400)
    exponent_theory = []
    exponent_theory2 = []

    for kappa in kappas_theory:
        if(kappa<beta_r):
            beta_kappa, E_kappa, Kd_kappa = get_p_properties(betas, Q0, Es, dE, kappa)
            beta_act = np.min([beta_r, beta_kappa])
            exponent_theory.append(kappa*lambda_B/(lambda_A*beta_act))
            exponent_theory2.append(1/(beta_act))
        else:
            beta_kappa, E_kappa, Kd_kappa = get_p_properties(betas, Q0, Es, dE, kappa)
            beta_act = np.min([beta_r, beta_kappa])
            exponent_theory.append(kappa*lambda_B/(lambda_A*beta_act*1))
            exponent_theory2.append(1/(beta_act*1))

    ax_exponents.plot(kappas_theory, exponent_theory, color = 'indigo', linestyle = '-', marker = '', linewidth = 3, ms = 14, alpha = 1)
    if(i_N_r==0):
        ax_exponents.plot(p, exponent_sim[-1], color = 'indigo', linestyle = '', marker = 'D', linewidth = 3, ms = 14, alpha = 1, label = r'$\zeta$')
    else:
        ax_exponents.plot(p, exponent_sim[-1], color = 'indigo', linestyle = '', marker = 'D', linewidth = 3, ms = 14, alpha = 1)
    ax_exponents.errorbar(x=p, y=exponent_sim[-1], yerr = std_exponent_sim[-1], ls = 'none', color = 'indigo', alpha = .6)
    #ax_exponents.vlines(beta_r, .33, .65, lw = 1, ls = '--', color = 'black')

    #ax_exponents.scatter(kappas_theory[np.array(exponent_theory)<0.60][-1], 0.60, s = 140, facecolors='none', edgecolors='indigo', marker = 'o', lw=2)
    #ax_exponents.errorbar(x = kappas_theory[np.array(exponent_theory)<0.60][-1], y = 0.60, yerr = 0.03, xerr = np.array([[kappas_theory[np.array(exponent_theory)<0.60][-1] - kappas_theory[np.array(exponent_theory)<(0.60-0.03)][-1], kappas_theory[np.array(exponent_theory)<(0.60+0.03)][-1] - kappas_theory[np.array(exponent_theory)<0.60][-1]]]).T, color = 'indigo', alpha = .6)
    print(kappas_theory[np.array(exponent_theory)<0.60][-1], kappas_theory[np.array(exponent_theory)<(0.60-0.03)][-1], kappas_theory[np.array(exponent_theory)<(0.60+0.03)][-1])

    #ax_exponents.scatter(kappas_theory[np.array(exponent_theory)<0.57][-1], 0.57, s = 140, facecolors='none', edgecolors='indigo', marker = '^', lw=2)
    #ax_exponents.errorbar(x = kappas_theory[np.array(exponent_theory)<0.57][-1], y = 0.57, yerr = 0.12, ls = 'none', color = 'indigo', alpha = .6)

    #ax_exponents_2 = ax_exponents.twinx()
    ax_exponents.plot(kappas_theory, exponent_theory2, color = 'olive', linestyle = '-', marker = '', linewidth = 3, ms = 14, alpha = 1)
    if(i_N_r==1):
        ax_exponents.plot(p, exponent_sim2[-1], color = 'olive', linestyle = '', marker = 'D', linewidth = 3, ms = 14, alpha = 1, label = r'$1/\beta^*$')
    else:
        ax_exponents.plot(p, exponent_sim2[-1], color = 'olive', linestyle = '', marker = 'D', linewidth = 3, ms = 14, alpha = .8)
    ax_exponents.errorbar(x=p, y=exponent_sim2[-1], yerr = std_exponent_sim2[-1], ls = 'none', color = 'olive', alpha = .6)

my_plot_layout(ax = ax_exponents, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_exponents.legend(fontsize = 28, title_fontsize = 30)
#ax_exponents.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#ax_exponents.set_ylim(bottom = 2e-2)
#ax_exponents.set_yticks([1, 0.1, 0.01, 0.001])
#ax_exponents.set_yticklabels([1, 0.1, 0.01])

#my_plot_layout(ax = ax_exponents_2, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_exponents.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
#ax_exponents.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#ax_exponents.set_ylim(bottom = 2e-2)
#ax_exponents.set_yticks([1, 0.1, 0.01, 0.001])
#ax_exponents.set_yticklabels([1, 0.1, 0.01])


fig_exponents.savefig('../../Figures/1_Dynamics/Ensemble/L%d/Exponent_N0_'%L+energy_model+'.pdf')

print('----END-----')




