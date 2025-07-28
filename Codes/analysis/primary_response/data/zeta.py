import sys
sys.path.append('../../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------
N_ens = 100
L_0 = 1e8
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

ps = [1, 3]
ps = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_red, my_blue2, my_green, my_gold, my_brown])
color_list = np.array([my_blue2, my_blue2, my_green, my_green, my_red, my_red, my_gold, my_gold, my_gold])

colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

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
print('beta_a = %.2f'%beta_pr, 'K_a = %.2e'%Kd_pr)

#t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
print('Loops...')

fig_exponents, ax_exponents = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

#--------------------------Loops--------------------------
for lambda_A in lambda_As:

    exponent_sim = []
    exponent_sim2 = []
    std_exponent_sim = []
    std_exponent_sim2 = []

    for i_p, p in enumerate((ps)):
        print('--------')
        print('p = %.2f...'%p)
        beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
        beta_act = np.min([beta_r, beta_p])

        #-----------------Loading data----------------------------
        file = open('../../../out/primary_immune_response'+'/processed_data_ranking_size_p-%.1f_lamA-%.1f_L0-%.0e.pkl'%(p, lambda_A, L_0), 'rb')
        data = pickle.load(file)
        return_data_type = 1

        n_first_clones = 100

        if(return_data_type):
            final_Nb = data[0]
            counts_final_Nb = data[1]
            trajectories = data[2]
            trajectories_rank = data[3]
        else:        
            print('Run "simulation_ensemble_ranking_size.py" script first ...')
        
        final_Nb = final_Nb/counts_final_Nb
        ranking = np.arange(1, n_first_clones+1)
        popt, pcov = curve_fit(my_linear_func, np.log(ranking[~np.isnan(final_Nb)][:100]), np.log(final_Nb[~np.isnan(final_Nb)][:100]))

        exponent_sim.append(-popt[1])
        std_exponent_sim.append(np.sqrt(pcov[1,1]))

        print('--------')

        beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
        beta_act = np.min([beta_r, beta_p])

        #-----------------Loading data----------------------------    
        file = open('../../../out/primary_immune_response'+'/processed_data_ranking_1_affinity_p-%.1f_lamA-%.1f_L0-%.0e.pkl'%(p, lambda_A, L_0), 'rb')
        data = pickle.load(file)
        return_data_type = 1

        n_first_clones = 100
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
        popt, pcov = curve_fit(my_linear_func, np.log(ranking[~np.isnan(final_E_log)][2:100]), np.log(final_E_log[~np.isnan(final_E_log)][2:100]))

        exponent_sim2.append(-popt[1])
        std_exponent_sim2.append(np.sqrt(pcov[1,1]))

    # ----- THEORY -------

    ps_theory_LS = np.linspace(1, beta_r,100)
    ps_theory_HS = np.linspace(beta_r, 5,100)
    exponent_theory_LS = []
    exponent_theory2_LS = []
    exponent_theory_HS = []
    exponent_theory2_HS = []

    for p in ps_theory_LS:
        beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
        beta_act = np.min([beta_r, beta_p])
        exponent_theory_LS.append(p*lambda_B/(lambda_A*beta_act))
        exponent_theory2_LS.append(0)
            
    for p in ps_theory_HS:
        beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
        beta_act = np.min([beta_r, beta_p])
        exponent_theory_HS.append(p*lambda_B/(lambda_A*beta_act*1))
        exponent_theory2_HS.append(1/(beta_act*1))

    exponent_theory_LS = np.array(exponent_theory_LS)
    exponent_theory_HS = np.array(exponent_theory_HS)
    exponent_sim = np.array(exponent_sim)

    ax_exponents.plot(ps_theory_LS, exponent_theory_LS, color = 'indigo', linestyle = '-', marker = '', linewidth = 3, ms = 14, alpha = 1, label = r'$\zeta$')#\textrm{ (size)}$')
    ax_exponents.plot(ps_theory_HS, exponent_theory_HS, color = 'indigo', linestyle = '-', marker = '', linewidth = 3, ms = 14, alpha = 1)
    ax_exponents.plot(ps, exponent_sim, color = 'indigo', linestyle = '', marker = 'D', linewidth = 3, ms = 14, alpha = 1)
    
    #ax_exponents.errorbar(x=ps, y=exponent_sim*(lambda_A/lambda_B), yerr = std_exponent_sim, ls = 'none', color = 'indigo', alpha = .6)
    #ax_exponents.vlines(beta_r, .33, .65, lw = 1, ls = '--', color = 'black')

    data_mean = 0.5
    data_std = 0.02

    ax_exponents.scatter(ps_theory_HS[np.array(exponent_theory_HS)<data_mean][-1], data_mean, s = 140, facecolors='none', edgecolors='indigo', marker = 'o', lw=2)
    ax_exponents.errorbar(x = ps_theory_HS[np.array(exponent_theory_HS)<data_mean][-1], y = data_mean, yerr = 2*data_std, xerr = np.array([[ps_theory_HS[np.array(exponent_theory_HS)<data_mean][-1] - ps_theory_HS[np.array(exponent_theory_HS)<(data_mean-2*data_std)][-1], ps_theory_HS[np.array(exponent_theory_HS)<(data_mean+2*data_std)][-1] - ps_theory_HS[np.array(exponent_theory_HS)<data_mean][-1]]]).T, color = 'indigo', alpha = .6)
    print(ps_theory_HS[np.array(exponent_theory_HS)<data_mean][-1], ps_theory_HS[np.array(exponent_theory_HS)<(data_mean-2*data_std)][-1], ps_theory_HS[np.array(exponent_theory_HS)<(data_mean+2*data_std)][-1])

    #ax_exponents.scatter(ps_theory[np.array(exponent_theory)<0.57][-1], 0.57, s = 140, facecolors='none', edgecolors='indigo', marker = '^', lw=2)
    #ax_exponents.errorbar(x = ps_theory[np.array(exponent_theory)<0.57][-1], y = 0.57, yerr = 0.12, ls = 'none', color = 'indigo', alpha = .6)

    #ax_exponents_2 = ax_exponents.twinx()
    # ax_exponents.plot(ps_theory_LS, np.array(exponent_theory2_LS), color = 'olive', linestyle = '-', marker = '', linewidth = 3, ms = 14, alpha = 1, label = r'$\beta^*$')#\textrm{ (geometric)}$')
    # ax_exponents.plot(ps_theory_HS, np.array(exponent_theory2_HS), color = 'olive', linestyle = '-', marker = '', linewidth = 3, ms = 14, alpha = 1)
    # ax_exponents.plot(ps, np.array(exponent_sim2), color = 'olive', linestyle = '', marker = 'D', linewidth = 3, ms = 14, alpha = 1)
    #ax_exponents.errorbar(x=ps, y=1/np.np.array(exponent_sim2), yerr = std_exponent_sim2, ls = 'none', color = 'olive', alpha = .6)

my_plot_layout(ax = ax_exponents, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
# ax_exponents.legend(fontsize = 24, title_fontsize = 30, loc = 4)
#ax_exponents.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
#ax_exponents.set_ylim(top = 2.4, bottom = .7)
# ax_exponents.set_ylim(bottom = .1)
#ax_exponents.set_yticks([1, 0.1, 0.01, 0.001])
#ax_exponents.set_yticklabels([1, 0.1, 0.01])

#my_plot_layout(ax = ax_exponents_2, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#ax_exponents.legend(fontsize = 32, title_fontsize = 34, title = r'$p$')
#ax_exponents.set_xlim(left = np.exp(E_m+2), right = np.exp(E_m+29))
#ax_exponents.set_ylim(bottom = 2e-2)
#ax_exponents.set_yticks([1, 0.1, 0.01, 0.001])
#ax_exponents.set_yticklabels([1, 0.1, 0.01])


fig_exponents.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2020/Exponent_'+energy_model+'_L0-%.0e.pdf'%(L_0))
plt.close()
print('----END-----')

