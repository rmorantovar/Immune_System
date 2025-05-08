import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
n_ens = 10
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
#E_m = -27.63
my_E_m = -25
C = 1e4
N_r = 1e8

print('min_E_PWM=%.2f'%(-25))
print('min_K_PWM=%.2e'%np.exp(-25))

time = np.linspace(T0, Tf, int((Tf-T0)/dT))

print('--------')
#----------------------------------------------------------------

energy_models = ['TCRen', 'Gaussian', 'MJ2']
colors_models = [my_blue, my_green, my_red]

#fig_em, ax_em = plt.subplots(1, 2, figsize = (20, 8))
for i_model, energy_model in enumerate(energy_models):
    print('--------')
    print('Model = ', energy_model)
    print('--------')

    Motif, M, Alphabet = get_motif_em(antigen = 'R', energy_model = energy_model, M = [], read_matrix = True, Text_files_path = Text_files_path)

    fig, ax = plt.subplots(figsize = (10, 8))
    fig2, ax2 = plt.subplots(figsize = (10, 8))
    fig3, ax3 = plt.subplots(figsize = (10, 8))
    fig4, ax4 = plt.subplots(figsize = (10, 8))
    fig5, ax5 = plt.subplots(figsize = (10, 8))
    fig6, ax6 = plt.subplots(figsize = (10, 8))
    fig7, ax7 = plt.subplots(figsize = (10, 8))

    beta_rs = []
    vars_l = []
    good_ls = []
    good_alphas = []
    K_rs = []
    K_ms = []
    Es_average = np.zeros(1000000-2)
    rho_average = np.zeros(1000000-2)
    counter = 0
    if energy_model=='TCRen':
        ls = range(21, 29)
        alphas = np.linspace(.5, 1.2, 30)
    if energy_model=='Gaussian':
        ls = range(18, 26)
        alphas = np.linspace(.5, 1.2, 30)
    if energy_model=='MJ2':
        ls = range(18, 26)
        alphas = np.linspace(.5, 1.2, 30)
    for l in ls:  
        print('--------')
        print('l=%d'%(l))
        #for alpha in tqdm(alphas):
        #print('--------')
        #print('alpha=%.2f'%(alpha))
        #print(Alphabet)
        #M-=np.mean(M)
        #print(np.std(M))
        #M/=np.std(M)   # NORMALIZING THE ENERGY MATRIX
        #M*=1.1
        
        # Loop over antigens 
        for a in tqdm(range(n_ens)):
            # --- Antigen ---
            antigen_ids = np.random.randint(0, 20, l) #L=20
            antigen = ''.join([Alphabet[index] for index in antigen_ids])
            #antigen = 'TACNSEYPNTTRAKCGRWYC' 
            #--------------------------Energy Motif--------------------------
            Motif, M, Alphabet = get_motif_em(antigen = antigen, energy_model = energy_model, M = M, read_matrix = False, Text_files_path = Text_files_path)
            E_0 = 0
            E_var = 0
            E_m = 0
            for i_l in np.arange(l):
                E_0+=np.mean(Motif[:, i_l])
                E_m+=np.min(Motif[:, i_l])
                E_var+=np.var(Motif[:, i_l])
            #print(E_m, E_0)
            alpha = -24/(E_m - E_0)
            #alpha = (-24-E_0)/(E_m - E_0)
            Motif*=alpha
            E_0_prime = 0
            E_var_prime = 0
            E_m_prime = 0
            for i_l in np.arange(l):
                E_0_prime+=np.mean(Motif[:, i_l])
                E_m_prime+=np.min(Motif[:, i_l])
                E_var_prime+=np.var(Motif[:, i_l])
                Motif[:,i_l]= Motif[:,i_l] - np.min(Motif[:,i_l])
            #print(E_var, E_var_prime)
            #print(E_m_prime - abs(E_0_prime - E_0), E_0_prime)
            #if((np.exp(E_m_prime- E_0_prime) > 3.7e-11/1.5) & (np.exp(E_m_prime- E_0_prime) < 3.7e-11*1.5)): #### K_m condition ####
            #--------------------------Entropy function--------------------------
            Es, dE, rho_0, betas = calculate_Q0(0.005, 100, 1000000, Motif, E_m_prime - E_0_prime, l)
            beta_r, E_r, K_r = get_repertoire_properties(betas, rho_0, Es, dE, N_r)
            if((K_r > 1e-7/1.5) & (K_r < 1e-7*1.5)): #### K^* condition ####
                print('motif found!')
                Es_average+=Es[:-1]
                rho_average+=np.log10(rho_0)
                counter+=1
                K_rs.append(K_r)
                K_ms.append(np.exp(E_m_prime - E_0_prime))
                beta_rs.append(beta_r)
                vars_l.append(E_var_prime/l)
                good_ls.append(l)
                good_alphas.append(alpha)
                ax.plot(np.exp(Es[:-1]), (rho_0), lw = 2, alpha = .05, color = colors_models[i_model])
            #print('Loops...')
            ####
            #print('--------')

    rho_average/=counter
    Es_average/=counter
    print('%d motif were found!'%counter)
    ax.plot(np.exp(Es_average), 10**rho_average, lw = 1.5, alpha = 1, color = colors_models[i_model], label = energy_model)

    data_k_r = np.histogram((K_rs), bins = np.logspace(-11, -4, 50), density = False)
    data_k_m = np.histogram((K_ms), bins = np.logspace(-15, -6, 50), density = False)
    data_beta_r = np.histogram((beta_rs), bins = np.linspace(1.5, 3.5, 20) density = False)
    data_alpha = np.histogram((good_alphas), bins = np.linspace(.5, 1.5, 20), density = False)
    data_ls = np.histogram((good_ls), bins = range(18, 28), density = False)
    data_vars = np.histogram((vars_l), density = False)
    #
    ax2.plot((data_k_r[1][:-1]), data_k_r[0], color = colors_models[i_model], ls = '--', marker = '', label = r'$K^*$')
    ax2.plot((data_k_m[1][:-1]), data_k_m[0], color = colors_models[i_model], ls = '-', marker = '', label = r'$K_m$')
    ax3.plot((data_beta_r[1][:-1]), data_beta_r[0], color = colors_models[i_model], ls = '-', marker = '')
    ax4.plot((data_alpha[1][:-1]), data_alpha[0], color = colors_models[i_model], ls = '-', marker = '')
    ax5.plot((data_ls[1][:-1]), data_ls[0], color = colors_models[i_model], ls = '-', marker = '')
    ax6.plot((data_vars[1][:-1]), data_vars[0], color = colors_models[i_model], ls = '-', marker = '')
    data2d = np.histogram2d(good_ls, good_alphas, bins = [ls, alphas])
    ls_data = data2d[1][:-1] + np.diff(data2d[1])
    alphas_data = data2d[2][:-1] + np.diff(data2d[2])
    df = pd.DataFrame(data = data2d[0], index = ls_data, columns = ['%.2f'%i for i in alphas_data])
    sns.heatmap(df, ax = ax7)

    #--------------------------Loops--------------------------
    ax.scatter(np.exp(np.mean(np.log(K_ms))), (1/(20**l)), marker = 'o', color = colors_models[i_model], s = 28)
    #ax.scatter(np.exp(np.mean(np.log(K_rs))), 10**rho_average[Es_average<(np.mean(np.log(K_rs)))][-1], marker = '*', color = colors_models[i_model], s = 28)

    my_plot_layout(ax=ax, xscale = 'log', yscale = 'log')
    ax.legend(fontsize = 20)
    #ax.set_xlim(left = 1e-13, right = 1e1)
    #ax.set_xlim(left = np.exp(my_E_m-1))
    fig.savefig('../../Figures/0_Shape_Space/energy_model_'+energy_model+'.pdf')

    my_plot_layout(ax=ax, xscale = 'log', yscale = 'linear')
    fig.savefig('../../Figures/0_Shape_Space/energy_model_2_'+energy_model+'.pdf')

    my_plot_layout(ax=ax2, xscale = 'log', yscale = 'linear')
    ax2.set_xlim(left = 1e-15, right = 1e-4)
    ax2.legend(fontsize = 20)
    #ax.set_xlim(left = np.exp(my_E_m-1))
    fig2.savefig('../../Figures/0_Shape_Space/K_r_K_m_'+energy_model+'.pdf')

    my_plot_layout(ax=ax3, xscale = 'linear', yscale = 'linear')
    ax3.set_xlim(left = 1.5, right = 3.5)
    #ax3.legend(fontsize = 20)
    #ax.set_xlim(left = np.exp(my_E_m-1))
    fig3.savefig('../../Figures/0_Shape_Space/beta_r_'+energy_model+'.pdf')

    my_plot_layout(ax=ax4, xscale = 'linear', yscale = 'linear')
    ax4.set_xlim(left = 0.4, right = 1.6)
    #ax4.legend(fontsize = 20)
    #ax4.set_xlim(left = np.exp(my_E_m-1))
    fig4.savefig('../../Figures/0_Shape_Space/alpha_'+energy_model+'.pdf')

    my_plot_layout(ax=ax5, xscale = 'linear', yscale = 'linear')
    ax5.set_xlim(left = 14, right = 27)
    #ax5.legend(fontsize = 20)
    #ax5.set_xlim(left = np.exp(my_E_m-1))
    fig5.savefig('../../Figures/0_Shape_Space/ls_'+energy_model+'.pdf')

    my_plot_layout(ax=ax6, xscale = 'linear', yscale = 'linear')
    ax6.set_xlim(left = 0, right = 3)
    #ax6.legend(fontsize = 20)
    #ax6.set_xlim(left = np.exp(my_E_m-1))
    fig6.savefig('../../Figures/0_Shape_Space/var_l_'+energy_model+'.pdf')

    #my_plot_layout(ax=ax7, xscale = 'linear', yscale = 'linear')
    ax7.tick_params(labelsize = 24)
    ax7.set_xticks(np.arange(0, len(alphas_data), 2) + 0.5)
    ax7.set_yticks(np.arange(0, len(ls_data)) + 0.5)
    ax7.set_xticklabels([r'$%.2f$'%i for i in alphas[::2]])
    ax7.set_yticklabels([r'$%d$'%i for i in ls_data]);
    #ax7.set_xlim(left = 14, right = 26)
    #ax7.set_ylim(bottom = 0.4, top = 1.6)
    #ax7.legend(fontsize = 20)
    #ax7.set_xlim(left = np.exp(my_E_m-1))
    fig7.savefig('../../Figures/0_Shape_Space/ls_alphas_'+energy_model+'.pdf')

    #fig_em.savefig('../../Figures/0_Shape_Space/motifs.pdf')

