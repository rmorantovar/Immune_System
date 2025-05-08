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

color_list = np.array([my_red, my_blue, myred])

l=21
print('--------')
print('l=%d'%(l))
print('--------')
#----------------------------------------------------------------

energy_models = ['TCRen']#, 'MJ2', 'Gaussian']
Alphabets = ['TCRen', 'MJ2']#, 'TCRen']
colors_models = [my_blue, my_green]

#fig_em, ax_em = plt.subplots(1, 2, figsize = (20, 8))
for i_model, energy_model in enumerate(energy_models):
    fig, ax = plt.subplots(figsize = (10, 8))
    fig2, ax2 = plt.subplots(figsize = (10, 8))
    fig3, ax3 = plt.subplots(figsize = (10, 8))

    print('Model = ', energy_model)
    Motif, M, Alphabet = get_motif('R', energy_model, Text_files_path)
    #print(Alphabet)
    #M-=np.mean(M)
    print(np.std(M))
    #M/=np.std(M)   # NORMALIZING THE ENERGY MATRIX
    #M*=1.1
    #mean_e = 0
    #var_e = 0
    #for i in range(20):
    #	mean_e+=np.mean(M[i, :])
    #	var_e+=np.var(M[i, :])
    #print(mean_e/20, np.mean(M))
    K_rs = []
    K_ms = []
    beta_rs = []
    #Es_average = np.linspace(-35, 5, 40*5)
    #de = np.diff(Es_average)[0]
    # rho_0_average = np.zeros_like(Es_average)
    # rho_0_average2 = np.zeros_like(Es_average)
    # rho_0_std = np.zeros_like(Es_average)
    # normalisations = np.zeros_like(Es_average)
    total_es = []
    total_rho_0 = []
    Es_average = np.zeros(1000000-2)
    rho_average = np.zeros(1000000-2)
    for a in tqdm(range(n_ens)):
        # --- Antigen ---
        antigen_ids = np.random.randint(0, 20, l) #L=20
        antigen = ''.join([Alphabet[index] for index in antigen_ids])
        #antigen = 'TACNSEYPNTTRAKCGRWYC' 
        #--------------------------Energy Motif--------------------------
        Motif, M, Alphabet = get_motif_em(antigen = antigen, energy_model = energy_model, M = M, read_matrix = False, Text_files_path = Text_files_path)
        #Change values by the minimum
        mean_e = 0
        var_e = 0
        for i_l in np.arange(l):
            #ax_em[i_model].scatter(i_l*np.ones_like(Motif[:,i_l]), Motif[:,i_l], color = colors_models[i_model], alpha = .5, marker = '.')
            Motif[:,i_l]= Motif[:,i_l] - np.min(Motif[:,i_l], axis=0)
            #ax_em[i_model].scatter(i_l*np.ones_like(Motif[:,i_l]), Motif[:,i_l], color = colors_models[i_model])
            mean_e+=np.mean(Motif[:, i_l])
            var_e+=np.var(Motif[:, i_l])
        Motif/=np.sqrt(var_e/l)
        #Motif*=1.1
        E_m = np.sum([np.min(Motif[:,i]) for i in range(len(Motif[0,:]))])
        E_peak = np.sum([np.mean(Motif[:,i]) for i in range(len(Motif[0,:]))])
        ####
        #--------------------------Entropy function--------------------------
        Es, dE, rho_0, betas = calculate_Q0(0.005, 100, 1000000, Motif, E_m - E_peak, l)
        rho_0/=rho_0[0]
        rho_0*=(1/20**l)
        Es_average+=Es[:-1]
        rho_average+=np.log10(rho_0)
        # Es_inter = np.linspace(np.min(Es), np.max(Es), 300)  
        # de = np.diff(Es_inter)[0]
        # rho_0_inter = np.ones(300)
        # for i_e, e in enumerate(Es_inter[:-1]):
        #     value_rho0 = np.mean(np.log(rho_0[(Es[:-1] > e) & (Es[:-1] < (e + de))]))
        #     rho_0_inter[i_e] = value_rho0
        #     # if(~np.isnan(value_rho0)):
        #     #     rho_0_average[i_e]+=value_rho0
        #     #     rho_0_average2[i_e]+=(value_rho0**2)
        #     #     normalisations[i_e]+=1
        # total_es+= list(Es_inter)
        # total_rho_0+= list(rho_0_inter)
        beta_r, E_r, K_r = get_repertoire_properties(betas, rho_0, Es, dE, N_r)
        K_rs.append(K_r)
        K_ms.append(np.exp(E_m - E_peak))
        beta_rs.append(beta_r)
        if(a%int(n_ens/10)==0):
            ax.plot(np.exp(Es[:-1]), np.log10(rho_0), lw = 2, alpha = .1, color = colors_models[i_model])
        #print('Loops...')
        ####
        #print('--------')
    
    ####  
    Es_average/=n_ens
    rho_average/=n_ens
    # rho_0_average/=normalisations
    # rho_0_average = np.exp(rho_0_average)

    # rho_0_average2/=normalisations
    # rho_0_average2 = np.exp(rho_0_average2)

    # rho_0_std = np.exp(np.sqrt(np.log(rho_0_average2) - np.log(rho_0_average)**2))
    
    #ax.plot(np.exp(Es_average), rho_0_average, lw = 2, alpha = 1, color = colors_models[i_model], label = energy_model)
    #ax.fill_between(np.exp(Es_average), rho_0_average/rho_0_std, rho_0_average*rho_0_std, color = colors_models[i_model], alpha = .4)

    h = np.histogram2d(total_es, total_rho_0, bins = (np.linspace(-30, 0, 100), np.linspace(-60, 1, 100)))
    data = h[0]
    es_array = np.linspace(-30, 0, 100)
    rho_array = np.linspace(-60, 1, 100)
    rho_0_final = [rho_array[np.argmax(data[i,:])] for i in range(len(es_array)-1)]
    es_final = [es_array[np.argmax(data[:,i])] for i in range(len(rho_array)-1)]
    ax.plot(np.exp(Es_average), rho_average, lw = 2, alpha = 1, color = colors_models[i_model], label = energy_model)
    #sns.heatmap(np.log(h[0]), cmap = 'Blues', ax = ax)

    data_k_r = np.histogram((K_rs), bins = np.logspace(-11, -4, 40), density = False)
    data_k_m = np.histogram((K_ms), bins = np.logspace(-15, -6, 40), density = False)
    data_beta_r = np.histogram((beta_rs), density = False)
    #
    ax2.plot((data_k_r[1][:-1]), data_k_r[0], color = colors_models[i_model], ls = '--', marker = '', label = r'$K_m$')
    ax2.plot((data_k_m[1][:-1]), data_k_m[0], color = colors_models[i_model], ls = '-', marker = '', label = r'$K^*$')
    ax3.plot((data_beta_r[1][:-1]), data_beta_r[0], color = colors_models[i_model], ls = '-', marker = '')
    #--------------------------Loops--------------------------
    ax.scatter(np.exp(np.mean(np.log(K_ms))), np.log10(1/(20**l)), marker = 'o', color = colors_models[i_model], s = 25)
    ax.scatter(np.exp(np.mean(np.log(K_rs))), rho_average[Es_average<(np.mean(np.log(K_rs)))][-1], marker = '*', color = colors_models[i_model], s = 25)

    my_plot_layout(ax=ax, xscale = 'log', yscale = 'linear')
    ax.legend(fontsize = 20)
    #ax.set_xlim(left = 1e-13, right = 1e1)
    #ax.set_xlim(left = np.exp(my_E_m-1))
    fig.savefig('../../Figures/0_Shape_Space/energy_model_'+energy_model+'.pdf')

    my_plot_layout(ax=ax2, xscale = 'log', yscale = 'linear')
    ax2.set_xlim(left = 1e-15, right = 1e-4)
    ax2.legend(fontsize = 20)
    #ax.set_xlim(left = np.exp(my_E_m-1))
    fig2.savefig('../../Figures/0_Shape_Space/K_r_K_m_'+energy_model+'.pdf')

    my_plot_layout(ax=ax3, xscale = 'linear', yscale = 'linear')
    #ax3.set_xlim(left = 1e-15, right = 1e-4)
    #ax3.legend(fontsize = 20)
    #ax.set_xlim(left = np.exp(my_E_m-1))
    fig3.savefig('../../Figures/0_Shape_Space/beta_r_'+energy_model+'.pdf')


#fig_em.savefig('../../Figures/0_Shape_Space/motifs.pdf')

