import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
n_ens = 200
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
my_E_ms = -25
C = 1e4
N_r = 1e8

print('min_E_PWM=%.2f'%(-25))
print('min_K_PWM=%.2e'%np.exp(-25))

time = np.linspace(T0, Tf, int((Tf-T0)/dT))

color_list = np.array([my_red, my_blue])

l=20
print('--------')
print('l=%d'%(l))
print('--------')
#----------------------------------------------------------------

energy_models = ['MJ2', 'TCRen']#, 'Gaussian']
Alphabets = ['MJ2', 'TCRen']#, 'TCRen']
colors_models = [my_blue, my_red]
fig, ax = plt.subplots(figsize = (10, 8))
fig2, ax2 = plt.subplots(figsize = (10, 8))

#fig_em, ax_em = plt.subplots(1, 2, figsize = (20, 8))
for i_e, energy_model in enumerate(energy_models):

    print('Model = ', energy_model)
    Motif, M, Alphabet = get_motif('R', energy_model, Text_files_path)
    print(Alphabet)
    #M-=np.mean(M)
    print(np.std(M))
    #M/=np.std(M)
    #M*=1.2
    #mean_e = 0
    #var_e = 0
    #for i in range(20):
    #	mean_e+=np.mean(M[i, :])
    #	var_e+=np.var(M[i, :])
    #print(mean_e/20, np.mean(M))
    
    Q0_avg = np.zeros(499998)
    Krs = []
    for a in tqdm(range(n_ens)):
        # --- Antigen ---
        antigen_ids = np.random.randint(0, 20, l) #L=20
        antigen = ''.join([Alphabet[index] for index in antigen_ids])
        #antigen = 'TACNSEYPNTTRAKCGRWYC' 
        ####
        #print(antigen)
        ####
        #--------------------------Energy Motif--------------------------
        Motif, M, Alphabet = get_motif_em(antigen = antigen, energy_model = energy_model, M = M, read_matrix = False, Text_files_path = Text_files_path)
        ####
        #print('min_E_PWM=%.2f'%(np.sum([np.min(Motif[:,i]) for i in range(len(Motif[0,:]))])))
        ####
        #Motif-=np.mean(Motif)
        #Motif/=np.std(Motif)
        ####
        #print('var motif = ',np.var(Motif))
        ####
        E_ms = np.sum([np.min(Motif[:,i]) for i in range(len(Motif[0,:]))])
        E_peak = np.sum([np.mean(Motif[:,i]) for i in range(len(Motif[0,:]))])
        ####
        #print('E_peak = ', E_peak)
        ####
        #print('min_E_PWM=%.2f'%(np.sum([np.min(Motif[:,i]) for i in range(len(Motif[0,:]))])))
        #print('min_K_PWM=%.2e'%np.exp(np.sum([np.min(Motif[:,i]) for i in range(len(Motif[0,:]))])))
        #Change values by the minimum
        mean_e = 0
        var_e = 0
        Motif_new = np.ones_like(Motif)
        for i_l in np.arange(l):
            #ax_em[i_e].scatter(i_l*np.ones_like(Motif[:,i_l]), Motif[:,i_l], color = colors_models[i_e], alpha = .5, marker = '.')
            Motif_new[:,i_l]= Motif[:,i_l] - np.min(Motif[:,i_l], axis=0)
            #ax_em[i_e].scatter(i_l*np.ones_like(Motif_new[:,i_l]), Motif_new[:,i_l], color = colors_models[i_e])
            mean_e+=np.mean(Motif[:, i_l])
            var_e+=np.var(Motif_new[:, i_l])
        #Motif_new/=np.std(Motif_new)
        ####
        #print('var motif 2 = ',var_e/l)
        ####
        #ax_em[i_e].plot(range(l), [np.sort(Motif[:,i_l])[0] for i_l in range(l)], color = colors_models[i_e])
        #ax_em[i_e].plot(range(l), [np.sort(Motif_new[:,i_l])[-1] for i_l in range(l)], color = colors_models[i_e])
        #ax_em[i_e].plot(range(l), [np.sort(Motif_new[:,i_l])[-2] for i_l in range(l)], color = colors_models[i_e], ls = '--')
        #ax_em[i_e].plot(range(l), [np.sort(Motif_new[:,i_l])[-8] for i_l in range(l)], color = colors_models[i_e], ls = '--')
        ####
        #print('Variances = ', np.var(Motif)*l, var_e, np.var(Motif_new)*l)
        ####
        #--------------------------Entropy function--------------------------
        Es, dE, Q0, betas = calculate_Q0(0.01, 100, 500000, Motif_new, my_E_ms, l)
        Q0_avg+=np.log(Q0)
        Kds = np.exp(Es[:-1])
        beta_r, E_r, K_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
        Krs.append(K_r)
        if(a%10==0):
            ax.plot(np.exp(Es[:-1]+0), Q0, lw = 2, alpha = .2, color = colors_models[i_e])
        #print('Loops...')
        ####
        #print('--------')
        ####
        #ax.plot(Es-E_peak+0, np.exp(-((Es-E_peak+0)/(np.std(Motif_new)*np.sqrt(l)))**2/2)/(np.sqrt(2*np.pi)*np.sqrt(l)*np.std(Motif_new)), color = 'k', ls = '--', alpha = .6)
        #ax.plot(Es-E_peak+0, np.exp(-((Es-E_peak+0)**2/(2*var_e)))/np.sqrt(2*np.pi*var_e), color = colors_models[i_e], ls = '--', alpha = 1)

    Q0_avg/=n_ens
    Q0_avg = np.exp(Q0_avg)
    ax.plot(np.exp(Es[:-1]+0), Q0_avg, lw = 2, alpha = 1, color = colors_models[i_e], label = energy_model)
    ax2.hist((Krs), bins = np.logspace(-11, 3, 80), color = colors_models[i_e], density = False, alpha = 1, label = energy_model)
    #ax.plot(np.exp(data_K_r[1][:-1]), data_K_r[0], color = colors_models[i_e])
    #--------------------------Loops--------------------------
    #my_plot_layout(ax=ax_em[i_e], xscale = 'linear', yscale = 'linear', title = energy_model)
    #ax_em[i_e].set_ylim(-7, 6)

my_plot_layout(ax=ax, xscale = 'log', yscale = 'linear')
ax.legend(fontsize = 20)
ax.set_xlim(1e-12, 1e1)
#ax.set_xlim(left = np.exp(my_E_ms-1))
fig.savefig('../../Figures/0_Shape_Space/energy_model.pdf')

my_plot_layout(ax=ax, xscale = 'log', yscale = 'log')
ax.legend(fontsize = 20)
ax.set_xlim(1e-12, 1e1)
#ax.set_xlim(left = np.exp(my_E_ms-1))
fig.savefig('../../Figures/0_Shape_Space/energy_model_log.pdf')

my_plot_layout(ax=ax2, xscale = 'log', yscale = 'linear')
ax2.set_xlim(1e-12, 2)
ax2.legend(fontsize = 20)
#ax.set_xlim(left = np.exp(my_E_ms-1))
fig2.savefig('../../Figures/0_Shape_Space/K_r.pdf')


#fig_em.savefig('../../Figures/0_Shape_Space/motifs.pdf')

