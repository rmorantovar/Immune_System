import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
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
#----------------------------------------------------------------

energy_models = ['MJ2', 'TCRen']
for energy_model in energy_models:
    print('Model = ', energy_model)
    Alphabet = np.loadtxt(Text_files_path+'Input_files/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)
    print(Alphabet)
    print('--------')
    print('Loops...')
    #--------------------------Loops--------------------------
    fig_E_ms, ax_E_ms = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_peak, ax_peak = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    fig_K_r, ax_K_r = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
    K_rs = []
    Q0_avg = np.zeros(399998)
    n_ens = 100
    for i in tqdm(range(n_ens)):
        antigen_ids = np.random.randint(0, 20, l) #L=20
        antigen = ''.join([Alphabet[index] for index in antigen_ids])

        #--------------------------Energy Motif--------------------------
        PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path)
        E_ms = np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
        E_peak = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
        #print('min_E_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
        #print('min_K_PWM=%.2e'%np.exp(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
        #Change values by the minimum
        for i_l in np.arange(l):
            PWM_data[:,i_l]-=np.min(PWM_data[:,i_l], axis=0)
        #--------------------------Entropy function--------------------------
        Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, l)
        Q0_avg+=np.log(Q0)
        Kds = np.exp(Es[:-1])
        beta_r, E_r, K_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
        K_rs.append(K_r*np.exp(-E_ms+my_E_ms))
        if(i%5==0):
            ax_E_ms.plot(np.exp(Es[:-1]-E_ms+my_E_ms), Q0, color = 'grey', lw = .4, alpha = .2)
            ax_peak.plot(np.exp(Es[:-1]-E_peak+0), Q0, color = 'grey', lw = .4, alpha = .2)

    Q0_avg/=n_ens
    Q0_avg = np.exp(Q0_avg)
    ax_E_ms.plot(np.exp(Es[:-1]-E_ms+my_E_ms), Q0_avg, color = 'black', lw = 2, ls = '-')
    my_antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
    #--------------------------Energy Motif--------------------------
    PWM_data, M, Alphabet = get_motif(my_antigen, energy_model, Text_files_path)
    E_ms = np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    E_peak = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    #print('min_E_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
    #print('min_K_PWM=%.2e'%np.exp(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
    #Change values by the minimum
    for i_l in np.arange(l):
        PWM_data[:,i_l]-=np.min(PWM_data[:,i_l], axis=0)
    #--------------------------Entropy function--------------------------
    Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, l)
    Kds = np.exp(Es[:-1])
    beta_r, E_r, K_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
    K_rs.append(K_r*np.exp(-E_ms+my_E_ms))
    ax_E_ms.plot(np.exp(Es[:-1]-E_ms+my_E_ms), Q0, color = 'black', lw = 2, ls = '--')
    ax_peak.plot(np.exp(Es[:-1]-E_peak+0), Q0, color = 'black', lw = 2, ls = '--')

    bins = np.logspace(np.log10(np.exp(my_E_ms)), np.log10(np.exp(my_E_ms+15)), 30)
    data_K_r = ax_K_r.hist(K_rs, bins = bins, color = 'grey')
    ax_K_r.vlines(K_r*np.exp(-Es[0]+my_E_ms), 0, np.max(data_K_r[0]), color = 'black')

    my_plot_layout(ax=ax_E_ms, xscale = 'log', yscale = 'log')
    ax_E_ms.set_xlim(left = np.exp(my_E_ms-1))
    fig_E_ms.savefig('../../Figures/0_Shape_Space/DoS_E_ms_' + energy_model + '.pdf')

    my_plot_layout(ax=ax_peak, xscale = 'log', yscale = 'linear')
    ax_peak.set_xlim(right = np.exp(0+1))
    fig_peak.savefig('../../Figures/0_Shape_Space/DoS_peak_' + energy_model + '.pdf')

    my_plot_layout(ax=ax_K_r, xscale = 'log')
    ax_K_r.set_xlim(left = np.exp(my_E_ms-1))
    fig_K_r.savefig('../../Figures/0_Shape_Space/K_r_' + energy_model + '.pdf')




