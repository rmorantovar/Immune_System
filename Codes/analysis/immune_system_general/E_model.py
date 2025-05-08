import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
n_ens = 100
T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*10
#N_c = 1e5
#E_m = -27.63
my_E_m = -24
C = 1e4
L_0 = 1e9
E_o = 4

print('min_E_PWM=%.2f'%(my_E_m))
print('min_K_PWM=%.2e'%np.exp(my_E_m))

time = np.linspace(T0, Tf, int((Tf-T0)/dT))

print('--------')
#----------------------------------------------------------------
data_ready = 0
energy_models = ['Gaussian', 'MJ2']# 'TCRen']

#fig_em, ax_em = plt.subplots(1, 2, figsize = (20, 8))
for i_model, energy_model in enumerate(energy_models):
    print('--------')
    print('Model = ', energy_model)
    print('--------')

    Motif, M, Alphabet = get_motif_em(antigen = 'R', energy_model = energy_model, M = [], read_matrix = True, Text_files_path = Text_files_path)

    E_trajectories = np.array([], dtype = object)
    rho_trajectories = np.array([], dtype = object)    
    beta_rs = []
    stds = []
    good_ls = []
    good_alphas = []
    K_rs = []
    K_ms = []
    Es_average = np.zeros(1000000-2)
    rho_average = np.zeros(1000000-2)
    counter = 0
    if energy_model=='TCRen':
        ls = range(14, 30)
    if energy_model=='Gaussian':
        ls = range(14, 22)
    if energy_model=='MJ2':
        ls = range(17, 27)
    for l in tqdm(ls):  

        # Loop over antigens 
        for a in range(n_ens):
            # --- Antigen ---
            antigen_ids = np.random.randint(0, 20, l) #L=20
            #antigen = ''.join([Alphabet[index] for index in antigen_ids])
            #antigen = 'TACNSEYPNTTRAKCGRWYC' 
            #--------------------------Energy Motif--------------------------
            #Motif, M, Alphabet = get_motif_em(antigen = antigen, energy_model = energy_model, M = M, read_matrix = False, Text_files_path = Text_files_path)
            Motif = get_motif_em_ids(antigen_ids = antigen_ids, energy_model = energy_model, M = M, read_matrix = False, Text_files_path = Text_files_path)
            # E_0 = 0
            # E_var = 0
            # E_m = 0
            # for i_l in np.arange(l):
            #     E_0+=np.mean(Motif[:, i_l])
            #     E_m+=np.min(Motif[:, i_l])
            #     E_var+=np.var(Motif[:, i_l])
            E_0 = np.sum(Motif.mean(axis = 0))
            E_m = np.sum(Motif.min(axis = 0))
            E_var = np.sum(Motif.var(axis = 0))
            #print(E_m, E_0)
            for j in range(5):
                r = np.random.rand(1) - 0.5
                alpha = (my_E_m + np.sign(r)*np.log(1+abs(r)) - E_o)/(E_m - E_0)
                #alpha = (my_E_m-E_0)/(E_m - E_0)
                Motif_new = Motif*alpha
                # E_0_prime = 0
                # E_var_prime = 0
                # E_m_prime = 0
                E_0_prime = np.sum(Motif_new.mean(axis = 0))
                E_m_prime = np.sum(Motif_new.min(axis = 0))
                E_var_prime = np.sum(Motif_new.var(axis = 0))
                for i_l in np.arange(l):
                    #E_0_prime+=np.mean(Motif_new[:, i_l])
                    #E_m_prime+=np.min(Motif_new[:, i_l])
                    #E_var_prime+=np.var(Motif_new[:, i_l])
                    Motif_new[:,i_l]= Motif_new[:,i_l] - np.min(Motif_new[:,i_l])
                #print(E_0_prime, E_m_prime, E_m_prime - E_0_prime)
                #if((np.exp(E_m_prime- E_0_prime) > 3.7e-11/1.5) & (np.exp(E_m_prime- E_0_prime) < 3.7e-11*1.5)): #### K_m condition ####
                #--------------------------Entropy function--------------------------
                Es, dE, rho_0, betas = calculate_Q0(0.005, 100, 1000000, Motif_new, E_m_prime - E_0_prime + E_o, l)
                beta_r, E_r, K_r = get_repertoire_properties(betas, rho_0, Es, dE, L_0)
                if((K_r > 8e-8/1.5) & (K_r < 8e-8*1.5)): #### K^* condition ####
                    #print('motif found!')
                    Es_average+=Es[:-1]
                    rho_average+=np.log10(rho_0)
                    counter+=1
                    K_rs.append(K_r)
                    K_ms.append(np.exp(E_m_prime - E_0_prime + E_o))
                    beta_rs.append(beta_r)
                    stds.append(np.sqrt(E_var_prime/l))
                    good_ls.append(l)
                    good_alphas.append(alpha)
                    if(a%10==0):
                        E_trajectories = np.append(E_trajectories, Es[:-1][::100])
                        rho_trajectories = np.append(rho_trajectories, rho_0[::100])

    rho_average/=counter
    Es_average/=counter 

    E_trajectories = np.append(E_trajectories, Es_average[::100])
    rho_trajectories = np.append(rho_trajectories, 10**rho_average[::100])

    df = pd.DataFrame({'K_r': K_rs, 'K_m':K_ms,'beta_r':beta_rs,'alpha':good_alphas,'l':good_ls,'var':stds})
    df.to_csv(Text_files_path + 'E_model/All_data_' + energy_model + '.csv')
    f = open(Text_files_path + 'E_model/trajectories_'+energy_model+'.pkl', 'wb')
    pickle.dump([E_trajectories, rho_trajectories], f, pickle.HIGHEST_PROTOCOL)
    
