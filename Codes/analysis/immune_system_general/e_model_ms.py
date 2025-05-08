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
sigmas = np.linspace(1, 1.8, 10)
ls = np.arange(15, 26, 1)
energy_models = ['MJ2', 'TCRen']#, 'Gaussian']
Alphabets = ['MJ2', 'TCRen']#, 'TCRen']
colors_models = [my_blue, my_red]

for i_e, energy_model in enumerate(energy_models):
    print('Model = ', energy_model)
    fig, ax = plt.subplots(figsize = (10, 8))
    K_star = dict()
    for i_sigma, sigma in enumerate(sigmas):
        print('sigma = ', sigma)
        K_star[sigma] = dict()
        Motif, M, Alphabet = get_motif('R', energy_model, Text_files_path)
        print(np.std(M))
        M/=np.std(M)
        M*=sigma
        for i_l, l in enumerate(ls):   
            print('l = ', l)    
            E_star = 0
            for a in tqdm(range(n_ens)):
                # --- Antigen ---
                antigen_ids = np.random.randint(0, 20, l) #L=20
                antigen = ''.join([Alphabet[index] for index in antigen_ids])
                #antigen = 'TACNSEYPNTTRAKCGRWYC' 
                #--------------------------Energy Motif--------------------------
                Motif, M, Alphabet = get_motif_em(antigen = antigen, energy_model = energy_model, M = M, read_matrix = False, Text_files_path = Text_files_path)
                E_ms = np.sum([np.min(Motif[:,i]) for i in range(len(Motif[0,:]))])
                E_peak = np.sum([np.mean(Motif[:,i]) for i in range(len(Motif[0,:]))])
                #Change values by the minimum
                for j in np.arange(l):
                    Motif[:,j]= Motif[:,j] - np.min(Motif[:,j], axis=0)
                #--------------------------Entropy function--------------------------
                Es, dE, Q0, betas = calculate_Q0(0.01, 100, 500000, Motif, my_E_ms, l)
                Kds = np.exp(Es[:-1])
                beta_r, E_r, K_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
                E_star+=np.log(K_r)

            E_star/=n_ens
            K_star[sigma][l] = np.log10(np.exp(E_star))
    df = pd.DataFrame(K_star)
    ax2 = sns.heatmap(ax = ax, data = df, center = -7, cmap = 'coolwarm', yticklabels = [r'$%d$'%l for l in ls], xticklabels = [r'$%.2f$'%s for s in sigmas])
    ax2.figure.axes[-1].yaxis.label.set_size(24)
    ax.tick_params(labelsize = 24)
    ax.set_title(r'$K^*$', fontsize = 24)
    #my_plot_layout(ax=ax, xscale = 'linear', yscale = 'linear', title = r'$K^*$')
    #ax.set_xlim(1e-12, 1e1)
    #ax.set_xlim(left = np.exp(my_E_ms-1))
    fig.savefig('../../Figures/0_Shape_Space/K_star_'+energy_model+'.pdf')


