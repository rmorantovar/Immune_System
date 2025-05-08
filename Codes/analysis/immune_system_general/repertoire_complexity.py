import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 1

T0 = 3
Tf = 10
Tf_sim = 6.5
#Tf = 10
dT = 0.01
lambda_A = 6
k_pr = 1 # min^-1
k_pr = k_pr*60 # hour^-1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = np.flip([1, 2, 3, 4])

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

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_red, my_blue2, my_green, my_gold, my_brown])

#colors_kappa = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_kappa = np.flip(['tab:blue','tab:green','tab:red'])
colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append([colors_kappa[i], colors_kappa[i], colors_kappa[i], colors_kappa[i]])

lambda_B = lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5
#N_c = 1e5
E_ms = -27.63
E_ms = -25.32
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
#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
#antigen = 'CMFILVWYAGTSQNEDHRKPFMRTPTRMCW' #L=30
L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
model = 'TCRen'
#model = 'MJ2'
Alphabet = np.loadtxt(Text_files_path+'Input_files/Alphabet_'+model+'.txt', dtype=bytes, delimiter='\t').astype(str)
print(Alphabet)

#--------------------------Loops--------------------------
K_naive = 1e-7
print('Loops...')
kappa = 2.0
fig_L, ax_L = plt.subplots()
fig_L2, ax_L2 = plt.subplots()
Ls = np.arange(10, 30, 2)
N_rs = np.logspace(1, 11, 11)
A = np.zeros((len(Ls), len(N_rs)))
A[0, -1] = 1
K_kappas = []
K_ms = []

for i_L, L in enumerate(tqdm(Ls)):
    K = 0
    K_0 = 0
    for i in range(20):
        #--------------------------Energy Motif--------------------------
        antigen = np.random.choice(Alphabet, L)
        PWM_data, M, Alphabet = get_motif(antigen, model, Text_files_path)
        K_0 += np.exp((np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
        #Change values by the minimum
        for i in np.arange(L):
            PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

        avg_E = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))]) + E_ms
        var_E = np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
        #--------------------------Entropy function--------------------------
        Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
        Kds = np.exp(Es[:-1])
        beta_1, E_1, Kd_1 = get_kappa_properties(betas, Q0, Es, dE, 1)

        #--------------------------Proofreading properties--------------------------
        beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)

        beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
        K += Kd_kappa
    K_0 /=20
    K /= 20
    K_kappas.append(K)
    K_ms.append(K_0)
    if(K<K_naive):
        for i_N_r, N_r in enumerate(N_rs):
            beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
            if(Kd_r>K_naive):
                #print(L, '%.0e'%N_r)
                A[i_L, i_N_r] = 1

ax_L.plot(Ls, K_kappas)
ax_L.set_yscale('log')
fig_L.savefig('../../Figures/0_Shape_Space/L.pdf')

ax_L2.plot(Ls, K_ms)
ax_L2.set_yscale('log')
fig_L2.savefig('../../Figures/0_Shape_Space/ms.pdf')

fig, ax = plt.subplots()
sns.heatmap(A, ax=ax, cmap = 'binary', cbar = True)
ax.set_xticklabels(['%.0e'%N_r for N_r in N_rs])
ax.set_yticklabels(Ls);
fig.savefig('../../Figures/0_Shape_Space/map.pdf')





