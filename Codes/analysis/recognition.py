import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------
N_ens = 1
N_rs = [1e7, 1e8, 1e9, 1e10]
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
E_ms = -25
C = 1e4
AA = 1

kappas = [2.2, 2.0, 1.8, 1.5]#, 1]
kappas = [1.4, 1.8, 2.2]
kappas = np.linspace(1, 4, 50)

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
color_list = np.array([my_blue2, my_green, my_brown, my_red, my_gold])

#colors_kappa = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_kappa = np.flip(['tab:blue','tab:green','tab:red'])
colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]

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
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, model, Text_files_path)
print('min_e_PWM=%.4f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
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
print('beta_a = %.2f'%beta_pr, 'K_a = %.2e'%Kd_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')
#--------------------------Loops--------------------------
print('Loops...')

fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.1, 'right':.92, 'bottom':.1, 'top': 0.96})
ax_2 = ax.twinx()
K_rec = dict()
beta_rec = dict()
col_map = 'cividis'
Rainbow = plt.cm.get_cmap('rainbow', 50)
Viridis = plt.cm.get_cmap('viridis', 50)
for N_r in N_rs:
    K_rec[N_r] = dict()
    beta_rec[N_r] = dict()
    print('________')
    print('N_r = %.0e'%N_r)
    fig_Q0, ax_Q0 = plt.subplots(figsize=(5,4), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
    #print('beta_r = %.1f'%beta_r, 'K_d_r = %.2e'%Kd_r)

    for i_kappa, kappa in enumerate(kappas):
        #print('--------')

        beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
        #print('beta_kappa = %.2f...'%kappa, 'K_kappa = %.2e...'%Kd_kappa)

        K_rec[N_r][kappa] = np.max([Kd_r, Kd_kappa])
        beta_rec[N_r][kappa] = np.min([beta_r, beta_kappa])

    ax.plot(kappas, K_rec[N_r].values(), color = Viridis(int(50*(1-abs((np.log(N_r)-np.min(np.log(N_rs)))/(np.max(np.log(N_rs))-np.min(np.log(N_rs))))))), ls = '--', lw = 4) 
    ax_2.plot(kappas, beta_rec[N_r].values(), color = Viridis(int(50*(1-abs((np.log(N_r)-np.min(np.log(N_rs)))/(np.max(np.log(N_rs))-np.min(np.log(N_rs))))))), lw = 4, label = r'$10^{%d}$'%(int(np.log10(N_r)))) 

my_plot_layout(ax=ax, yscale = 'log', xscale = 'linear', ticks_labelsize = 38)
#ax.set_xticks([])
#ax.set_xlim(right = 1e1) #use 1e-3 for other plots
#ax.set_ylim(bottom = 1e-13, top = 1.5)

my_plot_layout(ax=ax_2, yscale = 'linear', xscale = 'linear', ticks_labelsize = 38)
ax_2.legend(loc = 6, fontsize = 20)
#ax_2.set_xticks([])
#ax_2.set_xlim(right = 1e1) #use 1e-3 for other plots
#ax_2.set_ylim(bottom = 1e-13, top = 1.5)


fig.savefig('../../Figures/_Summary/recognition.pdf')




