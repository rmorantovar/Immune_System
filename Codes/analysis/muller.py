import sys
sys.path.append('../library/')
#from Immuno_models import*
from functions import*
plt.rcParams['text.usetex'] = True


Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

Matrix = 'MJ2'

if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()

#--------------- PARAMETERS ---------------------
N_ens = 1
N_r = 5e4
N_r = 5e5
N_r = 1e7
T0 = 3
Tf = 8
#Tf = 10
dT = 0.01
lambda_A = 6
#lambda_A = 8
k_pr = 0.1
#k_pr = 180 # hour^-1
k_pr = k_pr*24 #days^-1

thetas = [2.2, 2.0, 1.8, 1.5]#, 1]
thetas = [1.4, 1.8, 2.2]
thetas = [2, 1]

colors_theta = ['tab:blue','darkblue', 'olive', 'orange', 'darkred']
colors_theta = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
lambda_B = lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
#N_c = 1e5
E_ms = -27.63
C = 3e4

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'FTSENAYCGR'
antigen = 'TACNSEYPNTTK'
antigen = 'TACNSEYPNTTKCGRWYC'
#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'

L=len(antigen)

#----------------------------------------------------------------

#--------------------------Energy Motif--------------------------
PWM_data = get_motif(antigen, 'TCRen', Text_files_path)
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)


Es, dE, Q0, betas = calculate_Q0(0.01, 50, 200000, PWM_data, E_ms, L)
S = np.log(Q0)
Kds = np.exp(Es[:-1])

beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
print('beta_r = %.1f'%beta_r)
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
Kd_r = np.exp(E_r)

E_pr = Es[:-1][Kds<(k_pr/k_on)][-1]
Kd_pr = np.exp(E_pr)
beta_pr = betas[:-1][Kds<Kd_pr][-1]
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
#----------------------------------------------------------------

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
colors = ['tab:blue', 'tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = [0]
linear = 0

for i_theta, theta in enumerate(thetas):
    print('theta = %.2f...'%theta)
    for rep in range(6):
        fig_clones2, ax_clones2 = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
        beta_theta, E_theta, Kd_theta = get_theta_properties(betas, Q0, Es, dE, theta)

        parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, theta, linear, N_ens)+energy_model
        data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies%d.txt'%rep, sep = '\t', header=None)

        min_e_data = np.min(data[0])
        max_e_data = np.max(data[0])

        data_active = data.loc[data[1]==1]
        min_act_t = np.min(data_active[3])
        print('t_act_data: %.2f'%min_act_t)
        data_active = data_active.loc[data_active[3]<(min_act_t+1.3)]
        activation_times = np.array(data_active[3])
        print(len(activation_times))
        energies  = np.array(data_active[0])

        ar1, ar2 = np.histogram(activation_times, bins = time)

        data_N_active_linages = np.cumsum(ar1)

        #---- B cell linages ----
        clone_sizes = get_clones_sizes_C(int(data_N_active_linages[-1]), time, activation_times, lambda_B, C, dT, )
        print('Activated clones:',data_N_active_linages[-1], np.shape(clone_sizes))
        #---- Activation rate ----
        #-----------------------------Activation time------------------------
        t_act = get_t_act(time, N_r, Q0, k_on, k_pr, lambda_A, Es, dE, theta, N_c)
        #---------------------------------------------------------------- 
        #-----------------------------m(t)-----------------------------------
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*4)/N_A, Es, theta, lambda_A, N_c, dE)
        M_r = N_r*N_c*np.sum(Q0*p_a*dE)
        print(M_r)
        #----------------------------------------------------------------

        #---- Total Pop size ----
        total_pop = np.sum(clone_sizes, axis = 0)
        total_pop_active = total_pop - total_pop[0] + 1
        t_C = time[total_pop_active<(C-1)][-1] # Calculate time for reaching carrying capacity

        #-----t_C filter-------
        lim_size = 2
        clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
        print('Applying filter...')
        print('Activated clones:', np.shape(clone_sizes_C))
        #-------------------------

        total_pop_active = np.sum(clone_sizes_C, axis = 0) #C
        bcell_freqs = clone_sizes_C/total_pop_active
        entropy = -np.sum(bcell_freqs*np.log(bcell_freqs), axis = 0)

        #---- Stackplots ----
        colors_muller = []
        min_bell_freq = np.min(bcell_freqs[:,-1])
        for c in range(int(len(clone_sizes_C[:,0]))):
            if bcell_freqs[c, -1]>(20*min_bell_freq/theta):
                colors_muller.append(colors_theta[i_theta])
            else:
                colors_muller.append('silver')

        angles = np.random.random(len(energies_C))*2*np.pi #C
        energies_C = energies_C - (E_ms) 
        radious = ((energies_C/(30))) #C

        ax_clones2.stackplot(time, bcell_freqs, colors = colors_muller);

        my_plot_layout(ax = ax_clones2, ticks_labelsize=34)

        cumsum_freqs = np.cumsum(bcell_freqs, axis = 0)

        if(i_theta==2):
            for j in range(len(activation_times_C)):
                ax_clones2.plot(time, cumsum_freqs[j, :], linewidth = .5, color = 'black')

        ax_clones2.set_xlim(T0, Tf)
        ax_clones2.set_ylim(0, 1)
        ax_clones2.set_yticks([])
            
        fig_clones2.savefig('../../Figures/1_Dynamics/Trajectories/Muller/B_cell_clones_2_theta-%.2f_%d.pdf'%(theta, rep), dpi = 10)

