import sys
sys.path.append('../../my_lib/')
from funcs import*
plt.rcParams['text.usetex'] = True

text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

project = 'memory_response'
subproject = 'growth_Th'

output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
os.makedirs(output_plot, exist_ok=True)

#----------------------------------------------------------------
model = 'TCRen'

# --- parameters ---
U_c     = 10.0    # activation threshold
lambda_A = 1.     # exponential antigen growth rate
dt       = 1e-4    # time step

k0_list  = [0.01, 1.0]  # slow vs fast antigen accumulation per cell
dlambda_list = np.linspace(0.6, 1.4, 5)  # ratio of cell division rate to antigen growth rate
print(dlambda_list)
def simulate(k0, tau_B = 1, seed=0):
    p_div = dt / tau_B  # division probability for activated cells
    if lambda_A == 0:
        T_max = tau_B * 20  # total simulation time
    else:
        T_max = tau_B * 10  # total simulation time
    np.random.seed(int(seed))
    # each cell: [U, activated_flag]
    cells = [[0.0, False]]
    times, sizes, antigens = [], [], []
    t = 0.0
    while t < T_max:
        u_act = k0 * np.exp(lambda_A * t)  # per-cell accumulation rate at time t        
        new_cells = []
        for U, act in cells:
            # 1) antigen accumulation
            U += u_act * dt
            if not act:
                # 2) activation once threshold is reached
                if U >= U_c:
                    act = True
                    new_cells.append([U, act]) 
                else:
                    new_cells.append([U, act]) 
            else:
                # 3) division + dilution
                if np.random.rand() < p_div:
                    xi = np.random.rand()      # random split fraction
                    U1 = xi * U
                    U2 = U - U1
                    new_cells.append([U1, U1 > U_c])
                    new_cells.append([U2, U2 > U_c])
                else:
                    new_cells.append([U, act])                

        cells = new_cells
        times.append(t)
        sizes.append(len(cells))
        antigens.append(np.mean([U for U, act in cells]))
        t += dt

    return np.array(times), np.array(sizes), np.array(antigens)

# --- run for two regimes and plot ---
for dlambda in tqdm(dlambda_list):
    fig_cells, ax_cells = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_antigen, ax_antigen = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

    if lambda_A == 0:
        lambda_B = dlambda * 1   # cell division rate
        tau_B = 1/lambda_B     # mean division time once above threshold
    else:
        lambda_B = dlambda*lambda_A   # cell division rate
        tau_B = 1/lambda_B     # mean division time once above threshold
    
    for k0 in k0_list:
        t, B, A = simulate(k0, tau_B=tau_B, seed=datetime.now().timestamp())
        ax_cells.plot(t/tau_B, B, label=f"k0={k0}")
        ax_antigen.plot(t/tau_B, A, label=f"k0={k0}")

    ax_cells.plot(t/tau_B, np.exp((lambda_B)*(np.array(t)-1.7)), '--', color = 'k', label=r'$\lambda_B$')
    ax_cells.plot(t/tau_B, np.exp((lambda_A)*(np.array(t)-1.7)), '-', color = 'k', label=r'$\lambda_A$')

    if lambda_A > 0:
        ax_antigen.plot(t/tau_B, np.exp((lambda_B)*(np.array(t)-1.7)), '--', color = 'k', label=r'$\lambda_B$' )
        ax_antigen.plot(t/tau_B, np.exp((lambda_A - lambda_B)*(np.array(t)-1.7)), ':', color = 'k', label=r'$\lambda_A - \lambda_B$')
        ax_antigen.plot(t/tau_B, np.exp((lambda_A)*(np.array(t)-1.7)), '-', color = 'k', label=r'$\lambda_A$')

    ax_cells.set_yscale('log')
    ax_cells.set_xlabel(r'Time/$\tau_B$', fontsize=16)
    ax_cells.set_ylabel('Cell population size', fontsize=16)
    ax_cells.set_ylim(1, 1e4)
    ax_cells.legend(fontsize=14)
    fig_cells.savefig(output_plot + '/cells_tauB-%.2f_lamA-%.1f.pdf'%(tau_B, lambda_A))

    if lambda_A == 0:
        ax_antigen.set_yscale('linear')
    else:
        ax_antigen.set_yscale('log')
        ax_antigen.set_ylim(1, 1e4)
    ax_antigen.set_xlabel(r'Time/$\tau_B$', fontsize=16)
    ax_antigen.set_ylabel('Presented antigen', fontsize=16)
    ax_antigen.legend(fontsize=14)
    fig_antigen.savefig(output_plot + '/antigen_tauB-%.2f_lamA-%.1f.pdf'%(tau_B, lambda_A))