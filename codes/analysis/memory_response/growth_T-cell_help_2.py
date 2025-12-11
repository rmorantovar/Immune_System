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

fig_cells, ax_cells = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
fig_antigen, ax_antigen = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

# --- parameters ---
pi_c     = 10.0    # activation threshold
tau_B    = 0.5     # mean division time once above threshold
T_max    = tau_B * 10  # total simulation time
dt       = 1e-3

lambda_A = 1.     # exponential antigen growth rate

k0_list  = [0.01, 1.0]  # slow vs fast antigen accumulation per cell

def simulate(k0, seed=0):
    np.random.seed(seed)
    # each cell: [pi, activated_flag]
    cells = [[0.0, False]]
    times, sizes, antigens = [], [], []

    t = 0.0
    while t < T_max:
        k_t = k0 * np.exp(lambda_A * t)  # per-cell accumulation rate at time t
        p_div = dt / tau_B
        new_cells = []

        for pi, act in cells:
            # 1) antigen accumulation
            pi += k_t * dt

            # 2) activation once threshold is reached
            if not act and pi >= pi_c:
                act = True

            # 3) division + dilution
            if act and np.random.rand() < p_div:
                xi = np.random.rand()      # random split fraction
                pi1 = xi * pi
                pi2 = pi - pi1
                new_cells.append([pi1, pi1 >= pi_c])
                new_cells.append([pi2, pi2 >= pi_c])
            else:
                new_cells.append([pi, act])

        cells = new_cells
        times.append(t)
        sizes.append(len(cells))
        antigens.append(np.mean([pi for pi, act in cells]))
        t += dt

    return np.array(times), np.array(sizes), np.array(antigens)

# --- run for two regimes and plot ---
for k0 in k0_list:
    t, B, A = simulate(k0)
    ax_cells.plot(t/tau_B, B, label=f"k0={k0}")
    ax_antigen.plot(t/tau_B, A, label=f"k0={k0}")

ax_cells.plot(t/tau_B, np.exp((1/tau_B)*(np.array(t)-1.7)), '--', label='Exponential growth, rate=%.2f'%(1/tau_B), color = 'k' )
ax_cells.plot(t/tau_B, np.exp((lambda_A)*(np.array(t)-1.7)), '-', color = 'k' )

ax_antigen.plot(t/tau_B, np.exp((1/tau_B)*(np.array(t)-1.7)), '--', color = 'k' )
ax_antigen.plot(t/tau_B, np.exp((lambda_A - 1/tau_B)*(np.array(t)-1.7)), ':', color = 'k' )
ax_antigen.plot(t/tau_B, np.exp((lambda_A)*(np.array(t)-1.7)), '-', color = 'k' )

ax_cells.set_yscale('log')
ax_cells.set_xlabel(r'Time/$\tau_B$', fontsize=16)
ax_cells.set_ylabel('Cell population size', fontsize=16)
ax_cells.set_ylim(1, 1e4)
ax_cells.legend(fontsize=14)
fig_cells.savefig(output_plot + '/cells_tauB-%.2f_lambdaA-%.2f.pdf'%(tau_B, lambda_A))

ax_antigen.set_yscale('log')
ax_antigen.set_xlabel(r'Time/$\tau_B$', fontsize=16)
ax_antigen.set_ylabel('Cell population size', fontsize=16)
ax_antigen.set_ylim(1, 1e4)
ax_antigen.legend(fontsize=14)
fig_antigen.savefig(output_plot + '/antigen_tauB-%.2f_lambdaA-%.2f.pdf'%(tau_B, lambda_A))