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

# ---------- parameters ----------
u0        = 0.01      # base on-rate
lambda_A  = 2.0       # antigen growth rate
lambda_B = 3 # * np.log(2) #(days)^-1
w         = 1.0       # presentation per binding event
pi_c      = 4.0      # activation threshold
tau_B     = 1/lambda_B       # mean division time once activated
T_max     = 4.      # total simulation time
dt        = 1e-3      # time step
T0 = 0
Tf = 8
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1

k_on = 1e6*24*3600; #(M*days)^-1
b0 = 1e5*10
E_ms = -24
C = 1e4
p = 3

m_c = int(np.ceil(pi_c / w))  # min binding events to activate

rng = np.random.default_rng(time.time_ns() % (2**32 - 1))

for w in [0.5, 1.0, 2.0]:
    m_c = int(np.ceil(pi_c / w))  # min binding events to activate
    # ---------- state ----------
    # each cell: [N_bind, activated_flag]
    cells = [ [0, False] ]   # start from 1 naïve B cell

    times = []
    sizes = []
    antigens = []

    t = 0.0
    while t < T_max:
        u_on = u0 * 100*np.exp(lambda_A * t)      # time-dependent association rate
        p_div = dt / tau_B                    # division probability for activated cells
        next_cells = []
        for N, act in cells:
            # binding for all cells
            N += rng.poisson(u_on * dt)
            antigens.append(N)
            act = (N >= m_c)  # activation by threshold

            if act < p_div:
                # division with stochastic partition of complexes
                x1 = rng.random()
                N1 = int(np.round(x1 * N))
                N2 = N - N1
                next_cells.append([N1, N1 >= m_c])
                next_cells.append([N2, N2 >= m_c])
            else:
                next_cells.append([N, act])

        cells = next_cells
        times.append(t)
        sizes.append(len(cells))
        t += dt

    # example: print some output
    ax_cells.plot(times, sizes, label = 'w = '+str(w))
    ax_antigen.plot(times, antigens, label = 'w = '+str(w))


ax_cells.plot(times, np.exp(lambda_B*(np.array(times)-1.7)), '--', label='Exponential growth, rate='+str(lambda_B), color = 'k' )

ax_cells.set_yscale('log')
ax_cells.set_xlabel('Time (days)', fontsize=16)
ax_cells.set_ylabel('Cell population size', fontsize=16)
ax_cells.set_ylim(1, 1e6)
ax_cells.legend(fontsize=14)
fig_cells.savefig(output_plot + '/cells_tauB-%.2f_tauA-%.2f.pdf'%(tau_B, 1/lambda_A))

ax_antigen.set_yscale('log')
ax_antigen.set_xlabel('Time (days)', fontsize=16)
ax_antigen.set_ylabel('Antigen', fontsize=16)
ax_antigen.set_ylim(1, 1e6)
ax_antigen.legend(fontsize=14)
fig_antigen.savefig(output_plot + '/antigen_tauB-%.2f_tauA-%.2f.pdf'%(tau_B, 1/lambda_A))
