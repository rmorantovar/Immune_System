import sys
sys.path.append('../../lib/')
from funcs import*

def gillespie_poisson(n, mu, rho, c, max_time, event_callback=None):
    """
    Simulate n+1 Poisson processes using the Gillespie algorithm.
    - n processes with rate mu
    - 1 process with rate rho
    Each event triggers the event_callback function, if provided.
    Returns:
        - events: A list of lists, where each sublist contains the event times for one process.
        - all_event_log: A list of tuples (time, process_index) for all events in chronological order.
    """
    # Initialize times and events
    current_time = 0
    E0 = 1
    energies = [0] * n
    pbs = [1/(1+np.exp(i - E0)*c) for i in energies]
    t = [current_time]
    status_t = [np.prod(1-np.array(pbs))]
    t2 = [0]
    infections_t2 = [0]
    
    # events = [[] for _ in range(n + 1)]  # List to store events for each process
    # all_event_log = []  # Global log of all events (time, process_index)
    
    # Define rates for all processes
    rates = [mu] * n + [rho]
    
    while current_time < max_time:
        # Compute total rate (sum of all process rates)
        total_rate = sum(rates)
        if total_rate == 0:
            break
        
        # Sample the time to the next event
        delta_t = np.random.exponential(1 / total_rate)
        current_time += delta_t
        
        if current_time >= max_time:
            break
        
        # Determine which process gets the event
        cumulative_rates = np.cumsum(rates)
        random_value = np.random.uniform(0, total_rate)
        process_index = np.searchsorted(cumulative_rates, random_value)
        if process_index<n:
            energies[process_index]+=1
        else:
            r = np.random.rand()
            if r < status_t[-1]:
                infections_t2.append(infections_t2[-1]+1)
                t2.append(current_time)
                energies = [0] * n
            else:
                infections_t2.append(infections_t2[-1])
                t2.append(current_time)
        
        # Record the event for the selected process
        pbs = [1/(1+np.exp(i - E0)*c) for i in energies]
        status_t.append(np.prod(1-np.array(pbs)))
        t.append(current_time)
        # events[process_index].append(current_time)
        # all_event_log.append((current_time, process_index))
        
        # Trigger the callback function, if provided
        if event_callback is not None:
            event_callback(current_time, process_index)
    
    return t, status_t, t2, infections_t2

# Define a custom callback function
def my_event_callback(event_time, process_index):
    """
    Custom function to run whenever an event occurs.
    """
    # print(f"Event occurred in process {process_index + 1} at time {event_time:.2f}")

def main():
    # Setting up command-line argument parser
    parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
    parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
    parser.add_argument('--N_ens', type=int, default=10, help="Number of times to execute the process.")
    parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
    parser.add_argument('--N_evo', type=int, default = 0)
    parser.add_argument('--N_epi', type=int, default = 1)
    parser.add_argument('--L0', type=int, default=10**8, help="Number of random sequences.")
    parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
    parser.add_argument('--t_lim', type=float, default=8., help="Threshold for activation time.") # Use 8 for L0>1e6
    parser.add_argument('--E_lim', type=float, default=-6., help="Threshold for the sum of entries.") # Use -6 for L0>1e6
    parser.add_argument('--E_m', type=float, default=-24, help="Threshold for the sum of entries.")
    parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
    parser.add_argument('--p', type=float, default=4, help="# steps.")
    parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
    parser.add_argument('--lamA', type=float, default=6., help="Antigen growth rate.")
    parser.add_argument('--lamB', type=float, default=2., help="Antigen growth rate.")
    parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
    parser.add_argument('--random_antigen', type=int, default=0)
    # parser.add_argument('--antigen', type=str, default='TACNSEYPNTTRAKCGRWYR')
    parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR')
    parser.add_argument('--energy_model', type=str, default = 'TCRen')
    parser.add_argument('--seqs', type=int, default = 1)
    parser.add_argument('--one_WT', type=int, default = 1)
    parser.add_argument('--secondary', type=int, default = 0)
    parser.add_argument('--secondary_all', type=int, default = 0)
    parser.add_argument('--pro', type=str, default='epitope_complexity', help="project.")
    parser.add_argument('--subpro', type=str, default='minimal_ID', help="subproject.")
    parser.add_argument('--exp', type=int, default=2, help="experiment.")
    args = parser.parse_args()

    # ------------ PARAMETERS AND INPUTS ------------
    N_ant = args.N_ant
    N_ens = args.N_ens
    N_inf = args.N_inf
    N_evo = args.N_evo
    N_epi = args.N_epi
    L0 = args.L0
    l = args.l
    E_lim = args.E_lim
    t_lim = args.t_lim
    E_m = args.E_m
    if L0>=1e6:
        chunk_size = args.chunk_size
    else:
        chunk_size = args.L0
    p = args.p
    k_step = args.k_step
    lamA = args.lamA
    lamB = args.lamB
    n_jobs = args.n_jobs
    random_antigen = args.random_antigen
    antigen = args.antigen
    energy_model = args.energy_model
    
    seqs = args.seqs
    one_WT = args.one_WT
    secondary = args.secondary
    secondary_all = args.secondary_all

    if N_evo == -1:
        N_evo = 'R'

    dT = 0.002
    C = 1e4
    T = 12
    time_array = np.linspace(0, T, int((T-0)/dT))
    Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)

    project = args.pro
    subproject = args.subpro
    experiment = args.exp
    root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
    pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
    
    start_time = time.time()
    print('Starting simulation ...')

    # Plot the results
    fig, ax = plt.subplots(figsize=(10, 6))
    fig_c, ax_c = plt.subplots(figsize=(10, 6))
    colors = [my_blue, my_red, my_green, my_cyan, my_blue2]
    output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
    os.makedirs(output_plot, exist_ok=True)

    # Parameters
    # n = 2        # Number of processes at rate mu
    # mu = 20       # Rate for the first n processes
    rho = 1     # Rate for the last process
    max_time = 100  # Simulation time
    mus = np.logspace(-1, 0, 2)
    results = defaultdict(list)

    for mu in tqdm(mus):
        for n in range(1, 8):
            for c in ['0', 'n']:
                inf = 0
                for i in range(200):
                    if c == '0':
                        # Simulate the processes with the callback
                        times, energy_t, times_infections, infections = gillespie_poisson(n, mu, rho, 1, max_time, event_callback=my_event_callback)
                    else:
                        times, energy_t, times_infections, infections = gillespie_poisson(n, mu, rho, n, max_time, event_callback=my_event_callback)
                    
                    inf += infections[-1]
                results['mu'].append(mu)
                results['n'].append(n)
                results['c'].append(c)
                results['infections'].append(inf/200)

    results_df = pd.DataFrame(results)
    results_df['normalized_infections'] = results_df.groupby(['mu', 'c'])['infections'].transform(lambda x: x / x.iloc[0])
    results_df['normalized_infections2'] = results_df.groupby(['mu', 'n'])['normalized_infections'].transform(lambda x: x - x.iloc[0])
    for c in ['0', 'n']:
        results_c = results_df.loc[results_df['c']==c]
        if c == '0':
            for j, mu in enumerate(mus):
                results_mu = results_c.loc[results_c['mu']==mu]
                ax.plot(results_mu['n'], results_mu['normalized_infections'], color = colors[j], ls = '--', linewidth = 2)
        else:
            for j, mu in enumerate(mus):
                results_mu = results_c.loc[results_c['mu']==mu]
                ax.plot(results_mu['n'], results_mu['normalized_infections'], color = colors[j], ls = '-', linewidth = 2, label = r'$10^{%d}$'%np.log10(mu))
                ax_c.plot(results_mu['n'], results_mu['normalized_infections2'], color = colors[j], ls = '-', linewidth = 2, label = r'$10^{%d}$'%np.log10(mu))

    # ax.set_title("Simulated Poisson Processes (Gillespie Algorithm)")
    ax.tick_params(labelsize = 18)
    ax.set_xlabel(r"$n$", fontsize = 20)
    ax.set_ylabel(r"$\textrm{Fitness}$", fontsize = 20)
    ax.set_xscale('linear')
    ax.set_xticks(range(1, n))
    ax.legend(fontsize = 18, title = r'$\mu/\rho$', title_fontsize = 20)
    fig.savefig(output_plot+'/landscape.pdf')

    ax_c.tick_params(labelsize = 18)
    ax_c.set_xlabel(r"$n$", fontsize = 20)
    ax_c.set_ylabel(r"$\textrm{Cost}$", fontsize = 20)
    ax_c.set_xscale('log')
    # ax_c.set_xticks(range(1, n+1))
    ax_c.legend(fontsize = 18, title = r'$\mu/\rho$', title_fontsize = 20)
    fig_c.savefig(output_plot+'/Cost.pdf')

    # Print Final execution time
    end_time = time.time()
    print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()

