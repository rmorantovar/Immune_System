import sys
sys.path.append('../../my_lib/')
from funcs import*
import random

# ---------- ODE System for Full Repertoire ----------
def peptide_repertoire_ode(t, y, kin, kout, kon, kstep, koffs):
    N = len(koffs)
    p = y[0:N]
    pMHC = y[N:2*N]
    pMHCp = y[2*N:3*N]

    dp = (kin / N) - kout * p - kon * p + koffs * pMHC
    dpMHC = kon * p - koffs * pMHC - kstep * pMHC
    dpMHCp = kstep * pMHC - koffs * pMHCp

    return np.concatenate([dp, dpMHC, dpMHCp])

def run_ode_simulation(N=10000, kin=10, kout=0.1, kon=10, kstep=0.1):
    koffs = np.logspace(-4, 1, N)
    initial_state = np.zeros(3 * N)
    t_span = (0, 500)
    t_eval = np.linspace(*t_span, 100)

    sol = solve_ivp(
        peptide_repertoire_ode, t_span, initial_state,
        args=(kin, kout, kon, kstep, koffs), t_eval=t_eval
    )

    p_ss = sol.y[0:N, -1]
    pMHC_ss = sol.y[N:2*N, -1]
    pMHCp_ss = sol.y[2*N:3*N, -1]

    return koffs, p_ss, pMHC_ss, pMHCp_ss


# ---------- Binned Gillespie Model ----------
def make_reaction_functions():
    def import_peptide(state):
        s = state.copy(); s[0] += 1; return s
    def export_peptide(state):
        s = state.copy(); s[0] -= 1; return s
    def bind(state):
        s = state.copy(); s[0] -= 1; s[1] += 1; return s
    def unbind(state):
        s = state.copy(); s[0] += 1; s[1] -= 1; return s
    def step(state):
        s = state.copy(); s[1] -= 1; s[2] += 1; return s
    def unbind_pMHCp(state):
        s = state.copy(); s[2] -= 1; return s
    return [import_peptide, export_peptide, bind, unbind, step, unbind_pMHCp]

def make_propensity_functions(kin, kout, kon, koff, kstep):
    return [
        lambda state: kin,
        lambda state: kout * state[0],
        lambda state: kon * state[0],
        lambda state: koff * state[1],
        lambda state: kstep * state[1],
        lambda state: koff * state[2]
    ]

def gillespie_algorithm(initial_state, reactions, propensities, max_time):
    state = np.array(initial_state, dtype=int)
    time = 0.0
    times = [time]
    states = [state.copy()]
    
    while time < max_time:
        a = np.array([p(state) for p in propensities])
        a0 = np.sum(a)
        if a0 == 0: break
        r1, r2 = np.random.random(), np.random.random()
        tau = -np.log(r1) / a0
        reaction_index = np.searchsorted(np.cumsum(a), r2 * a0)
        state = reactions[reaction_index](state)
        time += tau
        times.append(time)
        states.append(state.copy())
    return times, np.array(states)

def run_binned_stochastic_simulation(koffs, nbins=50, kin_total=10, kout=0.1, kon=10, kstep=0.1):
    bin_edges = np.logspace(np.log10(koffs.min()), np.log10(koffs.max()), nbins + 1)
    counts, _ = np.histogram(koffs, bins=bin_edges)
    bin_koffs = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    results = []

    for count, koff_bin in zip(counts, bin_koffs):
        if count == 0: continue

        kin_bin = (kin_total * count) / len(koffs)
        reactions = make_reaction_functions()
        props = make_propensity_functions(kin_bin, kout, kon, koff_bin, kstep)
        times, states = gillespie_algorithm([0, 0, 0], reactions, props, max_time=400)

        n = len(states)
        mean_p = np.mean(states[int(n/2):, 0])
        mean_pMHC = np.mean(states[int(n/2):, 1])
        mean_pMHCp = np.mean(states[int(n/2):, 2])
        results.append((koff_bin, count, mean_p, mean_pMHC, mean_pMHCp))

    results = np.array(results)
    return results


# ---------- Plotting ----------
def plot_ode(ax, koffs, p, pMHC, pMHCp):

    ax.loglog(koffs, p, '.', label='p')
    ax.loglog(koffs, pMHC, '.', label='pMHC')
    ax.loglog(koffs, pMHCp, '.', label='pMHCp')
    

def plot_binned(ax, results):
    koff_bins, counts, p, pMHC, pMHCp = results.T
    weighted_p = counts * p
    weighted_pMHC = counts * pMHC
    weighted_pMHCp = counts * pMHCp

    ax.loglog(koff_bins, weighted_p/N, 'o-', label='p')
    ax.loglog(koff_bins, weighted_pMHC/N, 'o-', label='pMHC')
    ax.loglog(koff_bins, weighted_pMHCp/N, 'o-', label='pMHCp')
    


# ---------- Main ----------
if __name__ == "__main__":
    # Shared parameters
    N = 10000
    kin = 10
    kout = 0.1
    kon = 10
    kstep = 0.1
    koffs = np.logspace(-4, 1, N)
    fig, ax = plt.subplots(figsize=(8*1.62, 8), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

    # ---- Run ODE Simulation ----
    print("Running ODE simulation...")
    koffs_ode, p, pMHC, pMHCp = run_ode_simulation(N=N, kin=kin, kout=kout, kon=kon, kstep=kstep)
    plot_ode(ax, koffs_ode, p, pMHC, pMHCp)

    # ---- Run Binned SSA Simulation ----
    print("Running binned Gillespie simulation...")
    results_binned = run_binned_stochastic_simulation(koffs, nbins=50, kin_total=kin, kout=kout, kon=kon, kstep=kstep)
    plot_binned(ax, results_binned)

    # plt.xlabel('koff')
    # plt.ylabel('Concentration')
    # plt.legend()
    # plt.title("ODE Steady-State")
    # plt.grid(True, which='both', ls='--')
    # plt.tight_layout()

    # plt.xlabel('koff (bin avg)')
    # plt.ylabel('Weighted mean count')
    # plt.legend()
    # plt.grid(True, which='both', ls='--')
    # plt.title("Binned Stochastic Steady-State")
    # plt.tight_layout()

    my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    ax.legend(fontsize = 20, loc = 1)
    fig.savefig('/Users/robertomorantovar/Library/CloudStorage/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/presentation_MHC/sim_rep-%.2f.pdf'%kstep)

