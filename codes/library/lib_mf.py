
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import scipy.stats as stats
from dataclasses import dataclass
from cycler import cycler

N_Avg = 6.02214076e23  # Avogadro's number (molecules per mole)
k_BT = 1.380649e-23*293

my_red = np.array((228,75,41))/256.
my_purple = np.array((125,64,119))/256.
my_purple2 = np.array((116,97,164))/256.
my_green = np.array((125,165,38))/256.
my_blue = np.array((76,109,166))/256.
my_gold = np.array((215,139,45))/256.
my_brown = np.array((182,90,36))/256.
my_blue2 = np.array((80,141,188))/256.
my_yellow = np.array((246,181,56))/256.
my_yellow2 = np.array((242, 192, 65))/256.
my_green2 = np.array((158,248,72))/256.
my_cyan = 'tab:cyan'

antigen_color = my_yellow

my_green_a = np.array((159, 206, 99))/256.
my_green_b = np.array((79, 173, 91))/256.
my_green_c = np.array((94, 129, 63))/256.

plt.rcParams['text.usetex'] = True

# Option 1: via rcParams with a named colormap
plt.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.tab10.colors)

# Option 2: use tab20 if you need more distinct colors
# plt.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.tab20.colors)

# Option 3: apply a built-in style that uses tab colors
# plt.style.use('tableau-colorblind10')   # or 'default' which uses tab10


# ============================================================
# Parameters
# ============================================================

@dataclass
class Parameters:
    """All model parameters in one place."""

    # --- Antigen ---
    N_A0: float = 1.0           # initial antigen concentration
    lambda_A: float = 1.0       # antigen growth rate (>0: replication, <0: decay)
    delta_A: float = 1.0        # antigen clearance rate (set to 0 for pure exponential)

    # --- pMHC ---
    k_on: float = 1.0           # antigen encounter rate (affinity-independent)
    delta_pi: float = 1.0       # pMHC surface turnover rate (fast)
    pi_star: float = 1.0           # pMHC half-maximal threshold
    hill: float = 3.0             # Hill coefficient for T cell activation (controls sharpness of transition)

    # --- Binding ---
    eta: float = 1.0          # specificity parameter: psi(DG) = exp(-eta * DG)
    beta_star: float = 2.0      # density-of-states exponent

    # --- T cells ---
    N_T0: float = 1000.0         # initial T cell pool
    delta_T: float = 0.0        # T cell turnover rate (0 = constant pool)
    h0: float = 0.1         # T-B contact rate constant
    tau_eng: float = 0.1        # T-B engagement duration
    b0: float = 0.3            # intrinsic cell division rate
    K_T: float = 10.0           # saturation constant for T cell limitation
    Tcell_growth_factor: float = 1.0 # factor controlling T cell expansion upon activation (relative to b0)

    # --- B cells ---
    delta_B: float = 0.0        # B cell death rate

    # --- Repertoire ---
    DG_min: float = 0.0         # highest affinity (smallest DG)
    DG_max: float = 5.0        # lowest affinity in grid
    M: int = 20                # number of grid points
    Omega_0: float = 1.0        # density-of-states prefactor
    
    # --- Simulation options ---
    T_lim: bool = False         # whether to include T cell limitation in the demand function
    memory: bool = False        # whether to use memory-like initial conditions for N_B
    
# ============================================================
# Auxiliary functions
# ============================================================
_call_count = [0]  # list so we can mutate from inside the function

def build_repertoire_grid(p):
    """
    Grid mode: uniform grid in DG, with Omega weighting.
 
    Returns
    -------
    DG_arr : array of DG values (length M)
    psi_arr : array of psi(DG) values
    weights : array of Omega(DG) * dDG (used in demand integral)
    M : number of grid points
    """
    DG_arr = np.linspace(p.DG_min, p.DG_max, p.M)
    dDG = DG_arr[1] - DG_arr[0]
    psi_arr = psi(DG_arr, p.eta)
    Omega_arr = Omega(DG_arr, p.beta_star, p.Omega_0)
    weights = Omega_arr * dDG  # each bin represents this many clones
    return DG_arr, psi_arr, weights, p.M


def build_repertoire_stochastic(p, DG_max_sim=None, seed=None):
    """
    Stochastic mode: sample individual clone energies from Omega(DG).
 
    Samples from the truncated exponential distribution on [DG_min, DG_max_sim].
    The number of clones is determined by integrating Omega over the range.
 
    Parameters
    ----------
    p : Parameters
    DG_max_sim : float
        Upper bound of the simulation range. If None, uses p.DG_max.
    seed : int or None
        Random seed for reproducibility.
 
    Returns
    -------
    DG_arr : array of DG values (length L_sim), sorted
    psi_arr : array of psi(DG) values
    weights : array of ones (each entry is one clone)
    L_sim : number of clones
    """
    if DG_max_sim is None:
        DG_max_sim = p.DG_max
 
    rng = np.random.default_rng(seed)
    n_extra = 3
    # Replace DG_min with DG_floor in the sampling
    DG_floor = p.DG_min - n_extra / p.beta_star  # extend by n_extra expected clones
 
    # Total number of clones in [DG_min, DG_max_sim]:
    # L_sim = int(Omega_0 / beta_star * (exp(beta_star * DG_max_sim) - exp(beta_star * DG_min)))
    L_sim = int(p.Omega_0 / p.beta_star * (
        np.exp(p.beta_star * DG_max_sim) - np.exp(p.beta_star * DG_floor)
    ))
 
    if L_sim <= 0:
        raise ValueError(f"No clones in range [{p.DG_min}, {DG_max_sim}]. "
                         f"Check Omega_0 and beta_star.")
 
    # Sample from truncated exponential: P(DG) ~ exp(beta_star * DG)
    # on [DG_min, DG_max_sim].
    # Use inverse CDF method:
    # CDF(DG) = (exp(beta_star * DG) - exp(beta_star * DG_min)) /
    #           (exp(beta_star * DG_max_sim) - exp(beta_star * DG_min))
    # DG = (1/beta_star) * ln(u * (exp(beta_star*DG_max_sim) - exp(beta_star*DG_min))
    #                         + exp(beta_star*DG_min))
 
    u = rng.uniform(0, 1, size=L_sim)
    exp_min = np.exp(p.beta_star * DG_floor)
    exp_max = np.exp(p.beta_star * DG_max_sim)
    DG_arr = (1.0 / p.beta_star) * np.log(u * (exp_max - exp_min) + exp_min)
    DG_arr = np.sort(DG_arr)
 
    psi_arr = psi(DG_arr, p.eta)
    weights = np.ones(L_sim)  # each entry is one clone
 
    return DG_arr, psi_arr, weights, L_sim
 

def psi(DG, eta):
    """Binding/internalization probability. Boltzmann gate."""
    return (1 + np.exp(DG+0.5))**-eta


def G1(pi_vec, p):
    """B-T cell engagement probability."""
    # return pi_vec**p.hill / (pi_vec**p.hill + p.Theta**p.hill)
    return pi_vec


def Omega(DG, beta_star, Omega_0):
    """Density of states (number of clones per unit DG)."""
    return Omega_0 * np.exp(beta_star * DG)


def compute_demand(pi_vec, N_B_vec, weights, p):
    """
    Compute the demand function D(t).
    Grid mode:   D = tau_eng * h0 * sum_i Omega_i * dDG * vis_i * N_B_i
    Stochastic:  D = tau_eng * h0 * sum_i vis_i * N_B_i
    In both cases, 'weights' absorbs the difference:
      grid:       weights_i = Omega_i * dDG
      stochastic: weights_i = 1
    """
    activated = (N_B_vec > 2).astype(float)
    # return p.tau_eng * p.h0 * np.sum(weights * visibility * N_B_vec * activated)
    return np.sum(weights * G1(pi_vec, p) * N_B_vec)


def compute_h_T(pi_vec, N_T, N_B_vec, weights, p):
    """Compute help rate h_T."""
    # h_T = h0 * G1(pi_vec, p) * N_T_free/ (N_T_free + p.K_T)  # K_T is an arbitrary saturation constant to prevent unbounded growth of h_T

    D = compute_demand(pi_vec, N_B_vec, weights, p)
    if p.T_lim:
        N_T_free = N_T / (1.0 + D/p.K_T)  # T cell limitation: free T cells decrease as demand increases
    else:
        N_T_free = N_T 

    h_T = p.h0 * G1(pi_vec, p) * N_T_free #linear regime for \pi

    return h_T


def compute_h_B(pi_vec, N_B_vec, weights, p):
    """Compute help rate h_B."""
    
    D = compute_demand(pi_vec, N_B_vec, weights, p)
    
    h_B =  p.h0 *  D/(p.K_T + D) #linear regime

    return h_B


def compute_lambda_B(pi_vec, N_B_vec, N_T, weights, p):
    """
    Compute division rate for each clone.

    lambda_B = (1/h_T + 1/b0)^{-1}
    h_T = h0 * pi^h/(pi^h+Theta^h) * N_T_free/(N_T_free + K_T)
    """
    # h_T = p.h0 * G1(pi_vec, p) * N_T_free/ (N_T_free + p.K_T)  # K_T is an arbitrary saturation constant to prevent unbounded growth of h_T
    # h_T = p.h0 * G1(pi_vec, p) * N_T_free #linear regime   

    h_T = compute_h_T(pi_vec, N_T, N_B_vec, weights, p)

    # Avoid division by zero when h_T = 0
    with np.errstate(divide='ignore', invalid='ignore'):
        inv_lambda = np.where(h_T > 0, (1.0 / h_T) + (1.0 / p.b0), np.inf)
        lambda_B = np.where(inv_lambda < np.inf, 1.0 / inv_lambda, 0.0)

    return lambda_B


def compute_lambda_T(pi_vec, N_B_vec, weights, p):
    """
    Compute division rate for T cells.

    lambda_T = (1/h_B + 1/b0)^{-1}
    h_B = h0 * D, where D is the demand function.
    """
    # h_B = p.h0 * G1(pi_vec, p) * N_T_free/ (N_T_free + p.K_T)  # K_T is an arbitrary saturation constant to prevent unbounded growth of h_T
    h_B =  compute_h_B(pi_vec, N_B_vec, weights, p) #linear regime

    # Avoid division by zero when h_B = 0
    with np.errstate(divide='ignore', invalid='ignore'):
        inv_lambda = np.where(h_B > 0, (1.0 / h_B) + (1.0 / p.b0), np.inf)
        lambda_T = np.where(inv_lambda < np.inf, 1.0 / inv_lambda, 0.0)

    return lambda_T


def compute_N_B_tot(res, threshold=2.0):
    """Total number of activated B cells."""
    N_B = res['N_Bo'] + res['N_Ba']
    w = res['weights']
    activated = (N_B > threshold).astype(float)
    return np.sum(w[:, None] * N_B * activated, axis=0)
 
 
def compute_L_act(res, threshold=2.0):
    """Number of activated clones."""
    N_B = res['N_Bo'] + res['N_Ba']
    w = res['weights']
    activated = (N_B > threshold).astype(float)
    return np.sum(w[:, None] * activated, axis=0)
 

def compute_N1(res):
    N_B = res['N_Bo'] + res['N_Ba']
    return np.max(N_B, axis=0)


def find_t_D(res, D_threshold=1.0):
    """Find the time at which D(t) crosses D_threshold."""
    idx = np.where(res['D'] >= D_threshold)[0]
    if len(idx) > 0:
        return res['t'][idx[0]]
    return np.inf
 

def compute_zipf(res, time_index=-1, threshold=1.5):
    """
    Compute the Zipf (rank-size) distribution at a given time.
 
    Returns
    -------
    ranks : array, 1-indexed ranks
    sizes : array, clone sizes sorted descending
    """
    N_B = res['N_Bo'][:, time_index] + res['N_Ba'][:, time_index]
    w = res['weights']
 
    if res['mode'] == 'stochastic':
        # Each entry is one clone
        activated = N_B > threshold
        sizes = N_B[activated]
        sizes = np.sort(sizes)[::-1]
        ranks = np.arange(1, len(sizes) + 1)
    else:
        # Grid mode: each bin has w[i] clones of the same size
        # Expand into individual clones for ranking
        sizes_expanded = []
        for i in range(len(N_B)):
            if N_B[i] > threshold:
                n_clones = max(1, int(np.round(w[i])))
                sizes_expanded.extend([N_B[i]] * n_clones)
        sizes = np.array(sorted(sizes_expanded, reverse=True))
        ranks = np.arange(1, len(sizes) + 1)
 
    return ranks, sizes/sizes[0]  # normalize by largest clone size


def compute_potency(res, time_index=-1, threshold=2.0):
    """
    Compute the potency (sum of the clone size times exp(-E)) at a given time.
 
    Returns
    -------
    potency : float, total activated B cells normalized by exp(-DG) weighting
    """
    N_B = res['N_Bo'][:, time_index] + res['N_Ba'][:, time_index]
    DG = res['DG']
    w = res['weights']
 
    activated = N_B > threshold
    potency = np.sum(w[activated] * N_B[activated] * np.exp(-DG[activated]))
    
    return potency 


def compute_potency_t(res, threshold=2.0):
    """Ω-weighted potency at every time point, activated cells only."""
    N_B = res['N_Bo'] + res['N_Ba']            # (M, N_t)
    w = res['weights']                          # (M,)
    DG = res['DG']                              # (M,)
    activated = (N_B > threshold).astype(float)
    Z_B = w[:, None] * N_B * np.exp(-DG[:, None]) * activated
    return np.sum(Z_B, axis=0)                  # (N_t,)


def compute_yield(res, time_index=-1, threshold=2.0):
    """
    Compute the yield (sum of the clone size) at a given time.
 
    Returns
    -------
    yield : float, total activated B cells
    """
    N_B = res['N_Bo'][:, time_index] + res['N_Ba'][:, time_index]
    DG = res['DG']
    w = res['weights']
 
    activated = N_B > threshold
    yield_ = np.sum(w[activated] * N_B[activated])
    
    return yield_ 


def produce_memory(res, n_mem=int(1e4)):
    per_clone = res['N_Bo'][:, -1] + res['N_Ba'][:, -1]   # per-clone count
    w = res['weights']
    activated_idx = np.where(per_clone > 2.0)[0]           # threshold on per-clone
    counts = (per_clone * w)[activated_idx]                # cells per bin, for sampling
    total = counts.sum()
    if total <= 0:
        return np.array([]), np.array([])
    n_draw = min(n_mem, int(total))
    n_draw = (n_mem * total)/(n_mem + total)  # adjust n_draw to avoid sampling more than available
    # print(n_draw)
    probs = counts / total
    drawn = np.random.multinomial(n_draw, probs)
    nz = drawn > 0
    DG_memory = res['DG'][activated_idx[nz]]
    N_memory = drawn[nz]
    return DG_memory, N_memory

# ============================================================
# ODE system
# ============================================================

def pack_state_approx(N_A, pi_vec, N_B_vec, N_T, M):
    """Pack all state variables into a single vector for the integrator."""
    y = np.zeros(2 * M + 2)
    y[0] = N_A
    y[1:M+1] = pi_vec
    y[M+1:2*M+1] = N_B_vec
    y[2*M+1] = N_T
    return y


def pack_state_complete(N_A, pi_vec, N_Bo_vec, N_BT_vec, N_Ba_vec, N_To, N_Ta, M):
    """Pack all state variables into a single vector for the integrator."""
    y = np.zeros(4 * M + 3)
    y[0] = N_A
    y[1:M+1] = pi_vec
    y[M+1:2*M+1] = N_Bo_vec
    y[2*M+1:3*M+1] = N_BT_vec
    y[3*M+1:4*M+1] = N_Ba_vec
    y[4*M+1] = N_To
    y[4*M+2] = N_Ta

    return y


def pack_state_semicomplete(N_A, pi_vec, N_Bo_vec, N_Ba_vec, N_To, N_Ta, M):
    """Pack all state variables into a single vector for the integrator."""
    y = np.zeros(4 * M + 3)
    y[0] = N_A
    y[1:M+1] = pi_vec
    y[M+1:2*M+1] = N_Bo_vec
    y[2*M+1:3*M+1] = N_Ba_vec
    y[3*M+1] = N_To
    y[3*M+2] = N_Ta

    return y


def pack_state_null(N_A, N_Bo_vec, N_Ba_vec, M):
    """Pack all state variables into a single vector for the integrator."""
    y = np.zeros(4 * M + 3)
    y[0] = N_A
    y[1:M+1] = N_Bo_vec
    y[M+1:2*M+1] = N_Ba_vec

    return y


def unpack_state_approx(y, M):
    """Unpack the state vector."""
    N_A = y[0]
    pi_vec = y[1:M+1]
    N_B_vec = y[M+1:2*M+1]
    N_T = y[2*M+1]
    return N_A, pi_vec, N_B_vec, N_T


def unpack_state_complete(y, M):
    """Unpack the state vector."""
    N_A = y[0]
    pi_vec = y[1:M+1]
    N_Bo_vec = y[M+1:2*M+1]
    N_BT_vec = y[2*M+1:3*M+1]
    N_Ba_vec = y[3*M+1:4*M+1]
    N_To = y[4*M+1]
    N_Ta = y[4*M+2]
    return N_A, pi_vec, N_Bo_vec, N_BT_vec, N_Ba_vec, N_To, N_Ta


def unpack_state_semicomplete(y, M):
    """Unpack the state vector."""
    N_A = y[0]
    pi_vec = y[1:M+1]
    N_Bo_vec = y[M+1:2*M+1]
    N_Ba_vec = y[2*M+1:3*M+1]
    N_To = y[3*M+1]
    N_Ta = y[3*M+2]
    return N_A, pi_vec, N_Bo_vec, N_Ba_vec, N_To, N_Ta


def unpack_state_null(y, M):
    """Unpack the state vector."""
    N_A = y[0]
    N_Bo_vec = y[1:M+1]
    N_Ba_vec = y[M+1:2*M+1]
    return N_A, N_Bo_vec, N_Ba_vec


def rhs_approx(t, y, p, M, psi_vec, weights):
    """Right-hand side of the coupled ODE system."""
 
    N_A, pi_vec, N_B_vec, N_T = unpack_state_approx(y, M)
 
    # Ensure non-negativity
    N_A = max(N_A, 0.0)
    pi_vec = np.maximum(pi_vec, 0.0)
    N_B_vec = np.maximum(N_B_vec, 0.0)
    N_T = max(N_T, 0.0)
 
    # --- Antigen ---
    # S_A = p.lambda_A * N_A if p.lambda_A > 0 else 0.0
    # dN_A = S_A - p.delta_A * N_A
    pb = (1 + (1e-9/(1e6*24*3600*np.exp(2.0*t)/N_Avg)))**(-1)  # or whatever dependence you intend
    dN_A = (p.lambda_A * (1 - pb) - p.delta_A*pb) * N_A - 0.01 * N_A
 
    # --- Demand and free T cells ---
    D = compute_demand(pi_vec, N_B_vec, weights, p)
    # change here to test T cell limitation vs no limitation
    if p.T_lim:
        N_T_free = N_T / (1.0 + D/p.K_T)  # T cell limitation: free T cells decrease as demand increases
    else:
        N_T_free = N_T 
 
    # --- Division rate ---
    # lambda_B = compute_lambda_B(pi_vec, N_T_free, p)
    lambda_B = compute_lambda_B(pi_vec, N_B_vec, N_T, weights, p)
    # print(lambda_B)
 
    # --- pMHC dynamics ---
    dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_B * pi_vec
 
    # --- Clone size dynamics ---
    dN_B = (lambda_B - p.delta_B) * N_B_vec
 
    # --- T cell dynamics ---
    lambda_T = compute_lambda_T(pi_vec, N_B_vec, weights, p)
    # print(lambda_T)
    dN_T = lambda_T * N_T  - p.delta_T * N_T
    # dN_T = (- p.delta_T) * N_T
 
    return pack_state_approx(dN_A, dpi, dN_B, dN_T, M)


def rhs_complete(t, y, p, M, psi_vec, weights):
    _call_count[0] += 1
    if _call_count[0] % 1000 == 0:
        print(f"  t = {t:.4f}, calls = {_call_count[0]}")
        
    """Right-hand side of the coupled ODE system."""
 
    N_A, pi_vec, N_Bo_vec, N_BT_vec, N_Ba_vec, N_To, N_Ta = unpack_state_complete(y, M)
 
    # Ensure non-negativity
    N_A = max(N_A, 0.0)
    pi_vec = np.maximum(pi_vec, 0.0)
    N_Bo_vec = np.maximum(N_Bo_vec, 0.0)
    N_BT_vec = np.maximum(N_BT_vec, 0.0)
    N_Ba_vec = np.maximum(N_Ba_vec, 0.0)
    N_To = max(N_To, 0.0)
    N_Ta = max(N_Ta, 0.0)
 
    # --- Antigen ---
    # S_A = p.lambda_A * N_A if p.lambda_A > 0 else 0.0
    # dN_A = S_A - p.delta_A * N_A
    pb = (1 + (1e-9/(1e6*24*3600*np.exp(2.0*t)/N_Avg)))**(-1)  # or whatever dependence you intend
    dN_A = (p.lambda_A * (1 - pb) - p.delta_A*pb) * N_A - 0.01 * N_A
      
    # --- pMHC dynamics ---
    N_B_tot = N_Bo_vec + N_BT_vec + N_Ba_vec
    N_T_tot = N_To + N_Ta
    lambda_eff = np.where(N_B_tot > 0, p.b0 * N_Ba_vec / N_B_tot, 0.0)
   
    dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_eff * pi_vec
 
    # --- Free B-cell clones ---
    dN_Bo = - p.k_on * pi_vec**p.hill * N_Bo_vec * N_To/(N_To + p.K_T) + 2 * p.b0 * N_Ba_vec - p.delta_B * N_Bo_vec

    # --- T-B conjugates ---
    dN_BT = p.k_on * pi_vec**p.hill * N_Bo_vec * N_To/(N_To + p.K_T) - p.h0 * N_BT_vec

    # --- Activated B-cell clones ---
    dN_Ba = p.h0 * N_BT_vec - p.b0 * N_Ba_vec  - p.delta_B * N_Ba_vec
 
    # --- Free T cells ---
    dN_To = - p.k_on * np.sum(weights * pi_vec**p.hill * N_Bo_vec) * N_To/(N_To + p.K_T) + 2 * p.b0 * N_Ta - p.delta_T * N_To

    # --- Activated T cells ---
    dN_Ta = p.h0 * np.sum(weights * N_BT_vec) - p.b0 * N_Ta - p.delta_T * N_Ta 
 
    return pack_state_complete(dN_A, dpi, dN_Bo, dN_BT, dN_Ba, dN_To, dN_Ta, M)


def rhs_semicomplete(t, y, p, M, psi_vec, weights):
    # _call_count[0] += 1
    # if _call_count[0] % 1000 == 0:
    #     print(f"  t = {t:.4f}, calls = {_call_count[0]}")
        
    """Right-hand side of the coupled ODE system."""
 
    N_A, pi_vec, N_Bo_vec, N_Ba_vec, N_To, N_Ta = unpack_state_semicomplete(y, M)
 
    # Ensure non-negativity
    N_A = max(N_A, 0.0)
    pi_vec = np.maximum(pi_vec, 0.0)
    N_Bo_vec = np.maximum(N_Bo_vec, 0.0)
    N_Ba_vec = np.maximum(N_Ba_vec, 0.0)
    N_To = max(N_To, 0.0)
    N_Ta = max(N_Ta, 0.0)
 
    # --- Antigen ---
    # S_A = p.lambda_A * N_A if p.lambda_A > 0 else 0.0
    # dN_A = S_A - p.delta_A * N_A
    pb = (1 + (1e-9/(1e6*24*3600*np.exp(2.0*t)/N_Avg)))**(-1)  # or whatever dependence you intend
    dN_A = (p.lambda_A * (1 - pb) - p.delta_A*pb) * N_A - 0.01 * N_A
      
    # --- pMHC dynamics ---
    N_B_tot = N_Bo_vec + N_Ba_vec
    N_T_tot = N_To + N_Ta
    lambda_eff = np.where(N_B_tot > 0, p.b0 * N_Ba_vec / N_B_tot, 0.0)
   
    dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_eff * pi_vec
 
    # --- Free B-cell clones ---
    # dN_Bo = - p.h0 * (pi_vec**p.hill / (pi_vec**p.hill + p.pi_star**p.hill)) * N_Bo_vec * N_To/(N_To + p.K_T) + 2 * p.b0 * N_Ba_vec - p.delta_B * N_Bo_vec
    dN_Bo = - p.h0 * pi_vec**p.hill * N_Bo_vec * N_To/(N_To + p.K_T) + 2 * p.b0 * N_Ba_vec - p.delta_B * N_Bo_vec

    # --- Activated B-cell clones ---
    # dN_Ba = p.h0 * (pi_vec**p.hill / (pi_vec**p.hill + p.pi_star**p.hill)) * N_Bo_vec * N_To/(N_To + p.K_T) - p.b0 * N_Ba_vec  - p.delta_B * N_Ba_vec
    dN_Ba = p.h0 * pi_vec**p.hill * N_Bo_vec * N_To/(N_To + p.K_T) - p.b0 * N_Ba_vec  - p.delta_B * N_Ba_vec
 
    # --- Free T cells ---
    # dN_To = - p.h0 * np.sum(weights * (pi_vec**p.hill / (pi_vec**p.hill + p.pi_star**p.hill)) * N_Bo_vec) * N_To/(N_To + p.K_T) + p.Tcell_growth_factor * p.b0 * N_Ta - p.delta_T * N_To
    dN_To = - p.h0 * np.sum(weights * pi_vec**p.hill * N_Bo_vec) * N_To/(N_To + p.K_T) + p.Tcell_growth_factor * p.b0 * N_Ta - p.delta_T * N_To

    # --- Activated T cells ---
    # dN_Ta = p.h0 * np.sum(weights * (pi_vec**p.hill / (pi_vec**p.hill + p.pi_star**p.hill)) * N_Bo_vec) * N_To/(N_To + p.K_T) - p.b0 * N_Ta - p.delta_T * N_Ta 
    dN_Ta = p.h0 * np.sum(weights * pi_vec**p.hill * N_Bo_vec) * N_To/(N_To + p.K_T) - p.b0 * N_Ta - p.delta_T * N_Ta 
 
    return pack_state_semicomplete(dN_A, dpi, dN_Bo, dN_Ba, dN_To, dN_Ta, M)


def rhs_null(t, y, p, M, psi_vec, weights):
    # _call_count[0] += 1
    # if _call_count[0] % 1000 == 0:
    #     print(f"  t = {t:.4f}, calls = {_call_count[0]}")
        
    """Right-hand side of the coupled ODE system."""
 
    N_A, N_Bo_vec, N_Ba_vec = unpack_state_null(y, M)
 
    # Ensure non-negativity
    N_A = max(N_A, 0.0)
    N_Bo_vec = np.maximum(N_Bo_vec, 0.0)
    N_Ba_vec = np.maximum(N_Ba_vec, 0.0)
 
    # --- Antigen ---
    # S_A = p.lambda_A * N_A if p.lambda_A > 0 else 0.0
    # dN_A = S_A - p.delta_A * N_A
    pb = (1 + (1e-9/(1e6*24*3600*np.exp(2.0*t)/N_Avg)))**(-1)  # or whatever dependence you intend
    dN_A = (p.lambda_A * (1 - pb) - p.delta_A*pb) * N_A - 0.01 * N_A
      
    # --- antigen processing ---
    u_p = p.k_on * psi_vec * N_A
 
    # --- Free B-cell clones ---
    # dN_Bo = - p.h0 * (pi_vec**p.hill / (pi_vec**p.hill + p.pi_star**p.hill)) * N_Bo_vec * N_To/(N_To + p.K_T) + 2 * p.b0 * N_Ba_vec - p.delta_B * N_Bo_vec
    dN_Bo = - p.h0 * u_p * N_Bo_vec + 2 * p.b0 * N_Ba_vec - p.delta_B * N_Bo_vec

    # --- Activated B-cell clones ---
    # dN_Ba = p.h0 * (pi_vec**p.hill / (pi_vec**p.hill + p.pi_star**p.hill)) * N_Bo_vec * N_To/(N_To + p.K_T) - p.b0 * N_Ba_vec  - p.delta_B * N_Ba_vec
    dN_Ba = p.h0 * u_p * N_Bo_vec - p.b0 * N_Ba_vec  - p.delta_B * N_Ba_vec
 
    return pack_state_null(dN_A, dN_Bo, dN_Ba, M)

# ============================================================
# Simulation
# ============================================================

def run_simulation_approx(p=None, t_span=None, t_eval=None, mode='grid', DG_max_sim=None, seed=None):      
    """
    Run the mean-field simulation.
 
    Parameters
    ----------
    p : Parameters
    t_span : tuple (t_start, t_end)
    t_eval : array of output time points
    mode : 'grid' or 'stochastic'
    DG_max_sim : float (stochastic mode only)
        Upper bound for clone sampling. Default: p.DG_max.
    seed : int (stochastic mode only)
        Random seed.

    Returns
    -------
    result : dict with keys:
        't'       : time array
        'N_A'     : antigen concentration
        'pi'      : pMHC array (M x len(t))
        'N_B'     : clone size array (M x len(t))
        'N_T'     : T cell count
        'N_T_free': free T cells
        'h_T'     : T-cell help rates
        'h_B'     : B-cell help rate
        'D'       : demand function
        'lambda_B': division rate array (M x len(t))
        'lambda_T': division rate array (len(t))
        'DG_grid' : Delta G grid
        'Omega'   : density of states
        'params'  : parameters used
    """
    if p is None:
        p = Parameters()
    if t_span is None:
        t_span = (0.0, 10.0)
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500)
 
    # --- Build repertoire ---
    if mode == 'grid':
        DG_arr, psi_arr, weights, M = build_repertoire_grid(p)
    elif mode == 'stochastic':
        DG_arr, psi_arr, weights, M = build_repertoire_stochastic(
            p, DG_max_sim=DG_max_sim, seed=seed
        )
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'grid' or 'stochastic'.")
 
    # --- Initial conditions ---
    N_A_init = p.N_A0
    pi_init = np.zeros(M)  # no pMHC at t=0
    if p.memory:
        N_B_init = 1e2*np.exp(-p.eta * (p.b0 / p.lambda_A + 1) * DG_arr)  # memory
    else:
        N_B_init = np.ones(M)  # naive
    N_T_init = p.N_T0
 
    y0 = pack_state_approx(N_A_init, pi_init, N_B_init, N_T_init, M)
 
    # --- Integrate ---
    sol = solve_ivp(
        fun=lambda t, y: rhs_approx(t, y, p, M, psi_arr, weights),
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8,
        atol=1e-10,
        max_step=0.01,
    )
 
    if not sol.success:
        print(f"Warning: integration failed: {sol.message}")
 
    # --- Unpack ---
    t = sol.t
    N_A, pi_arr, N_B_arr, N_T_arr = unpack_state_approx(sol.y, M)
 
    # --- Derived quantities ---
    N_steps = len(t)
    D_arr = np.zeros(N_steps)
    N_T_free_arr = np.zeros(N_steps)
    lambda_B_arr = np.zeros((M, N_steps))
    lambda_T_arr = np.zeros(N_steps)
    h_T_arr = np.zeros((M, N_steps))
    h_B_arr = np.zeros(N_steps)
    

    for j in range(N_steps):
        D_arr[j] = compute_demand(pi_arr[:, j], N_B_arr[:, j], weights, p)
        # change here to test T cell limitation vs no limitation
        if p.T_lim:
            N_T_free_arr[j] = N_T_arr[j] / (1.0 + D_arr[j]/p.K_T)  # T cell limitation: free T cells decrease as demand increases
        else:
            N_T_free_arr[j] = N_T_arr[j]
        # lambda_B_arr[:, j] = compute_lambda_B(pi_arr[:, j], N_T_free_arr[j], p)
        h_T_arr[:, j] = compute_h_T(pi_arr[:, j], N_T_arr[j], N_B_arr[:, j], weights, p)
        h_B_arr[j] = compute_h_B(pi_arr[:, j], N_B_arr[:, j], weights, p)
        lambda_B_arr[:, j] = compute_lambda_B(pi_arr[:, j], N_B_arr[:, j], N_T_arr[j], weights, p)
        lambda_T_arr[j] = compute_lambda_T(pi_arr[:, j], N_B_arr[:, j], weights, p)

    return {
        't': t,
        'N_A': N_A,
        'pi': pi_arr,
        'N_B': N_B_arr,
        'N_T': N_T_arr,
        'N_T_free': N_T_free_arr,
        'h_T': h_T_arr,
        'h_B': h_B_arr,
        'D': D_arr,
        'lambda_B': lambda_B_arr,
        'lambda_T': lambda_T_arr,
        'DG': DG_arr,
        'weights': weights,
        'params': p,
        'mode': mode,
        'M': M,
    }


def run_simulation_complete(p=None, t_span=None, t_eval=None, mode='grid', DG_max_sim=None, seed=None):      
    """
    Run the mean-field simulation.
 
    Parameters
    ----------
    p : Parameters
    t_span : tuple (t_start, t_end)
    t_eval : array of output time points
    mode : 'grid' or 'stochastic'
    DG_max_sim : float (stochastic mode only)
        Upper bound for clone sampling. Default: p.DG_max.
    seed : int (stochastic mode only)
        Random seed.

    Returns
    -------
    result : dict with keys:
        't'       : time array
        'N_A'     : antigen concentration
        'pi'      : pMHC array (M x len(t))
        'N_Bo'    : free B-cell clones (M x len(t))
        'N_BT'    : T-B conjugates (M x len(t))
        'N_Ba'    : activated B-cell clones (M x len(t))
        'N_To'    : free T cells (len(t))
        'N_Ta'    : activated T cells (len(t))
        'h_T'     : T-cell help rates
        'h_B'     : B-cell help rate
        'D'       : demand function
        'lambda_B': division rate array (M x len(t))
        'lambda_T': division rate array (len(t))
        'DG_grid' : Delta G grid
        'Omega'   : density of states
        'params'  : parameters used
    """
    if p is None:
        p = Parameters()
    if t_span is None:
        t_span = (0.0, 10.0)
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500  )
 
    # --- Build repertoire ---
    if mode == 'grid':
        DG_arr, psi_arr, weights, M = build_repertoire_grid(p)
    elif mode == 'stochastic':
        DG_arr, psi_arr, weights, M = build_repertoire_stochastic(
            p, DG_max_sim=DG_max_sim, seed=seed
        )
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'grid' or 'stochastic'.")
 
    # --- Initial conditions ---
    N_A_init = p.N_A0
    pi_init = np.zeros(M)  # no pMHC at t=0
    if p.memory:
        N_Bo_init = 1e2*np.exp(-p.eta * (p.b0 / p.lambda_A + 1) * DG_arr)  # memory
    else:
        N_Bo_init = np.ones(M)  # naive
    N_BT_init = np.zeros(M)
    N_Ba_init = np.zeros(M)
    N_To_init = p.N_T0
    N_Ta_init = 0.0
 
    y0 = pack_state_complete(N_A_init, pi_init, N_Bo_init, N_BT_init, N_Ba_init, N_To_init, N_Ta_init, M)
    
    
    # --- Integrate ---
    #method='BDF', 'LSODA' o 'RK45', rtol=1e-6, atol=1e-8, max_step=0.1
    _call_count[0] = 0
    sol = solve_ivp(
        fun=lambda t, y: rhs_complete(t, y, p, M, psi_arr, weights),
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method='LSODA',
        rtol=1e-4,
        atol=1e-5
    )
    print(f"Total RHS calls: {_call_count[0]}")
    if not sol.success:
        print(f"Warning: integration failed: {sol.message}")
 
    # --- Unpack ---
    t = sol.t
    N_steps = len(t)
 
    N_A, pi_arr, N_Bo_arr, N_BT_arr, N_Ba_arr, N_To_arr, N_Ta_arr = unpack_state_complete(sol.y, M)
    
    # --- Derived quantities ---
    D_arr = np.zeros(N_steps)
    N_T_free_arr = np.zeros(N_steps)
    lambda_B_arr = np.zeros((M, N_steps))
    lambda_T_arr = np.zeros(N_steps)
    h_T_arr = np.zeros((M, N_steps))
    h_B_arr = np.zeros(N_steps)
    

    return {
        't': t,
        'N_A': N_A,
        'pi': pi_arr,
        'N_Bo': N_Bo_arr,
        'N_BT': N_BT_arr,
        'N_Ba': N_Ba_arr,
        'N_To': N_To_arr,
        'N_Ta': N_Ta_arr,
        'DG': DG_arr,
        'weights': weights,
        'params': p,
        'mode': mode,
        'M': M,
    }


def run_simulation_semicomplete(p=None, t_span=None, t_eval=None, mode='grid', DG_max_sim=None, seed=None):      
    """
    Run the mean-field simulation.
 
    Parameters
    ----------
    p : Parameters
    t_span : tuple (t_start, t_end)
    t_eval : array of output time points
    mode : 'grid' or 'stochastic'
    DG_max_sim : float (stochastic mode only)
        Upper bound for clone sampling. Default: p.DG_max.
    seed : int (stochastic mode only)
        Random seed.

    Returns
    -------
    result : dict with keys:
        't'       : time array
        'N_A'     : antigen concentration
        'pi'      : pMHC array (M x len(t))
        'N_Bo'    : free B-cell clones (M x len(t))
        'N_Ba'    : activated B-cell clones (M x len(t))
        'N_To'    : free T cells (len(t))
        'N_Ta'    : activated T cells (len(t))
        'h_T'     : T-cell help rates
        'h_B'     : B-cell help rate
        'D'       : demand function
        'lambda_B': division rate array (M x len(t))
        'lambda_T': division rate array (len(t))
        'DG_grid' : Delta G grid
        'Omega'   : density of states
        'params'  : parameters used
    """
    if p is None:
        p = Parameters()
    if t_span is None:
        t_span = (0.0, 10.0)
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500)
 
    # --- Build repertoire ---
    if mode == 'grid':
        DG_arr, psi_arr, weights, M = build_repertoire_grid(p)
    elif mode == 'stochastic':
        DG_arr, psi_arr, weights, M = build_repertoire_stochastic(
            p, DG_max_sim=DG_max_sim, seed=seed
        )
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'grid' or 'stochastic'.")
 
    # --- Initial conditions ---
    N_A_init = p.N_A0
    pi_init = np.zeros(M)  # no pMHC at t=0
    if p.memory:
        N_Bo_init = 1e3*np.exp(-p.eta * (p.b0 / p.lambda_A + 1) * DG_arr)  # memory
    else:
        N_Bo_init = np.ones(M)  # naive
    N_Ba_init = np.zeros(M)
    N_To_init = p.N_T0
    N_Ta_init = 0.0
 
    y0 = pack_state_semicomplete(N_A_init, pi_init, N_Bo_init, N_Ba_init, N_To_init, N_Ta_init, M)
    
    # --- Integrate ---
    #method='BDF', 'LSODA' o 'RK45', rtol=1e-6, atol=1e-8, max_step=0.1
    # _call_count[0] = 0
    sol = solve_ivp(
        fun=lambda t, y: rhs_semicomplete(t, y, p, M, psi_arr, weights),
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method='LSODA',
        rtol=1e-4,
        atol=1e-5
    )

    if not sol.success:
        print(f"Warning: integration failed: {sol.message}")
 
    # --- Unpack ---
    t = sol.t
    N_steps = len(t)
 
    N_A, pi_arr, N_Bo_arr, N_Ba_arr, N_To_arr, N_Ta_arr = unpack_state_semicomplete(sol.y, M)
    
    # --- Derived quantities ---
    # h_T_arr = np.zeros((M, N_steps))

    return {
        't': t,
        'N_A': N_A,
        'pi': pi_arr,
        'N_Bo': N_Bo_arr,
        'N_Ba': N_Ba_arr,
        'N_To': N_To_arr,
        'N_Ta': N_Ta_arr,
        'DG': DG_arr,
        'weights': weights,
        'params': p,
        'mode': mode,
        'M': M,
    }


def run_simulation_null(p=None, t_span=None, t_eval=None, mode='grid', DG_max_sim=None, seed=None):      
    """
    Run the mean-field simulation.
 
    Parameters
    ----------
    p : Parameters
    t_span : tuple (t_start, t_end)
    t_eval : array of output time points
    mode : 'grid' or 'stochastic'
    DG_max_sim : float (stochastic mode only)
        Upper bound for clone sampling. Default: p.DG_max.
    seed : int (stochastic mode only)
        Random seed.

    Returns
    -------
    result : dict with keys:
        't'       : time array
        'N_A'     : antigen concentration
        'N_Bo'    : free B-cell clones (M x len(t))
        'N_Ba'    : activated B-cell clones (M x len(t))
        'lambda_B': division rate array (M x len(t))
        'lambda_T': division rate array (len(t))
        'DG_grid' : Delta G grid
        'Omega'   : density of states
        'params'  : parameters used
    """
    if p is None:
        p = Parameters()
    if t_span is None:
        t_span = (0.0, 10.0)
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500)
 
    # --- Build repertoire ---
    if mode == 'grid':
        DG_arr, psi_arr, weights, M = build_repertoire_grid(p)
    elif mode == 'stochastic':
        DG_arr, psi_arr, weights, M = build_repertoire_stochastic(
            p, DG_max_sim=DG_max_sim, seed=seed
        )
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'grid' or 'stochastic'.")
 
    # --- Initial conditions ---
    N_A_init = p.N_A0
    if p.memory:
        N_Bo_init = 1e3*np.exp(-p.eta * (p.b0 / p.lambda_A + 1) * DG_arr)  # memory
    else:
        N_Bo_init = np.ones(M)  # naive
    N_Ba_init = np.zeros(M)
 
    y0 = pack_state_null(N_A_init, N_Bo_init, N_Ba_init, M)
    
    # --- Integrate ---
    #method='BDF', 'LSODA' o 'RK45', rtol=1e-6, atol=1e-8, max_step=0.1
    # _call_count[0] = 0
    sol = solve_ivp(
        fun=lambda t, y: rhs_null(t, y, p, M, psi_arr, weights),
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method='LSODA',
        rtol=1e-4,
        atol=1e-5
    )

    if not sol.success:
        print(f"Warning: integration failed: {sol.message}")
 
    # --- Unpack ---
    t = sol.t
    N_steps = len(t)
 
    N_A, N_Bo_arr, N_Ba_arr = unpack_state_null(sol.y, M)
    
    return {
        't': t,
        'N_A': N_A,
        'N_Bo': N_Bo_arr,
        'N_Ba': N_Ba_arr,
        'DG': DG_arr,
        'weights': weights,
        'params': p,
        'mode': mode,
        'M': M,
    }

# ============================================================
# Visualization
# ============================================================

def plot_T_cell_analysis(res):
    """Basic diagnostic plots."""

    t = res['t']
    DG = res['DG']
    p = res['params']

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # --- (a) Antigen ---
    ax = axes[0, 0]
    ax.semilogy(t, res['N_A'], color = 'k')
    ax.semilogy(t, np.exp(p.lambda_A*t))
    ax.set_xlabel('Time')
    ax.set_ylabel('$N_A(t)$')
    ax.set_ylim(top = 1e13)
    ax.set_title('Antigen')

    # --- (b) Demand ---
    ax = axes[0, 1]
    ax.semilogy(t, res['D'], color = 'k', label='simulation')
    ax.hlines(1.0, t[0], t[-1], ls = '--', label='D=1 threshold')
    ax.set_xlabel('Time')
    ax.set_ylabel('$D(t)$')
    ax.set_title('Demand function')
    ax.set_ylim(top=1e4)
    ax.set_ylim(bottom=1e-4)

    # --- (c) Free T cells ---
    ax = axes[0, 2]
    ax.semilogy(t, res['N_T_free']/res['N_T'], label='$N_T^{o}/N_T$', color = 'k')
    ax.set_xlabel('Time')
    ax.set_ylabel('free T cells/T cells')
    ax.set_title('T cells')
    ax.set_ylim(bottom=1e-4)
    ax.legend()

    # --- (d) Division rate at selected times ---
    ax = axes[1, 0]
    n_snapshots = 8
    time_indices = np.linspace(0, len(DG)-1, n_snapshots, dtype=int)
    colors = plt.cm.Greens_r(np.linspace(0, 1, n_snapshots))
    for idx, color in zip(time_indices, colors):
        ax.plot(t, res['lambda_B'][idx, :]/p.b0, label=f'$\\Delta G$={DG[idx]:.1f}', alpha=0.8, color=color)

    ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Time')
    ax.set_ylabel('$\\lambda_B/b_o$')
    ax.set_title('Division rate')
    ax.legend(fontsize=8)

    # # --- (e) Total B-cell population and activated B-cell clones ---
    ax = axes[1, 1]
    N_B_tot = compute_N_B_tot(res)
    Lact = compute_L_act(res)
    ax.semilogy(t, N_B_tot, label='$N_B^{\\rm tot}$', color = 'k')
    ax.semilogy(t, Lact, label='$L_{\\rm act}$', color = 'r')
    ax.set_xlabel('Time')
    ax.set_ylabel('Activated B cells / clones')
    ax.set_title('B cell activation')
    ax.legend() 
    ax.set_xlabel('Time')
    ax.set_ylabel('Total activated B cells')
    ax.set_title('Total activated B cells')
    ax.set_ylim(top=1e5)

    # --- (f) Clone sizes at final time ---
    ax = axes[1, 2]
    N_B_final = res['N_B'][:, -1]
    if not p.memory:
        ax.semilogy(DG, N_B_final, marker='o', ls = '', label='simulation', ms = 4, markerfacecolor = my_red, markeredgecolor='k', markeredgewidth=0.5)
    else:
        ax.semilogy(DG, N_B_final, marker='o', ls = '', label='simulation', ms = 4, markerfacecolor = my_blue, markeredgecolor='k', markeredgewidth=0.5)
    ax.semilogy(DG, np.max(N_B_final)*np.exp(-p.eta*(p.b0/p.lambda_A + 1) * (DG - np.min(DG))), color=my_red, linestyle='--', label='naive')
    ax.semilogy(DG, np.max(N_B_final)*np.exp(-2*p.eta*(p.b0/p.lambda_A + 1) * (DG - np.min(DG))), color=my_blue, linestyle='--', label='memory')
    ax.set_xlabel('$\\Delta G$')
    ax.set_ylabel('$N_B(t_{\\rm final}, \\Delta G)$')
    ax.set_title('Clone size distribution')
    ax.set_ylim(bottom=0.9)
    # ax.set_xlim(left=-0.1, right=DG[N_B_final > 2.0][-1] + 0.5)
    ax.legend(loc = 3)
    
    axin1 = ax.inset_axes([0.6, 0.6, 0.38, 0.38])
    if res['mode'] == 'grid':
        axin1.semilogy(DG, res['weights'])
    else:
        axin1.hist(res['DG'], bins=10, density=False, color=my_green, alpha=0.7)
        axin1.set_yscale('log')
    axin1.set_xlabel('$\\Delta G$')
    axin1.set_ylabel('$\\Omega_0(\\Delta G)$')

    # # --- (e) Clone size heatmap ---
    # ax = axes[1, 1]
    # log_NB = np.log10(np.maximum(res['N_B'], 1.0))
    # im = ax.pcolormesh(t, DG, log_NB, shading='auto', cmap='viridis')
    # fig.colorbar(im, ax=ax, label='$\\log_{10} N_B$')
    # ax.set_xlabel('Time')
    # ax.set_ylabel('$\\Delta G$')
    # ax.set_title('Clone size (log)')
    # ax.set_xlim(right = 8)

    # # Overlay theoretical front
    # Gamma = p.h0 * p.N_T0 * p.k_on * p.N_A0 / (p.delta_pi * p.Theta)
    # if Gamma > 0 and p.lambda_A > 0:
    #     DG_front = (p.lambda_A / p.eta) * t - (1.0 / p.eta) * np.log(p.lambda_A / Gamma) # correct to the right off-set of the moving front!!!!
    #     DG_front_clipped = np.clip(DG_front, p.DG_min, p.DG_max) #change here p.DG_max for the actual front position, not the grid limit!!!!
    #     valid = DG_front >= p.DG_min
    #     ax.plot(t[valid], DG_front_clipped[valid], 'r--', linewidth=2,
    #             label='$\\Delta G_{\\rm mf}(t)$')
    #     ax.legend()

    
    
    plt.tight_layout()
    return fig


def plot_results(res):
    """Basic diagnostic plots."""

    t = res['t']
    DG = res['DG']
    p = res['params']

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # --- (a) Antigen ---
    ax = axes[0, 0]
    ax.semilogy(t, res['N_A'], color = 'k')
    ax.semilogy(t, np.exp(p.lambda_A*t))
    ax.set_xlabel('Time')
    ax.set_ylabel('$N_A(t)$')
    ax.set_ylim(top = 1e13)
    ax.set_title('Antigen')

    # --- (b) Free T cells ---
    ax = axes[0, 1]
    ax.semilogy(t, res['N_T_free']/res['N_T'], label='$N_T^{o}/N_T$', color = 'k')
    ax.set_xlabel('Time')
    ax.set_ylabel('free T cells/T cells')
    ax.set_title('T cells')
    ax.legend()

    # --- (c) Demand ---
    ax = axes[0, 2]
    ax.semilogy(t, res['D'], color = 'k', label='simulation')
    ax.hlines(1.0, t[0], t[-1], ls = '--', label='D=1 threshold')
    ax.set_xlabel('Time')
    ax.set_ylabel('$D(t)$')
    ax.set_title('Demand function')
    ax.set_ylim(bottom=1e-4)

    # --- (d) Clone sizes at final time ---
    ax = axes[1, 0]
    N_B_final = res['N_B'][:, -1]
    if not p.memory:
        ax.semilogy(DG, N_B_final, marker='o', ls = '', label='simulation', ms = 4, markerfacecolor = my_red, markeredgecolor='k', markeredgewidth=0.5)
    else:
        ax.semilogy(DG, N_B_final, marker='o', ls = '', label='simulation', ms = 4, markerfacecolor = my_blue, markeredgecolor='k', markeredgewidth=0.5)
    ax.semilogy(DG, np.max(N_B_final)*np.exp(-p.eta*(p.b0/p.lambda_A + 1) * (DG - np.min(DG))), color=my_red, linestyle='--', label='naive')
    ax.semilogy(DG, np.max(N_B_final)*np.exp(-2*p.eta*(p.b0/p.lambda_A + 1) * (DG - np.min(DG))), color=my_blue, linestyle='--', label='memory')
    ax.set_xlabel('$\\Delta G$')
    ax.set_ylabel('$N_B(t_{\\rm final}, \\Delta G)$')
    ax.set_title('Clone size distribution')
    ax.set_ylim(bottom=0.9)
    # ax.set_xlim(left=-0.1, right=DG[N_B_final > 2.0][-1] + 0.5)
    ax.legend(loc = 3)
    
    axin1 = ax.inset_axes([0.6, 0.6, 0.38, 0.38])
    if res['mode'] == 'grid':
        axin1.semilogy(DG, res['weights'])
    else:
        axin1.hist(res['DG'], bins=10, density=False, color=my_green, alpha=0.7)
        axin1.set_yscale('log')
    axin1.set_xlabel('$\\Delta G$')
    axin1.set_ylabel('$\\Omega_0(\\Delta G)$')

    # # --- (e) Clone size heatmap ---
    # ax = axes[1, 1]
    # log_NB = np.log10(np.maximum(res['N_B'], 1.0))
    # im = ax.pcolormesh(t, DG, log_NB, shading='auto', cmap='viridis')
    # fig.colorbar(im, ax=ax, label='$\\log_{10} N_B$')
    # ax.set_xlabel('Time')
    # ax.set_ylabel('$\\Delta G$')
    # ax.set_title('Clone size (log)')
    # ax.set_xlim(right = 8)

    # # Overlay theoretical front
    # Gamma = p.h0 * p.N_T0 * p.k_on * p.N_A0 / (p.delta_pi * p.Theta)
    # if Gamma > 0 and p.lambda_A > 0:
    #     DG_front = (p.lambda_A / p.eta) * t - (1.0 / p.eta) * np.log(p.lambda_A / Gamma) # correct to the right off-set of the moving front!!!!
    #     DG_front_clipped = np.clip(DG_front, p.DG_min, p.DG_max) #change here p.DG_max for the actual front position, not the grid limit!!!!
    #     valid = DG_front >= p.DG_min
    #     ax.plot(t[valid], DG_front_clipped[valid], 'r--', linewidth=2,
    #             label='$\\Delta G_{\\rm mf}(t)$')
    #     ax.legend()

    # # --- (e) Total B-cell population and activated B-cell clones ---
    ax = axes[1, 1]
    N_B_tot = compute_N_B_tot(res)
    Lact = compute_L_act(res)
    ax.semilogy(t, N_B_tot, label='$N_B^{\\rm tot}$', color = 'k')
    ax.semilogy(t, Lact, label='$L_{\\rm act}$', color = 'r')
    ax.set_xlabel('Time')
    ax.set_ylabel('Activated B cells / clones')
    ax.set_title('B cell activation')
    ax.legend() 
    ax.set_xlabel('Time')
    ax.set_ylabel('Total activated B cells')
    ax.set_title('Total activated B cells')

    
    # --- (f) Division rate at selected times ---
    ax = axes[1, 2]
    n_snapshots = 8
    time_indices = np.linspace(0, len(DG)-1, n_snapshots, dtype=int)
    colors = plt.cm.Greens_r(np.linspace(0, 1, n_snapshots))
    for idx, color in zip(time_indices, colors):
        ax.plot(t, res['lambda_B'][idx, :]/p.b0, label=f'$\\Delta G$={DG[idx]:.1f}', alpha=0.8, color=color)

    ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Time')
    ax.set_ylabel('$\\lambda_B/b_o$')
    ax.set_title('Division rate')
    ax.legend(fontsize=8)

    plt.tight_layout()
    return fig


def plot_diagnostics(res):
    """Basic diagnostic plots."""
    t = res['t']
    p = res['params']
 
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
 
    # (a) Antigen
    ax = axes[0, 0]
    ax.semilogy(t, res['N_A'])
    ax.set_xlabel('Time'); ax.set_ylabel('$N_A$'); ax.set_title('Antigen')
 
    # (b) T cells
    ax = axes[0, 1]
    ax.plot(t, res['N_T_free'], label='$N_T^{free}$')
    ax.plot(t, res['N_T'], '--', alpha=0.5, label='$N_T$')
    ax.set_xlabel('Time'); ax.set_ylabel('T cells'); ax.set_title('T cell pool')
    ax.legend()
 
    # (c) Demand
    ax = axes[0, 2]
    ax.semilogy(t, res['D'] + 1e-10)
    ax.axhline(1, color='k', linestyle=':', alpha=0.5)
    ax.set_xlabel('Time'); ax.set_ylabel('$D(t)$'); ax.set_title('Demand')
 
    # (d) N_B_tot and L_act
    ax = axes[1, 0]
    N_B_tot = compute_N_B_tot(res)
    L_act = compute_L_act(res)
    ax.semilogy(t, N_B_tot + 0.1, label='$\\bar{N}$')
    ax.semilogy(t, L_act + 0.1, label='$L_{act}$')
    ax.set_xlabel('Time'); ax.set_ylabel('Count'); ax.set_title('B cell response')
    ax.legend()
 
    # (e) Clone size vs DG at final time
    ax = axes[1, 1]
    DG = res['DG']
    N_B_final = res['N_B'][:, -1]
    if res['mode'] == 'stochastic':
        ax.scatter(DG, N_B_final, s=1, alpha=0.3)
    else:
        ax.semilogy(DG, N_B_final)
    ax.set_yscale('log')
    ax.set_xlabel('$\\Delta G$'); ax.set_ylabel('$N_B$')
    ax.set_title('Clone sizes (final)')
 
    # (f) Zipf plot at final time
    ax = axes[1, 2]
    ranks, sizes = compute_zipf(res, time_index=-1)
    if len(ranks) > 0:
        ax.loglog(ranks, sizes, '.', markersize=2)
        ax.set_xlabel('Rank $k$'); ax.set_ylabel('$N_k$')
        ax.set_title('Zipf plot (final)')
 
    plt.tight_layout()
    return fig


def plot_zip_comparison(res_grid, res_stochastic):
    """Compare Zipf plots for grid vs stochastic modes."""
    fig, ax = plt.subplots(figsize=(6, 4))
 
    # Grid mode
    ranks_g, sizes_g = compute_zipf(res_grid, time_index=-1)
    ax.loglog(ranks_g, sizes_g, 'o', label='Grid', alpha=0.5)
 
    # Stochastic mode
    ranks_s, sizes_s = compute_zipf(res_stochastic, time_index=-1)
    ax.loglog(ranks_s, sizes_s, '.', label='Stochastic', alpha=0.5)
 
    ax.set_xlabel('Rank $k$')
    ax.set_ylabel('Normalized clone size $N_k/N_1$')
    ax.set_title('Zipf plot comparison')
    ax.legend()
    plt.tight_layout()
    return fig