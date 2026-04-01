"""
Mean-field simulation of T-cell-mediated B-cell clonal expansion.

Simulates the coupled system:
    dN_A/dt = S_A(t) - delta_A * N_A
    dpi/dt  = k_on * psi(DG) * N_A - delta_pi * pi - lambda_B * pi
    dN_B/dt = (lambda_B - delta_B) * N_B
    dN_T/dt = S_T(t) - delta_T * N_T

with:
    lambda_B = (1/h_T + 1/b_0)^{-1}
    h_T = gamma * pi/(pi + Theta) * N_T_free
    N_T_free = N_T / (1 + D)
    D = tau_eng * gamma * sum_i Omega_i * pi_i/(pi_i + Theta) * N_B_i * dDG

State variables:
    N_A(t)          : scalar, antigen concentration
    pi(t, DG_i)     : per-cell pMHC for clone i
    N_B(t, DG_i)    : clone size for clone i
    N_T(t)          : scalar, total T cells
"""

import os

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass

N_Avg = 6.02214076e23  # Avogadro's number (molecules per mole)
# ============================================================
# Parameters
# ============================================================

@dataclass
class Parameters:
    """All model parameters in one place."""

    # --- Antigen ---
    N_A0: float = 1.0           # initial antigen concentration
    lambda_A: float = 1.0       # antigen growth rate (>0: replication, <0: decay)
    delta_A: float = 0.0        # antigen clearance rate (set to 0 for pure exponential)

    # --- pMHC ---
    k_on: float = 1.0           # antigen encounter rate (affinity-independent)
    delta_pi: float = 1.0       # pMHC surface turnover rate (fast)
    Theta: float = 1.0           # pMHC half-maximal threshold

    # --- Binding ---
    sigma: float = 1.0          # specificity parameter: psi(DG) = exp(-sigma * DG)
    beta_star: float = 2.0      # density-of-states exponent

    # --- T cells ---
    N_T0: float = 1000.0         # initial T cell pool
    delta_T: float = 0.0        # T cell turnover rate (0 = constant pool)
    gamma: float = 0.1         # T-B contact rate constant
    tau_eng: float = 0.1        # T-B engagement duration
    b_0: float = 0.3            # intrinsic cell division rate

    # --- B cells ---
    delta_B: float = 0.0        # B cell death rate

    # --- Repertoire ---
    DG_min: float = 0.0         # highest affinity (smallest DG)
    DG_max: float = 5.0        # lowest affinity in grid
    M: int = 20                # number of grid points
    Omega_0: float = 1.0        # density-of-states prefactor


# ============================================================
# Auxiliary functions
# ============================================================

def psi(DG, sigma):
    """Binding/internalization probability. Boltzmann gate."""
    return (1 + np.exp(DG))**-sigma


def Omega(DG, beta_star, Omega_0):
    """Density of states (number of clones per unit DG)."""
    return Omega_0 * np.exp(beta_star * DG)


def compute_demand(pi_vec, N_B_vec, Omega_vec, dDG, p):
    """
    Compute the demand function D(t).

    D = tau_eng * gamma * sum_i Omega_i * [pi_i/(pi_i + Theta)] * N_B_i * dDG
    """
    visibility = pi_vec / (pi_vec + p.Theta)
    integrand = Omega_vec * visibility * N_B_vec
    # return p.tau_eng * p.gamma * np.sum(integrand) * dDG
    return 0


def compute_lambda_B(pi_vec, N_T_free, p):
    """
    Compute division rate for each clone.

    lambda_B = (1/h_T + 1/b_0)^{-1}
    h_T = gamma * pi/(pi+Theta) * N_T_free
    """
    visibility = pi_vec / (pi_vec + p.Theta)
    h_T = p.gamma * visibility * N_T_free

    # Avoid division by zero when h_T = 0
    with np.errstate(divide='ignore', invalid='ignore'):
        inv_lambda = np.where(h_T > 0, (1.0 / h_T) + (1.0 / p.b_0), np.inf)
        lambda_B = np.where(inv_lambda < np.inf, 1.0 / inv_lambda, 0.0)

    return lambda_B


# ============================================================
# ODE system
# ============================================================

def pack_state(N_A, pi_vec, N_B_vec, N_T, M):
    """Pack all state variables into a single vector for the integrator."""
    y = np.zeros(2 * M + 2)
    y[0] = N_A
    y[1:M+1] = pi_vec
    y[M+1:2*M+1] = N_B_vec
    y[2*M+1] = N_T
    return y


def unpack_state(y, M):
    """Unpack the state vector."""
    N_A = y[0]
    pi_vec = y[1:M+1]
    N_B_vec = y[M+1:2*M+1]
    N_T = y[2*M+1]
    return N_A, pi_vec, N_B_vec, N_T


def rhs(t, y, p, DG_grid, Omega_vec, psi_vec, dDG):
    """Right-hand side of the coupled ODE system."""

    M = p.M
    N_A, pi_vec, N_B_vec, N_T = unpack_state(y, M)

    # Ensure non-negativity
    N_A = max(N_A, 0.0)
    pi_vec = np.maximum(pi_vec, 0.0)
    N_B_vec = np.maximum(N_B_vec, 0.0)
    N_T = max(N_T, 0.0)

    # --- Antigen ---
    # S_A = p.lambda_A * N_A if p.lambda_A > 0 else 0.0
    # dN_A = S_A - p.delta_A * N_A
    pb = (1 + (1e-9/(1e6*24*3600*np.exp(2.0*t)/N_Avg)))**(-1)  # or whatever dependence you intend
    dN_A = (p.lambda_A * (1 - pb) - 2*pb) * N_A - p.delta_A * N_A

    # --- Demand and free T cells ---
    D = compute_demand(pi_vec, N_B_vec, Omega_vec, dDG, p)
    N_T_free = N_T / (1.0 + D)

    # --- Division rate ---
    lambda_B = compute_lambda_B(pi_vec, N_T_free, p)

    # --- pMHC dynamics ---
    dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_B * pi_vec

    # --- Clone size dynamics ---
    dN_B = (lambda_B - p.delta_B) * N_B_vec

    # --- T cell dynamics ---
    dN_T = -p.delta_T * N_T  # constant pool if delta_T = 0

    return pack_state(dN_A, dpi, dN_B, dN_T, M)


# ============================================================
# Simulation
# ============================================================

def run_simulation(p=None, t_span=None, t_eval=None):
    """
    Run the mean-field simulation.

    Parameters
    ----------
    p : Parameters
        Model parameters. If None, uses defaults.
    t_span : tuple
        (t_start, t_end). If None, uses (0, 10).
    t_eval : array
        Time points for output. If None, uses 500 points.

    Returns
    -------
    result : dict with keys:
        't'       : time array
        'N_A'     : antigen concentration
        'pi'      : pMHC array (M x len(t))
        'N_B'     : clone size array (M x len(t))
        'N_T'     : T cell count
        'N_T_free': free T cells
        'D'       : demand function
        'lambda_B': division rate array (M x len(t))
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

    # --- Grid ---
    DG_grid = np.linspace(p.DG_min, p.DG_max, p.M)
    dDG = DG_grid[1] - DG_grid[0]

    # --- Precompute static arrays ---
    psi_vec = psi(DG_grid, p.sigma)
    Omega_vec = Omega(DG_grid, p.beta_star, p.Omega_0)

    # --- Initial conditions ---
    N_A_init = p.N_A0
    pi_init = np.zeros(p.M)  # no pMHC at t=0
    N_B_init = np.ones(p.M)  # one founder cell per clone
    N_T_init = p.N_T0

    y0 = pack_state(N_A_init, pi_init, N_B_init, N_T_init, p.M)

    # --- Integrate ---
    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, p, DG_grid, Omega_vec, psi_vec, dDG),
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

    # --- Unpack solution ---
    t = sol.t
    N_steps = len(t)

    N_A = sol.y[0, :]
    pi_arr = sol.y[1:p.M+1, :]           # shape (M, N_steps)
    N_B_arr = sol.y[p.M+1:2*p.M+1, :]    # shape (M, N_steps)
    N_T_arr = sol.y[2*p.M+1, :]

    # --- Compute derived quantities ---
    D_arr = np.zeros(N_steps)
    N_T_free_arr = np.zeros(N_steps)
    lambda_B_arr = np.zeros((p.M, N_steps))

    for j in range(N_steps):
        D_arr[j] = compute_demand(pi_arr[:, j], N_B_arr[:, j],
                                  Omega_vec, dDG, p)
        N_T_free_arr[j] = N_T_arr[j] / (1.0 + D_arr[j])
        lambda_B_arr[:, j] = compute_lambda_B(pi_arr[:, j],
                                              N_T_free_arr[j], p)

    return {
        't': t,
        'N_A': N_A,
        'pi': pi_arr,
        'N_B': N_B_arr,
        'N_T': N_T_arr,
        'N_T_free': N_T_free_arr,
        'D': D_arr,
        'lambda_B': lambda_B_arr,
        'DG_grid': DG_grid,
        'Omega': Omega_vec,
        'params': p,
    }


# ============================================================
# Visualization
# ============================================================

def plot_results(res):
    """Basic diagnostic plots."""

    t = res['t']
    DG = res['DG_grid']
    p = res['params']

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # --- (a) Antigen ---
    ax = axes[0, 0]
    ax.semilogy(t, res['N_A'])
    ax.semilogy(t, np.exp(p.lambda_A*t))
    ax.set_xlabel('Time')
    ax.set_ylabel('$N_A(t)$')
    ax.set_ylim(top = 1e13)
    ax.set_title('Antigen')

    # --- (b) Free T cells ---
    ax = axes[0, 1]
    ax.plot(t, res['N_T_free'], label='$N_T^{\\rm free}$')
    ax.plot(t, res['N_T'], '--', label='$N_T$', alpha=0.5)
    ax.set_xlabel('Time')
    ax.set_ylabel('T cells')
    ax.set_title('T cell pool')
    ax.legend()

    # --- (c) Demand ---
    ax = axes[0, 2]
    ax.semilogy(t, res['D'])
    ax.set_xlabel('Time')
    ax.set_ylabel('$D(t)$')
    ax.set_title('Demand function')

    # --- (d) Clone sizes at final time ---
    ax = axes[1, 0]
    N_B_final = res['N_B'][:, -1]
    ax.semilogy(DG, N_B_final, 'o-', label='simulation')
    ax.semilogy(DG, N_B_final[0]*np.exp(-p.sigma*(p.b_0/p.lambda_A + 1) * DG), 'r--', label='theory')
    ax.set_xlabel('$\\Delta G$')
    ax.set_ylabel('$N_B(t_{\\rm final}, \\Delta G)$')
    ax.set_title('Clone size distribution')
    ax.legend()

    # --- (e) Clone size heatmap ---
    ax = axes[1, 1]
    log_NB = np.log10(np.maximum(res['N_B'], 1.0))
    im = ax.pcolormesh(t, DG, log_NB, shading='auto', cmap='viridis')
    fig.colorbar(im, ax=ax, label='$\\log_{10} N_B$')
    ax.set_xlabel('Time')
    ax.set_ylabel('$\\Delta G$')
    ax.set_title('Clone size (log)')

    # Overlay theoretical front
    Gamma = p.gamma * p.N_T0 * p.k_on * p.N_A0 / (p.delta_pi * p.Theta)
    if Gamma > 0 and p.lambda_A > 0:
        DG_front = (p.lambda_A / p.sigma) * t \
                   - (1.0 / p.sigma) * np.log(p.lambda_A / Gamma)
        DG_front_clipped = np.clip(DG_front, p.DG_min, p.DG_max)
        valid = DG_front >= p.DG_min
        ax.plot(t[valid], DG_front_clipped[valid], 'r--', linewidth=2,
                label='$\\Delta G_{\\rm mf}(t)$')
        ax.legend()

    # --- (f) Division rate at selected times ---
    ax = axes[1, 2]
    n_snapshots = 5
    time_indices = np.linspace(0, len(DG)-1, n_snapshots, dtype=int)
    for idx in time_indices:
        ax.plot(t, res['lambda_B'][idx, :],
                label=f'$\\Delta G$={t[idx]:.1f}', alpha=0.8)
    ax.set_xlabel('Time')
    ax.set_ylabel('$\\lambda_B$')
    ax.set_title('Division rate')
    ax.legend(fontsize=8)

    plt.tight_layout()
    return fig


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_dynamics'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    p = Parameters(
        N_A0=1.0,
        lambda_A=6.0,
        delta_A=0.1,
        k_on=1e6*1e6*24*3600/N_Avg,
        delta_pi=1.0,
        Theta=100.0,
        sigma=2.0,
        beta_star=2.0,
        N_T0=1e5,
        delta_T=0.0,
        gamma=1.0,
        tau_eng=0.01,
        b_0=2.0,
        delta_B=0.0,
        DG_min=0.0,
        DG_max=6.0,
        M=10,
        Omega_0=1.0,
    )
    T = 15
    res = run_simulation(p=p, t_span=(0, T), t_eval=np.linspace(0, T, 500))

    fig = plot_results(res)
    fig.savefig(os.path.join(output_plot, 'ep_meanfield_results.png'), dpi=150, bbox_inches='tight')
    print(f"Saved: {os.path.join(output_plot, 'ep_meanfield_results.png')}")

    # Print some diagnostics
    print(f"\nFront velocity (theory): v = lambda_A / sigma = {p.lambda_A / p.sigma:.3f}")
    Gamma = p.gamma * p.N_T0 * p.k_on * p.N_A0 / (p.delta_pi * p.Theta)
    print(f"Gamma = {Gamma:.4f}")
    print(f"Final demand D = {res['D'][-1]:.3f}")
    print(f"Final N_T_free = {res['N_T_free'][-1]:.3f}")
    print(f"Max clone size = {res['N_B'][:, -1].max():.2e}")