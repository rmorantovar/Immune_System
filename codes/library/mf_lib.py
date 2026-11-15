"""
lib_mf.py -- Mean-field immune-activation models.

Refactored into a single modular engine:

    * `Parameters`         : all model parameters (dataclass).
    * repertoire builders  : grid / stochastic Delta-G sampling.
    * `StateLayout`        : automatic pack/unpack of named state blocks.
    * `MODELS`             : registry of the four models, selected by name:
                                'approx', 'complete', 'semicomplete', 'null'.
    * `run_simulation`     : ONE generic runner. Pick the model with
                             `run_simulation('semicomplete', ...)`.
    * analysis helpers     : potency, yield, N_B_tot, memory seeding, ...

The four models differ only in their state vector, their right-hand side and
their initial conditions -- all of which live in a small `ModelSpec`.  All the
plumbing (repertoire build, IC packing, `solve_ivp` call, unpacking, results
dict) is shared, so adding or editing a model touches one place only.

Plotting lives in `mf_plotting.py`; driving runs/sweeps lives in
`run_meanfield.py`.
"""

from dataclasses import dataclass, field, replace
from typing import Callable, List, Tuple, Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.stats import norm   # add near the top imports

# Physical constants (kept for parameter expressions in the driver).
N_Avg = 6.02214076e23          # Avogadro's number (molecules per mole)
k_BT = 1.380649e-23 * 293


# ============================================================
# Parameters
# ============================================================
@dataclass
class Parameters:
    """All model parameters in one place."""

    # --- Antigen ---
    N_A0: float = 1.0            # initial antigen concentration
    lambda_A: float = 1.0        # antigen growth rate (>0 replication, <0 decay)
    delta_A: float = 1.0         # antigen clearance rate (0 -> pure exponential)
    Z_c: float = 1e4             # neutralization-capacity threshold for feedback
    lambda_innate: float = 1.8   # innate immune response rate (external drive)
    threshold_innate: float = 2e3  # threshold for innate-response activation

    # --- pMHC ---
    k_on: float = 1.0            # antigen encounter rate (affinity-independent)
    delta_pi: float = 1.0        # pMHC surface turnover rate (fast)
    pi_star: float = 1.0         # pMHC half-maximal threshold
    hill: float = 3.0            # Hill coefficient for activation sharpness

    # --- Binding repertoire ---
    eta: float = 1.0             # specificity: psi(DG) = exp(-eta * DG)
    beta_star: float = 2.0       # density-of-states exponent
    sigma_E: float = 4.0         # Gaussian width in the density of states
    DG_min: float = 0.0          # highest affinity (smallest DG)
    DG_max: float = 5.0          # lowest affinity in grid
    M: int = 20                  # number of grid points
    omega_0: float = 1.0         # density-of-states prefactor

    # --- T cells ---
    N_T0: float = 1000.0         # initial T-cell pool
    delta_T: float = 0.0         # T-cell turnover rate (0 = constant pool)
    h0: float = 0.1              # T-B contact rate constant
    tau_eng: float = 0.1         # T-B engagement duration
    b0: float = 0.3              # intrinsic cell-division rate
    K_T: float = 10.0            # saturation constant for T-cell limitation
    Tcell_growth_factor: float = 1.0  # T-cell expansion factor upon activation

    # --- B cells ---
    delta_B: float = 0.0         # B-cell death rate
    n_mem: float = 1e4           # number of memory cells to sample from primary

    # --- Simulation options ---
    T_lim: bool = False          # include T-cell limitation in the demand function
    memory: bool = False         # use memory-like initial conditions for N_B
    h0_mem: Optional[float] = None  # memory-population h0 in the 'boost' models (None -> h0)
    t_switch: float = 6.0        # 'two_stage' models: antigen loading for t < t_switch, expansion after
    h_max: float = 1e3           # 'two_stage' models: max T-cell engagement rate per B cell (per day)

    # --- Explicit clones (mode='clones') ---
    clone_DG: Optional[tuple] = None   # DG of each clone, e.g. (0.0, 1.5)
    clone_N0: Optional[tuple] = None   # initial size of each clone (default 1 each)


# ============================================================
# Repertoire
# ============================================================
def psi(DG, eta):
    """Binding / internalization probability. Boltzmann gate."""
    return (1 + np.exp(DG + 0.5)) ** -eta


def Omega_0(DG, beta_star, omega_0, sigma_E):
    """Density of states (number of clones per unit DG)."""
    return omega_0 * np.exp(beta_star * DG - DG ** 2 / (2 * sigma_E ** 2))


def G1(pi_vec, p):
    """B-T cell engagement probability (linear regime)."""
    return pi_vec


def build_repertoire_grid(p):
    """Uniform grid in DG with Omega weighting. Returns (DG, psi, weights, M)."""
    DG_arr = np.linspace(p.DG_min, p.DG_max, p.M)
    dDG = DG_arr[1] - DG_arr[0]
    psi_arr = psi(DG_arr, p.eta)
    Omega_arr = Omega_0(DG_arr, p.beta_star, p.omega_0, p.sigma_E)
    weights = Omega_arr * dDG              # each bin represents this many clones
    return DG_arr, psi_arr, weights, p.M


def build_repertoire_stochastic(p, DG_max_sim=None, seed=None):
    """
    Sample individual clone energies from Omega_0(DG) on [DG_floor, DG_max_sim].
    Returns (DG, psi, weights=ones, L_sim).
    """
    if DG_max_sim is None:
        DG_max_sim = p.DG_max

    rng = np.random.default_rng(seed)
    n_extra = 3
    DG_floor = p.DG_min - n_extra / p.beta_star   # extend by n_extra expected clones

    L_sim = int(p.omega_0 / p.beta_star * (
        np.exp(p.beta_star * DG_max_sim) - np.exp(p.beta_star * DG_floor)
    ))
    if L_sim <= 0:
        raise ValueError(f"No clones in range [{p.DG_min}, {DG_max_sim}]. "
                         f"Check omega_0 and beta_star.")

    # Inverse-CDF sampling of the truncated exponential P(DG) ~ exp(beta_star*DG).
    u = rng.uniform(0, 1, size=L_sim)
    exp_min = np.exp(p.beta_star * DG_floor)
    exp_max = np.exp(p.beta_star * DG_max_sim)
    DG_arr = (1.0 / p.beta_star) * np.log(u * (exp_max - exp_min) + exp_min)
    DG_arr = np.sort(DG_arr)

    psi_arr = psi(DG_arr, p.eta)
    weights = np.ones(L_sim)
    return DG_arr, psi_arr, weights, L_sim


def build_repertoire_clones(p):
    """
    A handful of explicit clones at p.clone_DG (e.g. a two-clone competition).
    Returns (DG, psi, weights=ones, M=len(clone_DG)).
    """
    if p.clone_DG is None:
        raise ValueError("mode='clones' needs p.clone_DG (one DG per clone).")
    DG_arr = np.asarray(p.clone_DG, dtype=float)
    psi_arr = psi(DG_arr, p.eta)
    weights = np.ones(len(DG_arr))
    return DG_arr, psi_arr, weights, len(DG_arr)


def build_context(p, mode='grid', DG_max_sim=None, seed=None):
    """Everything the RHS/IC need that only depends on the repertoire + params."""
    if mode == 'grid':
        DG, psi_arr, weights, M = build_repertoire_grid(p)
    elif mode == 'stochastic':
        DG, psi_arr, weights, M = build_repertoire_stochastic(p, DG_max_sim, seed)
    elif mode == 'clones':
        DG, psi_arr, weights, M = build_repertoire_clones(p)
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'grid', 'stochastic' or 'clones'.")
    # naive initial clone sizes (all ones unless explicit clones say otherwise)
    if mode == 'clones' and p.clone_N0 is not None:
        N0 = np.asarray(p.clone_N0, dtype=float)
        if N0.shape != DG.shape:
            raise ValueError("clone_N0 and clone_DG must have the same length.")
    else:
        N0 = np.ones(M)
    return dict(DG=DG, psi=psi_arr, weights=weights,
                neut=weights * np.exp(-DG), M=M, mode=mode, N0=N0)


# ============================================================
# Generic state layout (automatic pack / unpack)
# ============================================================
class StateLayout:
    """
    Map named state blocks onto a flat vector for solve_ivp.

    blocks : list of (name, kind) where kind is 'scalar' or 'vector'.
             'scalar' occupies 1 slot, 'vector' occupies M slots.
    """

    def __init__(self, blocks: List[Tuple[str, str]], M: int):
        self.blocks = blocks
        self.M = M
        self.slices = {}
        i = 0
        for name, kind in blocks:
            n = 1 if kind == 'scalar' else M
            self.slices[name] = (i, i + n, kind)
            i += n
        self.size = i

    def pack(self, d: dict) -> np.ndarray:
        y = np.zeros(self.size)
        for name, (a, b, kind) in self.slices.items():
            if kind == 'scalar':
                y[a] = d[name]
            else:
                y[a:b] = d[name]
        return y

    def unpack(self, y: np.ndarray) -> dict:
        """
        Works for a 1-D state vector (returns floats / (M,) arrays) and for the
        2-D solution array of shape (size, N_t) (returns (N_t,) / (M, N_t)).
        """
        out = {}
        for name, (a, b, kind) in self.slices.items():
            out[name] = y[a] if kind == 'scalar' else y[a:b]
        return out


# ============================================================
# Model specifications
# ============================================================
@dataclass
class ModelSpec:
    """One model = state layout + rhs + initial conditions + solver defaults."""
    name: str
    blocks: List[Tuple[str, str]]
    rhs: Callable            # rhs(t, state_dict, p, ctx) -> dict of derivatives
    init: Callable           # init(p, ctx, memory_seed) -> dict of states
    solver: dict = field(default_factory=lambda: dict(
        method='LSODA', rtol=1e-4, atol=1e-5))
    postprocess: Optional[Callable] = None   # postprocess(res, p, ctx) -> None
    breaks: Optional[Callable] = None        # breaks(p) -> times where the rhs switches law;
                                             # integrated piecewise, ctx['segment'] = segment index

    def layout(self, M: int) -> StateLayout:
        return StateLayout(self.blocks, M)


# --- shared helpers used by several models --------------------------------
def compute_demand(pi_vec, N_B_vec, weights, p):
    """Demand function D(t) (used by the 'approx' model)."""
    return np.sum(weights * G1(pi_vec, p) * N_B_vec)


def compute_h_T(pi_vec, N_T, N_B_vec, weights, p):
    D = compute_demand(pi_vec, N_B_vec, weights, p)
    N_T_free = N_T / (1.0 + D / p.K_T) if p.T_lim else N_T
    return p.h0 * G1(pi_vec, p) * N_T_free


def compute_h_B(pi_vec, N_B_vec, weights, p):
    D = compute_demand(pi_vec, N_B_vec, weights, p)
    return p.h0 * D / (p.K_T + D)


def compute_lambda_B(pi_vec, N_B_vec, N_T, weights, p):
    h_T = compute_h_T(pi_vec, N_T, N_B_vec, weights, p)
    with np.errstate(divide='ignore', invalid='ignore'):
        inv = np.where(h_T > 0, (1.0 / h_T) + (1.0 / p.b0), np.inf)
        return np.where(inv < np.inf, 1.0 / inv, 0.0)


def compute_lambda_T(pi_vec, N_B_vec, weights, p):
    h_B = compute_h_B(pi_vec, N_B_vec, weights, p)
    with np.errstate(divide='ignore', invalid='ignore'):
        inv = np.where(h_B > 0, (1.0 / h_B) + (1.0 / p.b0), np.inf)
        return np.where(inv < np.inf, 1.0 / inv, 0.0)


# -------------------------------------------------------------------------
# approx model
# -------------------------------------------------------------------------
def _init_approx(p, ctx, memory_seed=None):
    M, DG = ctx['M'], ctx['DG']
    if p.memory:
        N_B = 1e2 * np.exp(-p.eta * (p.b0 / p.lambda_A + 1) * DG)
    else:
        N_B = np.ones(M)
    return dict(N_A=p.N_A0, pi=np.zeros(M), N_B=N_B, N_T=p.N_T0)


def _rhs_approx(t, s, p, ctx):
    psi_vec, weights = ctx['psi'], ctx['weights']
    N_A = max(s['N_A'], 0.0)
    pi_vec = np.maximum(s['pi'], 0.0)
    N_B_vec = np.maximum(s['N_B'], 0.0)
    N_T = max(s['N_T'], 0.0)

    pb = (1 + (1e-9 / (1e6 * 24 * 3600 * np.exp(2.0 * t) / N_Avg))) ** (-1)
    dN_A = (p.lambda_A * (1 - pb) - p.delta_A * pb) * N_A - 0.01 * N_A

    lambda_B = compute_lambda_B(pi_vec, N_B_vec, N_T, weights, p)
    dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_B * pi_vec
    dN_B = (lambda_B - p.delta_B) * N_B_vec
    lambda_T = compute_lambda_T(pi_vec, N_B_vec, weights, p)
    dN_T = lambda_T * N_T - p.delta_T * N_T
    return dict(N_A=dN_A, pi=dpi, N_B=dN_B, N_T=dN_T)


def _post_approx(res, p, ctx):
    """Derived quantities the diagnostics plots for the approx model expect."""
    t = res['t']
    weights = ctx['weights']
    M, N_steps = ctx['M'], len(t)
    pi_arr, N_B_arr, N_T_arr = res['pi'], res['N_B'], res['N_T']

    D = np.zeros(N_steps)
    N_T_free = np.zeros(N_steps)
    lambda_B = np.zeros((M, N_steps))
    lambda_T = np.zeros(N_steps)
    h_T = np.zeros((M, N_steps))
    h_B = np.zeros(N_steps)
    for j in range(N_steps):
        D[j] = compute_demand(pi_arr[:, j], N_B_arr[:, j], weights, p)
        N_T_free[j] = N_T_arr[j] / (1.0 + D[j] / p.K_T) if p.T_lim else N_T_arr[j]
        h_T[:, j] = compute_h_T(pi_arr[:, j], N_T_arr[j], N_B_arr[:, j], weights, p)
        h_B[j] = compute_h_B(pi_arr[:, j], N_B_arr[:, j], weights, p)
        lambda_B[:, j] = compute_lambda_B(pi_arr[:, j], N_B_arr[:, j], N_T_arr[j], weights, p)
        lambda_T[j] = compute_lambda_T(pi_arr[:, j], N_B_arr[:, j], weights, p)
    res.update(D=D, N_T_free=N_T_free, lambda_B=lambda_B,
               lambda_T=lambda_T, h_T=h_T, h_B=h_B)


# -------------------------------------------------------------------------
# complete model (explicit T-B conjugate N_BT)
# -------------------------------------------------------------------------
def _init_complete(p, ctx, memory_seed=None):
    M, DG = ctx['M'], ctx['DG']
    if p.memory:
        N_Bo = 1e2 * np.exp(-p.eta * (p.b0 / p.lambda_A + 1) * DG)
    else:
        N_Bo = np.ones(M)
    return dict(N_A=p.N_A0, pi=np.zeros(M), N_Bo=N_Bo,
                N_BT=np.zeros(M), N_Ba=np.zeros(M), N_To=p.N_T0, N_Ta=0.0)


def _rhs_complete(t, s, p, ctx):
    psi_vec, weights = ctx['psi'], ctx['weights']
    N_A = max(s['N_A'], 0.0)
    pi_vec = np.maximum(s['pi'], 0.0)
    N_Bo = np.maximum(s['N_Bo'], 0.0)
    N_BT = np.maximum(s['N_BT'], 0.0)
    N_Ba = np.maximum(s['N_Ba'], 0.0)
    N_To = max(s['N_To'], 0.0)
    N_Ta = max(s['N_Ta'], 0.0)

    pb = (1 + (1e-9 / (1e6 * 24 * 3600 * np.exp(p.b0 * t) / N_Avg))) ** (-1)
    dN_A = (p.lambda_A * (1 - pb) - p.delta_A * pb) * N_A - 0.01 * N_A

    N_B_tot = N_Bo + N_BT + N_Ba
    with np.errstate(divide='ignore', invalid='ignore'):
        lambda_eff = np.where(N_B_tot > 0, p.b0 * N_Ba / N_B_tot, 0.0)
    dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_eff * pi_vec

    dN_Bo = (- p.k_on * pi_vec ** p.hill * N_Bo * N_To / (N_To + p.K_T)
             + 2 * p.b0 * N_Ba - p.delta_B * N_Bo)
    dN_BT = (p.k_on * pi_vec ** p.hill * N_Bo * N_To / (N_To + p.K_T)
             - p.h0 * N_BT)
    dN_Ba = p.h0 * N_BT - p.b0 * N_Ba - p.delta_B * N_Ba
    dN_To = (- p.k_on * np.sum(weights * pi_vec ** p.hill * N_Bo)
             * N_To / (N_To + p.K_T) + 2 * p.b0 * N_Ta - p.delta_T * N_To)
    dN_Ta = p.h0 * np.sum(weights * N_BT) - p.b0 * N_Ta - p.delta_T * N_Ta
    return dict(N_A=dN_A, pi=dpi, N_Bo=dN_Bo, N_BT=dN_BT,
                N_Ba=dN_Ba, N_To=dN_To, N_Ta=dN_Ta)


# -------------------------------------------------------------------------
# shared antigen block for semicomplete / null
# -------------------------------------------------------------------------
def _dN_A_with_feedback(t, N_A, p, Z):
    """
    Antigen ODE. When p.memory is off, antigen is suppressed by an external
    innate drive; when on, by the neutralization capacity Z(t).
    """
    if not p.memory:
        innate_proxy = np.exp(p.lambda_innate * t)
        frac = (1.0 + p.threshold_innate / (innate_proxy + 1e-30)) ** (-1)
    else:
        frac = (1.0 + p.Z_c / (max(Z, 0.0) + 1e-30)) ** (-1)
    return (p.lambda_A * (1 - frac) - p.delta_A * frac) * N_A - 0.01 * N_A


def _memory_ic(p, ctx, memory_seed):
    """Shared naive/memory initial B-cell clone sizes (semicomplete & null)."""
    M, DG = ctx['M'], ctx['DG']
    if p.memory:
        if memory_seed is not None:
            return memory_seed.copy()
        alpha = p.eta * (p.b0 / p.lambda_A + 1)
        return 1e4 * np.exp(-alpha * DG)          # analytic fallback
    return ctx['N0'].copy()


# -------------------------------------------------------------------------
# semicomplete model (implicit conjugate; pMHC + T cells + activation latch a)
# -------------------------------------------------------------------------
def _init_semicomplete(p, ctx, memory_seed=None):
    M = ctx['M']
    return dict(N_A=p.N_A0, pi=np.zeros(M),
                N_Bo=_memory_ic(p, ctx, memory_seed), N_Ba=np.zeros(M),
                a=np.zeros(M), N_To=p.N_T0, N_Ta=0.0)


def _rhs_semicomplete(t, s, p, ctx):
    psi_vec, weights, neut = ctx['psi'], ctx['weights'], ctx['neut']
    N_A = max(s['N_A'], 0.0)
    pi_vec = np.maximum(s['pi'], 0.0)
    N_Bo = np.maximum(s['N_Bo'], 0.0)
    N_Ba = np.maximum(s['N_Ba'], 0.0)
    a_vec = s['a']
    N_To = max(s['N_To'], 0.0)
    N_Ta = max(s['N_Ta'], 0.0)

    pi_star = (p.b0 / p.h0) ** (1 / p.hill)
    g = pi_vec ** p.hill / (pi_vec ** p.hill + pi_star ** p.hill)
    da = (10.0 * p.b0) * (1.0 - a_vec) * g               # ratchet: up only

    N_B_tot = N_Bo + N_Ba
    with np.errstate(divide='ignore', invalid='ignore'):
        lambda_eff = np.where(N_B_tot > 0, p.b0 * N_Ba / N_B_tot, 0.0)

    Z = np.sum(neut * (N_Bo + N_Ba) * a_vec)
    dN_A = _dN_A_with_feedback(t, N_A, p, Z)

    dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_eff * pi_vec
    dN_Bo = (- p.h0 * pi_vec ** p.hill * N_Bo * N_To / (N_To + p.K_T)
             + 2 * p.b0 * N_Ba - p.delta_B * N_Bo)
    dN_Ba = (p.h0 * pi_vec ** p.hill * N_Bo * N_To / (N_To + p.K_T)
             - p.b0 * N_Ba - p.delta_B * N_Ba)
    dN_To = (- p.h0 * np.sum(weights * pi_vec ** p.hill * N_Bo)
             * N_To / (N_To + p.K_T)
             + p.Tcell_growth_factor * p.b0 * N_Ta - p.delta_T * N_To)
    dN_Ta = (p.h0 * np.sum(weights * pi_vec ** p.hill * N_Bo)
             * N_To / (N_To + p.K_T) - p.b0 * N_Ta - p.delta_T * N_Ta)
    return dict(N_A=dN_A, pi=dpi, N_Bo=dN_Bo, N_Ba=dN_Ba,
                a=da, N_To=dN_To, N_Ta=dN_Ta)


# -------------------------------------------------------------------------
# boost model: semicomplete with two B-cell populations on the same DG grid,
# naive (suffix '') and memory ('_m'), sharing antigen and the T-cell pool.
# Each population has its own pMHC load, activation latch and h0
# (naive: p.h0, memory: p.h0_mem). With an empty memory seed it reduces
# exactly to the semicomplete model.
# -------------------------------------------------------------------------
BOOST_POPS = (('', 'h0'), ('_m', 'h0_mem'))


def _init_boost(p, ctx, memory_seed=None):
    M = ctx['M']
    N_Bo_m = np.zeros(M) if memory_seed is None else memory_seed.copy()
    # Bins without memory stay frozen: integrating an exactly-empty population
    # makes N_Ba/N_B_tot jump between 0 and 1 and stalls the solver.
    ctx['live_m'] = N_Bo_m > 0
    return dict(N_A=p.N_A0,
                pi=np.zeros(M), N_Bo=ctx['N0'].copy(), N_Ba=np.zeros(M), a=np.zeros(M),
                pi_m=np.zeros(M), N_Bo_m=N_Bo_m, N_Ba_m=np.zeros(M), a_m=np.zeros(M),
                N_To=p.N_T0, N_Ta=0.0)


def _rhs_boost(t, s, p, ctx):
    psi_vec, weights, neut = ctx['psi'], ctx['weights'], ctx['neut']
    N_A = max(s['N_A'], 0.0)
    N_To = max(s['N_To'], 0.0)
    N_Ta = max(s['N_Ta'], 0.0)
    T_free = N_To / (N_To + p.K_T)

    d, flux_T, Z = {}, 0.0, 0.0
    for sfx, h0_name in BOOST_POPS:
        h0 = getattr(p, h0_name)
        h0 = p.h0 if h0 is None else h0
        pi_vec = np.maximum(s['pi' + sfx], 0.0)
        N_Bo = np.maximum(s['N_Bo' + sfx], 0.0)
        N_Ba = np.maximum(s['N_Ba' + sfx], 0.0)
        a_vec = s['a' + sfx]

        pi_star = (p.b0 / h0) ** (1 / p.hill)
        g = pi_vec ** p.hill / (pi_vec ** p.hill + pi_star ** p.hill)
        N_B_tot = N_Bo + N_Ba
        with np.errstate(divide='ignore', invalid='ignore'):
            lambda_eff = np.where(N_B_tot > 0, p.b0 * N_Ba / N_B_tot, 0.0)
        act = h0 * pi_vec ** p.hill * N_Bo * T_free          # per-clone activation flux
        live = ctx['live_m'] if sfx == '_m' else 1.0
        act = act * live

        d['a' + sfx] = (10.0 * p.b0) * (1.0 - a_vec) * g * live
        d['pi' + sfx] = (p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_eff * pi_vec) * live
        d['N_Bo' + sfx] = (- act + 2 * p.b0 * N_Ba - p.delta_B * N_Bo) * live
        d['N_Ba' + sfx] = (act - p.b0 * N_Ba - p.delta_B * N_Ba) * live
        flux_T += np.sum(weights * act)
        Z += np.sum(neut * N_B_tot * a_vec)

    d['N_A'] = _dN_A_with_feedback(t, N_A, p, Z)
    d['N_To'] = - flux_T + p.Tcell_growth_factor * p.b0 * N_Ta - p.delta_T * N_To
    d['N_Ta'] = flux_T - p.b0 * N_Ta - p.delta_T * N_Ta
    return d


# -------------------------------------------------------------------------
# two-stage models: antigen loading and T-cell help are separated in time.
#   t <  t_switch : B cells take up antigen (pMHC load pi grows), no division;
#   t >= t_switch : no antigen uptake; B cells engage T cells, divide, and the
#                   loaded pMHC is diluted by division until pi < pi_c.
# T cells form a finite pool of N_T0 cells: each conjugate (N_BT) occupies one
# T cell for tau_eng, so engagement slows as T cells run out. This is the only
# coupling between clones (and between naive and memory in 'two_stage_boost').
# h(pi) = h0 pi^hill / (1 + h0 pi^hill / h_max): a B cell cannot engage faster than
# it meets T cells; without the cap the loaded pMHC right after t_switch gives
# rates ~1e8/day and an intractably stiff system.
# Per population: resting N_Bo -(h(pi) T_free/N_T0)-> engaged N_BT
#                 -(1/tau_eng)-> licensed N_Ba -(b0)-> 2 resting N_Bo.
# -------------------------------------------------------------------------
# Engagement right after t_switch is very fast (h0 pi^hill >> b0): stiff, so BDF;
# atol well below the smallest per-clone memory sizes (~1e-4 cells).
TWO_STAGE_SOLVER = dict(method='BDF', rtol=1e-4, atol=1e-10)
TWO_STAGE_POPS = {'two_stage': (('', 'h0'),),
                  'two_stage_boost': (('', 'h0'), ('_m', 'h0_mem'))}


def _two_stage_pop_ic(p, ctx, N_Bo, sfx=''):
    M = ctx['M']
    return {'pi' + sfx: np.zeros(M), 'N_Bo' + sfx: N_Bo, 'N_BT' + sfx: np.zeros(M),
            'N_Ba' + sfx: np.zeros(M), 'a' + sfx: np.zeros(M)}


def _init_two_stage(p, ctx, memory_seed=None):
    ctx['pops'] = TWO_STAGE_POPS['two_stage']
    return dict(N_A=p.N_A0, **_two_stage_pop_ic(p, ctx, _memory_ic(p, ctx, memory_seed)))


def _init_two_stage_boost(p, ctx, memory_seed=None):
    ctx['pops'] = TWO_STAGE_POPS['two_stage_boost']
    N_Bo_m = np.zeros(ctx['M']) if memory_seed is None else memory_seed.copy()
    ctx['live_m'] = N_Bo_m > 0            # bins without memory stay frozen (see _init_boost)
    return dict(N_A=p.N_A0, **_two_stage_pop_ic(p, ctx, ctx['N0'].copy()),
                **_two_stage_pop_ic(p, ctx, N_Bo_m, '_m'))


def _rhs_two_stage(t, s, p, ctx):
    psi_vec, weights, neut = ctx['psi'], ctx['weights'], ctx['neut']
    N_A = max(s['N_A'], 0.0)
    loading = ctx['segment'] == 0 if 'segment' in ctx else t < p.t_switch
    engaged = sum(np.sum(weights * np.maximum(s['N_BT' + sfx], 0.0)) for sfx, _ in ctx['pops'])
    T_free = max(p.N_T0 - engaged, 0.0) / p.N_T0            # fraction of T cells free

    d, Z = {}, 0.0
    for sfx, h0_name in ctx['pops']:
        h0 = getattr(p, h0_name)
        h0 = p.h0 if h0 is None else h0
        live = ctx['live_m'] if sfx == '_m' else 1.0
        pi_vec = np.maximum(s['pi' + sfx], 0.0)
        N_Bo = np.maximum(s['N_Bo' + sfx], 0.0)
        N_BT = np.maximum(s['N_BT' + sfx], 0.0)
        N_Ba = np.maximum(s['N_Ba' + sfx], 0.0)
        a_vec = s['a' + sfx]

        N_B_tot = N_Bo + N_BT + N_Ba
        pi_star = (p.b0 / h0) ** (1 / p.hill)
        g = pi_vec ** p.hill / (pi_vec ** p.hill + pi_star ** p.hill)
        if loading:
            dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec
            eng = rel = div = 0.0
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                lambda_eff = np.where(N_B_tot > 0, p.b0 * N_Ba / N_B_tot, 0.0)
            dpi = - p.delta_pi * pi_vec - lambda_eff * pi_vec   # dilution by division
            h = h0 * pi_vec ** p.hill
            eng = h / (1.0 + h / p.h_max) * N_Bo * T_free      # engagement capped at h_max
            rel = N_BT / p.tau_eng
            div = p.b0 * N_Ba

        d['pi' + sfx] = dpi * live
        d['N_Bo' + sfx] = (- eng + 2 * div - p.delta_B * N_Bo) * live
        d['N_BT' + sfx] = (eng - rel - p.delta_B * N_BT) * live
        d['N_Ba' + sfx] = (rel - div - p.delta_B * N_Ba) * live
        d['a' + sfx] = (10.0 * p.b0) * (1.0 - a_vec) * g * live
        Z += np.sum(neut * N_B_tot * a_vec)

    d['N_A'] = _dN_A_with_feedback(t, N_A, p, Z)
    return d


def _post_two_stage(res, p, ctx):
    """Total cells per clone (resting + engaged + licensed) and the free T-cell fraction."""
    engaged = 0.0
    for sfx, _ in ctx['pops']:
        res['N_B' + sfx] = res['N_Bo' + sfx] + res['N_BT' + sfx] + res['N_Ba' + sfx]
        engaged = engaged + np.sum(ctx['weights'][:, None] * res['N_BT' + sfx], axis=0)
    res['T_free'] = np.maximum(p.N_T0 - engaged, 0.0) / p.N_T0


# -------------------------------------------------------------------------
# null model (no pMHC / no T cells; antigen drives activation directly)
# -------------------------------------------------------------------------
def _init_null(p, ctx, memory_seed=None):
    M = ctx['M']
    return dict(N_A=p.N_A0, N_Bo=_memory_ic(p, ctx, memory_seed),
                N_Ba=np.zeros(M), a=np.zeros(M))


def _rhs_null(t, s, p, ctx):
    psi_vec, neut = ctx['psi'], ctx['neut']
    N_A = max(s['N_A'], 0.0)
    N_Bo = np.maximum(s['N_Bo'], 0.0)
    N_Ba = np.maximum(s['N_Ba'], 0.0)
    a_vec = np.clip(s['a'], 0.0, 1.0)

    u_p = p.k_on * psi_vec * N_A            # antigen processing -> activation drive
    u_p_c = 1.0                             # analog of pi_c = b0/h0
    g = u_p ** p.hill / (u_p ** p.hill + u_p_c ** p.hill)
    da = 1.0 * (1.0 - a_vec) * g            # activation latch

    Z = np.sum(neut * (N_Bo + N_Ba) * a_vec)
    dN_A = _dN_A_with_feedback(t, N_A, p, Z)

    dN_Bo = - p.h0 * u_p * N_Bo + 2 * p.b0 * N_Ba - p.delta_B * N_Bo
    dN_Ba = p.h0 * u_p * N_Bo - p.b0 * N_Ba - p.delta_B * N_Ba
    return dict(N_A=dN_A, N_Bo=dN_Bo, N_Ba=dN_Ba, a=da)


# -------------------------------------------------------------------------
# registry
# -------------------------------------------------------------------------
MODELS = {
    'approx': ModelSpec(
        name='approx',
        blocks=[('N_A', 'scalar'), ('pi', 'vector'),
                ('N_B', 'vector'), ('N_T', 'scalar')],
        rhs=_rhs_approx, init=_init_approx,
        solver=dict(method='RK45', rtol=1e-8, atol=1e-10, max_step=0.01),
        postprocess=_post_approx,
    ),
    'complete': ModelSpec(
        name='complete',
        blocks=[('N_A', 'scalar'), ('pi', 'vector'), ('N_Bo', 'vector'),
                ('N_BT', 'vector'), ('N_Ba', 'vector'),
                ('N_To', 'scalar'), ('N_Ta', 'scalar')],
        rhs=_rhs_complete, init=_init_complete,
    ),
    'semicomplete': ModelSpec(
        name='semicomplete',
        blocks=[('N_A', 'scalar'), ('pi', 'vector'), ('N_Bo', 'vector'),
                ('N_Ba', 'vector'), ('a', 'vector'),
                ('N_To', 'scalar'), ('N_Ta', 'scalar')],
        rhs=_rhs_semicomplete, init=_init_semicomplete,
    ),
    'null': ModelSpec(
        name='null',
        blocks=[('N_A', 'scalar'), ('N_Bo', 'vector'),
                ('N_Ba', 'vector'), ('a', 'vector')],
        rhs=_rhs_null, init=_init_null,
    ),
    'boost': ModelSpec(
        name='boost',
        blocks=[('N_A', 'scalar'),
                ('pi', 'vector'), ('N_Bo', 'vector'), ('N_Ba', 'vector'), ('a', 'vector'),
                ('pi_m', 'vector'), ('N_Bo_m', 'vector'), ('N_Ba_m', 'vector'), ('a_m', 'vector'),
                ('N_To', 'scalar'), ('N_Ta', 'scalar')],
        rhs=_rhs_boost, init=_init_boost,
    ),
    'two_stage': ModelSpec(
        name='two_stage',
        blocks=[('N_A', 'scalar'), ('pi', 'vector'), ('N_Bo', 'vector'),
                ('N_BT', 'vector'), ('N_Ba', 'vector'), ('a', 'vector')],
        rhs=_rhs_two_stage, init=_init_two_stage, postprocess=_post_two_stage,
        breaks=lambda p: [p.t_switch],
        solver=TWO_STAGE_SOLVER,
    ),
    'two_stage_boost': ModelSpec(
        name='two_stage_boost',
        blocks=[('N_A', 'scalar'),
                ('pi', 'vector'), ('N_Bo', 'vector'), ('N_BT', 'vector'), ('N_Ba', 'vector'), ('a', 'vector'),
                ('pi_m', 'vector'), ('N_Bo_m', 'vector'), ('N_BT_m', 'vector'), ('N_Ba_m', 'vector'),
                ('a_m', 'vector')],
        rhs=_rhs_two_stage, init=_init_two_stage_boost, postprocess=_post_two_stage,
        breaks=lambda p: [p.t_switch],
        solver=TWO_STAGE_SOLVER,
    ),
}


# ============================================================
# Generic runner
# ============================================================
def run_simulation(model='semicomplete', p=None, t_span=None, t_eval=None,
                   mode='grid', DG_max_sim=None, seed=None, memory_seed=None,
                   drive=None, **solver_overrides):
    """
    Run any registered model.

    Parameters
    ----------
    model : 'approx' | 'complete' | 'semicomplete' | 'null' | 'boost'
    p : Parameters
    t_span : (t0, t1)
    t_eval : output times (default 500 points over t_span)
    mode : 'grid' | 'stochastic' | 'clones' (explicit p.clone_DG / p.clone_N0)
    memory_seed : optional per-clone N_Bo initial condition (memory phase;
                  the memory population in the 'boost' model)
    drive : optional dict(t=..., N_A=...). N_A is then not integrated but
            replayed from this trajectory (cubic spline in log N_A), so the
            response has no feedback on antigen. None -> unchanged behaviour.
    solver_overrides : override the model's default solve_ivp kwargs.

    Returns
    -------
    res : dict with the state arrays plus 't', 'DG', 'weights', 'params',
          'mode', 'M', 'model'.
    """
    if model not in MODELS:
        raise ValueError(f"Unknown model '{model}'. Choose from {list(MODELS)}.")
    spec = MODELS[model]

    if p is None:
        p = Parameters()
    if t_span is None:
        t_span = (0.0, 10.0)
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500)

    ctx = build_context(p, mode=mode, DG_max_sim=DG_max_sim, seed=seed)
    return _integrate(spec, p, ctx, spec.init(p, ctx, memory_seed), t_span, t_eval,
                      drive, mode, model, solver_overrides)


def _integrate(spec, p, ctx, init_state, t_span, t_eval, drive, mode, model, solver_overrides,
               strict=False):
    """Integrate one ModelSpec from init_state; shared by run_simulation and run_response.
    strict: raise if the solver fails (run_simulation only warns, as before)."""
    layout = spec.layout(ctx['M'])
    y0 = layout.pack(init_state)
    solver_kw = dict(spec.solver)
    solver_kw.update(solver_overrides)

    if drive is None:
        def fun(t, y):
            return layout.pack(spec.rhs(t, layout.unpack(y), p, ctx))
    else:
        log_NA = np.log(np.maximum(np.asarray(drive['N_A'], dtype=float), 1e-300))
        N_A_of = lambda t, f=CubicSpline(drive['t'], log_NA): np.exp(f(t))   # smooth: kinks stall the solver
        i_NA = layout.slices['N_A'][0]

        def fun(t, y):
            s = layout.unpack(y)
            s['N_A'] = N_A_of(t)
            dy = layout.pack(spec.rhs(t, s, p, ctx))
            dy[i_NA] = 0.0                     # N_A is read from the drive, not integrated
            return dy

    breaks = [] if spec.breaks is None else [b for b in spec.breaks(p) if t_span[0] < b < t_span[1]]
    if not breaks:
        sol = solve_ivp(fun, t_span, y0, t_eval=t_eval, **solver_kw)
        if not sol.success:
            if strict:
                raise RuntimeError(f"integration failed ({model}): {sol.message}")
            print(f"Warning: integration failed ({model}): {sol.message}")
        t_out, y_out = sol.t, sol.y
    else:
        # Integrate each segment separately so the solver never steps across a switch.
        edges = [t_span[0], *breaks, t_span[1]]
        t_parts, y_parts, y = [], [], y0
        for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
            ctx['segment'] = k
            last = k == len(edges) - 2
            te = t_eval[(t_eval >= a) & ((t_eval <= b) if last else (t_eval < b))]
            sol = solve_ivp(fun, (a, b), y, t_eval=te if last else np.append(te, b), **solver_kw)
            if not sol.success:
                raise RuntimeError(f"integration failed ({model}, segment {k}): {sol.message}")
            y = sol.y[:, -1]
            t_parts.append(sol.t[:len(te)])
            y_parts.append(sol.y[:, :len(te)])
        t_out, y_out = np.concatenate(t_parts), np.concatenate(y_parts, axis=1)

    res = dict(t=t_out, DG=ctx['DG'], weights=ctx['weights'],
               params=p, mode=mode, M=ctx['M'], model=model)
    res.update(layout.unpack(y_out))
    if drive is not None:
        res['N_A'] = N_A_of(t_out)
    if spec.postprocess is not None:
        spec.postprocess(res, p, ctx)
    return res


# ============================================================
# Unified response models (driven by run_meanfield.py)
# ============================================================
# A response is a set of B-cell populations on the same DG grid, 'naive'
# (state suffix '') and/or 'memory' (suffix '_m'), sharing antigen and T cells.
# Primary = naive only; recall = memory, with or without the naive repertoire.
# Each population has its own pMHC load, activation latch and h0
# (naive: p.h0, memory: p.h0_mem, default p.h0). Memory bins with no cells are
# frozen (integrating an exactly-empty population stalls the solver).
#
#   model   : 'continuous' -- loading, T-cell help and division at the same time
#                             (the semicomplete equations; T cells consumed/recycled)
#             'two_stage'  -- loading for t < t_switch, then expansion without
#                             uptake, limited by a finite pool of N_T0 T cells
#   antigen : 'innate'   -- antigen cut by the innate response (autonomous)
#             'feedback' -- antigen cut by the response's potency Z
#             'replay'   -- antigen replayed from a given trajectory (drive)
POP_SUFFIX = {'naive': '', 'memory': '_m'}
ANTIGEN_LAWS = ('innate', 'feedback', 'replay')


def _antigen_rate(t, N_A, p, Z, antigen):
    """dN_A/dt for the innate or feedback law (same arithmetic as _dN_A_with_feedback)."""
    if antigen == 'feedback':
        frac = (1.0 + p.Z_c / (max(Z, 0.0) + 1e-30)) ** (-1)
    else:
        innate_proxy = np.exp(p.lambda_innate * t)
        frac = (1.0 + p.threshold_innate / (innate_proxy + 1e-30)) ** (-1)
    return (p.lambda_A * (1 - frac) - p.delta_A * frac) * N_A - 0.01 * N_A


def _rhs_continuous(t, s, p, ctx):
    psi_vec, weights, neut = ctx['psi'], ctx['weights'], ctx['neut']
    N_A = max(s['N_A'], 0.0)
    N_To = max(s['N_To'], 0.0)
    N_Ta = max(s['N_Ta'], 0.0)
    T_term = N_To / (N_To + p.K_T)

    d, flux_T, Z = {}, 0.0, 0.0
    for sfx, h0, live in ctx['pops']:
        pi_vec = np.maximum(s['pi' + sfx], 0.0)
        N_Bo = np.maximum(s['N_Bo' + sfx], 0.0)
        N_Ba = np.maximum(s['N_Ba' + sfx], 0.0)
        a_vec = s['a' + sfx]

        pi_star = (p.b0 / h0) ** (1 / p.hill)
        g = pi_vec ** p.hill / (pi_vec ** p.hill + pi_star ** p.hill)
        N_B_tot = N_Bo + N_Ba
        with np.errstate(divide='ignore', invalid='ignore'):
            lambda_eff = np.where(N_B_tot > 0, p.b0 * N_Ba / N_B_tot, 0.0)
        act = h0 * pi_vec ** p.hill * N_Bo * N_To / (N_To + p.K_T)

        dd = {'a': (10.0 * p.b0) * (1.0 - a_vec) * g,
              'pi': p.k_on * psi_vec * N_A - p.delta_pi * pi_vec - lambda_eff * pi_vec,
              'N_Bo': - act + 2 * p.b0 * N_Ba - p.delta_B * N_Bo,
              'N_Ba': act - p.b0 * N_Ba - p.delta_B * N_Ba}
        if live is not None:
            dd = {k: v * live for k, v in dd.items()}
            flux_T += h0 * np.sum(weights * pi_vec ** p.hill * N_Bo * live) * T_term
        else:
            flux_T += h0 * np.sum(weights * pi_vec ** p.hill * N_Bo) * T_term
        d.update({k + sfx: v for k, v in dd.items()})
        Z += np.sum(neut * N_B_tot * a_vec)

    d['N_A'] = _antigen_rate(t, N_A, p, Z, ctx['antigen'])
    d['N_To'] = - flux_T + p.Tcell_growth_factor * p.b0 * N_Ta - p.delta_T * N_To
    d['N_Ta'] = flux_T - p.b0 * N_Ta - p.delta_T * N_Ta
    return d


def _rhs_two_stage_pops(t, s, p, ctx):
    psi_vec, weights, neut = ctx['psi'], ctx['weights'], ctx['neut']
    N_A = max(s['N_A'], 0.0)
    loading = ctx['segment'] == 0 if 'segment' in ctx else t < p.t_switch
    engaged = sum(np.sum(weights * np.maximum(s['N_BT' + sfx], 0.0)) for sfx, _, _ in ctx['pops'])
    T_free = max(p.N_T0 - engaged, 0.0) / p.N_T0            # fraction of T cells free

    d, Z = {}, 0.0
    for sfx, h0, live in ctx['pops']:
        pi_vec = np.maximum(s['pi' + sfx], 0.0)
        N_Bo = np.maximum(s['N_Bo' + sfx], 0.0)
        N_BT = np.maximum(s['N_BT' + sfx], 0.0)
        N_Ba = np.maximum(s['N_Ba' + sfx], 0.0)
        a_vec = s['a' + sfx]

        N_B_tot = N_Bo + N_BT + N_Ba
        pi_star = (p.b0 / h0) ** (1 / p.hill)
        g = pi_vec ** p.hill / (pi_vec ** p.hill + pi_star ** p.hill)
        if loading:
            dpi = p.k_on * psi_vec * N_A - p.delta_pi * pi_vec
            eng = rel = div = 0.0
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                lambda_eff = np.where(N_B_tot > 0, p.b0 * N_Ba / N_B_tot, 0.0)
            dpi = - p.delta_pi * pi_vec - lambda_eff * pi_vec   # dilution by division
            h = h0 * pi_vec ** p.hill
            eng = h / (1.0 + h / p.h_max) * N_Bo * T_free      # engagement capped at h_max
            rel = N_BT / p.tau_eng
            div = p.b0 * N_Ba

        dd = {'pi': dpi,
              'N_Bo': - eng + 2 * div - p.delta_B * N_Bo,
              'N_BT': eng - rel - p.delta_B * N_BT,
              'N_Ba': rel - div - p.delta_B * N_Ba,
              'a': (10.0 * p.b0) * (1.0 - a_vec) * g}
        if live is not None:
            dd = {k: v * live for k, v in dd.items()}
        d.update({k + sfx: v for k, v in dd.items()})
        Z += np.sum(neut * N_B_tot * a_vec)

    d['N_A'] = _antigen_rate(t, N_A, p, Z, ctx['antigen'])
    return d


# atol far below per-clone memory sizes, which can be << 1 cell (sparse memory pools)
RESPONSE_MODELS = {
    'continuous': dict(pop_blocks=('pi', 'N_Bo', 'N_Ba', 'a'),
                       shared=(('N_To', 'scalar'), ('N_Ta', 'scalar')),
                       shared_init=lambda p: dict(N_To=p.N_T0, N_Ta=0.0),
                       rhs=_rhs_continuous, solver=dict(method='LSODA', rtol=1e-4, atol=1e-10),
                       breaks=None),
    'two_stage':  dict(pop_blocks=('pi', 'N_Bo', 'N_BT', 'N_Ba', 'a'),
                       shared=(), shared_init=lambda p: {},
                       rhs=_rhs_two_stage_pops, solver=TWO_STAGE_SOLVER,
                       breaks=lambda p: [p.t_switch]),
}


def _post_response(res, p, ctx):
    """Per-population totals N_B<sfx> (all B-cell blocks) and the free T-cell fraction."""
    vec_blocks = [b for b in RESPONSE_MODELS[res['model']]['pop_blocks'] if b.startswith('N_B')]
    for name in res['pops']:
        sfx = POP_SUFFIX[name]
        res['N_B' + sfx] = sum(res[b + sfx] for b in vec_blocks)
    if res['model'] == 'two_stage':
        engaged = sum(np.sum(ctx['weights'][:, None] * res['N_BT' + POP_SUFFIX[n]], axis=0)
                      for n in res['pops'])
        res['T_free'] = np.maximum(p.N_T0 - engaged, 0.0) / p.N_T0
    else:
        res['T_free'] = res['N_To'] / (res['N_To'] + p.K_T)


def run_response(model, p, naive=True, memory_seed=None, antigen='innate', drive=None,
                 t_span=(0.0, 10.0), t_eval=None, mode='grid', DG_max_sim=None, seed=None,
                 **solver_overrides):
    """
    Run one response of a unified model.

    model       : 'continuous' | 'two_stage'
    naive       : include the naive repertoire (initial size ctx['N0'] per clone)
    memory_seed : per-clone memory sizes on the grid (None -> no memory population)
    antigen     : 'innate' | 'feedback' | 'replay' ('replay' needs drive=dict(t, N_A))
    Returns the usual results dict plus 'pops' (population names), 'antigen',
    N_B / N_B_m (per-clone totals) and T_free.
    """
    if model not in RESPONSE_MODELS:
        raise ValueError(f"Unknown model '{model}'. Choose from {list(RESPONSE_MODELS)}.")
    if antigen not in ANTIGEN_LAWS:
        raise ValueError(f"Unknown antigen law '{antigen}'. Choose from {ANTIGEN_LAWS}.")
    if (antigen == 'replay') != (drive is not None):
        raise ValueError("antigen='replay' needs drive=dict(t=..., N_A=...), and only it.")
    rm = RESPONSE_MODELS[model]
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500)

    ctx = build_context(p, mode=mode, DG_max_sim=DG_max_sim, seed=seed)
    ctx['antigen'] = antigen
    pops, init = [], dict(N_A=p.N_A0)
    if naive:
        pops.append(('naive', ctx['N0'].copy(), p.h0, None))
    if memory_seed is not None:
        seed_m = np.asarray(memory_seed, dtype=float).copy()
        pops.append(('memory', seed_m, p.h0 if p.h0_mem is None else p.h0_mem, seed_m > 0))
    if not pops:
        raise ValueError("A response needs naive=True and/or a memory_seed.")
    ctx['pops'] = [(POP_SUFFIX[name], h0, live) for name, _, h0, live in pops]
    blocks = [('N_A', 'scalar')]
    for name, N_Bo, _, _ in pops:
        sfx = POP_SUFFIX[name]
        blocks += [(b + sfx, 'vector') for b in rm['pop_blocks']]
        init.update({b + sfx: (N_Bo if b == 'N_Bo' else np.zeros(ctx['M'])) for b in rm['pop_blocks']})
    blocks += list(rm['shared'])
    init.update(rm['shared_init'](p))

    spec = ModelSpec(name=model, blocks=blocks, rhs=rm['rhs'], init=None,
                     solver=rm['solver'], breaks=rm['breaks'])
    res = _integrate(spec, p, ctx, init, t_span, t_eval, drive, mode, model, solver_overrides,
                     strict=True)
    res.update(pops=[name for name, _, _, _ in pops], antigen=antigen)
    _post_response(res, p, ctx)
    return res


def pop_view(res, pop):
    """
    One population of a run_response result as a single-population results dict
    (unsuffixed keys N_Bo, N_Ba, pi, a, N_B; params with h0 = that population's h0),
    so the analysis functions below (potency, specificity, memory seeding, ...)
    apply to it unchanged. pop='all' sums the populations' cells (pi, a of naive).
    """
    p = res['params']
    if pop == 'all':
        view = pop_view(res, res['pops'][0])
        for other in res['pops'][1:]:
            view['N_B'] = view['N_B'] + res['N_B' + POP_SUFFIX[other]]
        return view
    sfx = POP_SUFFIX[pop]
    h0 = p.h0 if (pop == 'naive' or p.h0_mem is None) else p.h0_mem
    view = {k: res[k] for k in ('t', 'DG', 'weights', 'mode', 'M', 'N_A', 'T_free')}
    view.update(params=replace(p, h0=h0), model=res['model'], pop=pop,
                **{b: res[b + sfx] for b in RESPONSE_MODELS[res['model']]['pop_blocks']},
                N_B=res['N_B' + sfx])
    return view


# ============================================================
# Analysis
# ============================================================
def _N_B(res):
    """Total per-clone B cells, handling models with/without N_Ba.
    Models that store their own total ('approx', 'two_stage': N_Bo + N_BT + N_Ba) use it."""
    if 'N_B' in res:
        return res['N_B']
    return res['N_Bo'] + res['N_Ba']


def compute_N_B_tot(res, threshold=2.0):
    """Total number of activated B cells over time."""
    N_B = _N_B(res)
    thr = N_B[:, 0] + 1
    activated = (N_B > thr[:, None]).astype(float)
    return np.sum(res['weights'][:, None] * N_B * activated, axis=0)


def compute_L_act(res, threshold=2.0):
    """Number of activated clones over time."""
    N_B = _N_B(res)
    activated = (N_B > threshold).astype(float)
    return np.sum(res['weights'][:, None] * activated, axis=0)


def compute_N1(res):
    return np.max(_N_B(res), axis=0)


def compute_potency_t(res, threshold=2.0):
    """
    Omega-weighted, neutralization-weighted potency Z(t).

    Dispatches on the model: the 'null' model uses the integrated activation
    latch `a`; models with an explicit pMHC use the pi > pi_star criterion.
    """
    N_tot = _N_B(res)
    neut = (res['weights'] * np.exp(-res['DG']))[:, None]
    if res.get('model') == 'null':
        return np.sum(neut * N_tot * res['a'], axis=0)
    p = res['params']
    pi_star = (p.b0 / p.h0) ** (1 / p.hill)
    ever_active = np.maximum.accumulate(res['pi'] > pi_star, axis=1).astype(float)
    return np.sum(neut * N_tot * ever_active, axis=0)


# Kept for backward compatibility with existing calls.
def compute_potency_null_t(res, threshold=2.0):
    N_tot = _N_B(res)
    neut = (res['weights'] * np.exp(-res['DG']))[:, None]
    return np.sum(neut * N_tot * res['a'], axis=0)


def compute_potency(res, time_index=-1, threshold=2.0):
    """Potency at a single time index (sum of clone size * exp(-DG))."""
    N_B = _N_B(res)[:, time_index]
    DG, w = res['DG'], res['weights']
    activated = N_B > threshold
    return np.sum(w[activated] * N_B[activated] * np.exp(-DG[activated]))


def compute_yield(res, time_index=-1, threshold=2.0):
    """Yield at a single time index (total activated B cells)."""
    N_B = _N_B(res)[:, time_index]
    w = res['weights']
    activated = N_B > threshold
    return np.sum(w[activated] * N_B[activated])


def find_t_c(res, Z_c=None):
    """Index/time when potency Z(t) first crosses Z_c. idx=None if never."""
    Z = compute_potency_t(res)
    if Z_c is None:
        Z_c = res['params'].Z_c
    hit = np.where(Z >= Z_c)[0]
    if len(hit) == 0:
        return None, np.inf
    return hit[0], res['t'][hit[0]]


def find_t_star(res, Z_c=None):
    """Index/time when potency Z(t) first crosses Z_c. idx=None if never."""
    Z = compute_potency_t(res)
    if Z_c is None:
        Z_c = res['params'].Z_c
    hit = np.where(Z >= 0)[0]
    if len(hit) == 0:
        return None, np.inf
    return hit[0], res['t'][hit[0]]


def compute_potency_P(res, t_c_star=None):
    """Potency P = -log[N_A(t_c)]/lambda_A (natural log). nan if Z_c unreached."""
    # Z = compute_potency_t(res)
    Z_c = res['params'].Z_c
    # P = Z[res['t']>t_c_star][0] if t_c_star is not None else np.inf
    idc, t_c = find_t_c(res, Z_c)
    # idstar, t_star = find_t_star(res, Z_c)
    # if Z_c is None:
    #     Z_c = res['params'].Z_c
    P = res['params'].lambda_A*(4.4-t_c) if idc is not None else np.nan
    return t_c_star, P


def compute_specificity(res, threshold=2.0):
    r"""
    D_KL(Omega || Omega_0): responding density vs naive density (nats).

    Omega(DG) = weights * final N_B, restricted to the ACTIVATED population using
    the SAME indicator as the potency:
        null model : the integrated activation latch a(t_end)
        pMHC models: clones that ever crossed pi > pi_star
    This keeps specificity and potency consistent and removes the artefact where a
    weakly-proliferating null response (few clones above a fixed N_B cutoff) looks
    spuriously narrow/specific. `threshold` is kept for signature compatibility.
    """
    N_B = _N_B(res)[:, -1]
    w = res['weights']
    if res.get('model') == 'null':
        active = np.clip(res['a'][:, -1], 0.0, 1.0)
    else:
        p = res['params']
        pi_star = (p.b0 / p.h0) ** (1 / p.hill)
        active = np.maximum.accumulate(res['pi'] > pi_star, axis=1)[:, -1].astype(float)
    Omega = w * N_B * active
    total = Omega.sum()
    if total <= 0:
        return np.nan
    p_ = Omega / total
    q = w / w.sum()
    mask = p_ > 0
    return float(np.sum(p_[mask] * np.log(p_[mask] / q[mask])))


def dkl_analytic(mu, DG_mf, DG_star=-16.0, sigma2=7.0):
    """
    Analytic D_KL of a double-truncated, tilted Gaussian responding distribution
    vs a naive N(0, sigma2). Truncation to [DG_star, DG_mf]:
        DG_star : high-affinity edge (negative here, = -beta* sigma2)
        DG_mf   : front position (absolute) = DG_star + front_DG_mf(res)
        mu      : tilt of the responding Gaussian
    D = delta^2/2 + delta*m1 - log Z,  delta=mu/sigma, m1=g/Z,
    g = phi(z') - phi(z''),  Z = Phi(z'') - Phi(z'),
    z' = (DG_star - mu)/sigma,  z'' = (DG_mf - mu)/sigma.
    """
    sigma = np.sqrt(sigma2)
    zp  = (DG_star - mu) / sigma
    zpp = (DG_mf   - mu) / sigma
    Z = norm.cdf(zpp) - norm.cdf(zp)
    if Z <= 0:
        return np.nan
    g  = norm.pdf(zp) - norm.pdf(zpp)
    m1 = g / Z
    delta = mu / sigma
    return delta ** 2 / 2 + delta * m1 - np.log(Z)


def entropy_analytic(mu, DG_mf, DG_star=-16.0, sigma2=7.0):
    """
    Differential entropy (nats) of the same double-truncated, tilted Gaussian as
    dkl_analytic: N(mu, sigma2) restricted to [DG_star, DG_mf].
    H = log(sqrt(2 pi e) sigma Z) + (z' phi(z') - z'' phi(z'')) / (2 Z).
    """
    sigma = np.sqrt(sigma2)
    zp  = (DG_star - mu) / sigma
    zpp = (DG_mf   - mu) / sigma
    Z = norm.cdf(zpp) - norm.cdf(zp)
    if Z <= 0:
        return np.nan
    return (np.log(np.sqrt(2 * np.pi * np.e) * sigma * Z)
            + (zp * norm.pdf(zp) - zpp * norm.pdf(zpp)) / (2 * Z))


def find_t_D(res, D_threshold=1.0):
    """Time at which the demand D(t) first crosses D_threshold (approx model)."""
    idx = np.where(res['D'] >= D_threshold)[0]
    return res['t'][idx[0]] if len(idx) else np.inf


def compute_zipf(res, time_index=-1, threshold=2.0):
    """Zipf (rank-size) distribution at a given time, normalized by the top clone."""
    N_B = _N_B(res)[:, time_index]
    w = res['weights']
    if res['mode'] == 'stochastic':
        sizes = np.sort(N_B[N_B > threshold])[::-1]
    else:
        expanded = []
        for i in range(len(N_B)):
            if N_B[i] > threshold:
                expanded.extend([N_B[i]] * max(1, int(np.round(w[i]))))
        sizes = np.array(sorted(expanded, reverse=True))
    ranks = np.arange(1, len(sizes) + 1)
    return ranks, sizes / sizes[0]


def DG_mf(res, threshold=2.0):
    """Front position from the sim: the largest activated DG (max DG whose final
    N_B exceeds threshold), in the sim frame == Fig.4C's (DG - DG*)."""
    NB = _N_B(res)[:, -1]
    act = NB > threshold
    return res['DG'][act].max() if act.any() else np.nan


def DG_min(res, threshold=2.0):
    """Minimum DG from the sim: the smallest activated DG (min DG whose final
    N_B exceeds threshold), in the sim frame == Fig.4C's (DG - DG*)."""
    NB = _N_B(res)[:, -1]
    act = NB > threshold
    return res['DG'][act].min() if act.any() else np.nan


# ============================================================
# Memory formation
# ============================================================
def produce_memory(res, n_mem=int(1e4)):
    """Sample a memory pool (DG, counts) from the final activated population."""
    per_clone = _N_B(res)[:, -1]
    w = res['weights']
    activated_idx = np.where(per_clone > 2.0)[0]
    counts = (per_clone * w)[activated_idx]
    total = counts.sum()
    if total <= 0:
        return np.array([]), np.array([])
    n_draw = (n_mem * total) / (n_mem + total)     # harmonic-mean draw size
    probs = counts / total
    drawn = np.random.multinomial(int(n_draw), probs)
    nz = drawn > 0
    return res['DG'][activated_idx[nz]], drawn[nz]


def memory_seed_from_primary(res_primary, p, n_mem=int(1e4), threshold=2.0):
    """
    Subsample the primary's activated population into a length-M grid IC.

    threshold : final per-clone size above which a clone counts as activated
                (scalar or per-clone array, e.g. 2 * initial size).

    Returns (expanded_cells, formed_memory) where formed_memory is the per-clone
    initial condition to pass to `run_simulation(..., memory_seed=formed_memory)`.
    """
    per_clone = _N_B(res_primary)[:, -1]
    w = res_primary['weights']
    activated = per_clone > threshold
    expanded_cells = (per_clone * w) * activated
    total = expanded_cells.sum()
    if total <= 0:
        raise ValueError("Primary produced no activated cells to sample from.")

    n_draw = (n_mem * total) / (n_mem + total)
    probs = expanded_cells / total
    drawn = np.random.multinomial(int(n_draw), probs)

    formed_memory = np.zeros_like(per_clone)
    nz = w > 0
    formed_memory[nz] = drawn[nz] / w[nz]
    return expanded_cells, formed_memory


def memory_seed_expected(res_primary, n_mem=int(1e4), threshold=2.0):
    """
    Deterministic version of memory_seed_from_primary: the expected multinomial
    draw instead of a random one. Returns (expanded_cells, formed_memory) in the
    same units (cells per DG bin; per-clone memory size).
    """
    per_clone = _N_B(res_primary)[:, -1]
    w = res_primary['weights']
    expanded_cells = (per_clone * w) * (per_clone > threshold)
    total = expanded_cells.sum()
    if total <= 0:
        raise ValueError("Primary produced no activated cells to sample from.")
    n_draw = (n_mem * total) / (n_mem + total)
    formed_memory = np.zeros_like(per_clone)
    nz = w > 0
    formed_memory[nz] = n_draw * expanded_cells[nz] / total / w[nz]
    return expanded_cells, formed_memory
