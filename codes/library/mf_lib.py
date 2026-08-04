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

from dataclasses import dataclass, field
from typing import Callable, List, Tuple, Optional

import numpy as np
from scipy.integrate import solve_ivp

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


def build_context(p, mode='grid', DG_max_sim=None, seed=None):
    """Everything the RHS/IC need that only depends on the repertoire + params."""
    if mode == 'grid':
        DG, psi_arr, weights, M = build_repertoire_grid(p)
    elif mode == 'stochastic':
        DG, psi_arr, weights, M = build_repertoire_stochastic(p, DG_max_sim, seed)
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'grid' or 'stochastic'.")
    return dict(DG=DG, psi=psi_arr, weights=weights,
                neut=weights * np.exp(-DG), M=M, mode=mode)


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
    return np.ones(M)


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
}


# ============================================================
# Generic runner
# ============================================================
def run_simulation(model='semicomplete', p=None, t_span=None, t_eval=None,
                   mode='grid', DG_max_sim=None, seed=None, memory_seed=None,
                   **solver_overrides):
    """
    Run any registered model.

    Parameters
    ----------
    model : 'approx' | 'complete' | 'semicomplete' | 'null'
    p : Parameters
    t_span : (t0, t1)
    t_eval : output times (default 500 points over t_span)
    mode : 'grid' | 'stochastic'
    memory_seed : optional per-clone N_Bo initial condition (memory phase)
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
    layout = spec.layout(ctx['M'])

    y0 = layout.pack(spec.init(p, ctx, memory_seed))
    solver_kw = dict(spec.solver)
    solver_kw.update(solver_overrides)

    def fun(t, y):
        return layout.pack(spec.rhs(t, layout.unpack(y), p, ctx))

    sol = solve_ivp(fun, t_span, y0, t_eval=t_eval, **solver_kw)
    if not sol.success:
        print(f"Warning: integration failed ({model}): {sol.message}")

    res = dict(t=sol.t, DG=ctx['DG'], weights=ctx['weights'],
               params=p, mode=mode, M=ctx['M'], model=model)
    res.update(layout.unpack(sol.y))
    if spec.postprocess is not None:
        spec.postprocess(res, p, ctx)
    return res


# ============================================================
# Analysis
# ============================================================
def _N_B(res):
    """Total per-clone B cells, handling models with/without N_Ba."""
    if 'N_Ba' in res:
        return res['N_Bo'] + res['N_Ba']
    return res['N_B']


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


def memory_seed_from_primary(res_primary, p, n_mem=int(1e4)):
    """
    Subsample the primary's activated population into a length-M grid IC.

    Returns (expanded_cells, formed_memory) where formed_memory is the per-clone
    initial condition to pass to `run_simulation(..., memory_seed=formed_memory)`.
    """
    per_clone = _N_B(res_primary)[:, -1]
    w = res_primary['weights']
    activated = per_clone > 2.0
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
