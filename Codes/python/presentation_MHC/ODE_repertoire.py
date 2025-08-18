import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---------- 1) Binning over log k with a user-supplied log-density Ω0(k) ----------
def make_log_bins(omega0_logdens: Callable[[np.ndarray], np.ndarray],
                  kmin: float, kmax: float, B: int, log_base: float = np.e,
                  N_total: int = 10_000) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build B bins over log(k) in [kmin, kmax], return:
      k_b (centers), w_b (weights summing to 1 over bins), N_b (integer counts summing to N_total)
    omega0_logdens is a *log-density*: density per d(ln k) if log_base=e (default).
    If your Ω0 is per d(log10 k), set log_base=10 (the Δlog factor handles it).
    """
    if not (kmin > 0 and kmax > kmin):
        raise ValueError("Require 0 < kmin < kmax.")
    # Edges equally spaced in the chosen log
    log = np.log if log_base == np.e else (lambda x: np.log10(x))
    dlog = (log(kmax) - log(kmin)) / B
    edges_log = log(kmin) + dlog * np.arange(B + 1)
    # Geometric centers
    centers = np.exp(edges_log[:-1] + 0.5 * dlog) if log_base == np.e else 10 ** (edges_log[:-1] + 0.5 * dlog)

    # Midpoint rule for ∫ Ω0(k) d(log k) over each bin
    dens_mid = omega0_logdens(centers)
    dens_mid = np.clip(dens_mid, 0.0, np.inf)  # guard
    # Each bin mass ∝ Ω0(k_mid) * Δlog
    mass = dens_mid * dlog
    if np.all(mass == 0):
        raise ValueError("Omega0(log-density) is zero everywhere on [kmin, kmax].")
    w = mass / mass.sum()

    # Integer counts per bin that sum to N_total
    N_float = N_total * w
    N_int = np.floor(N_float).astype(int)
    # Distribute the remainder by largest fractional parts (“largest remainder” method)
    remainder = N_total - N_int.sum()
    if remainder > 0:
        frac = N_float - N_int
        take = np.argsort(frac)[::-1][:remainder]
        N_int[take] += 1

    return centers, w, N_int

# ---------- 2) Parameters ----------
@dataclass
class ModelParams:
    alpha: int                    # number of maturation steps (m_1..m_alpha)
    k_step: float                 # forward step rate
    k_m_plus: float               # motor on-rate (per motor)
    k_er_plus: float = 0.0        # optional source from A0; set 0 for closed system
    k_er_minus: float = 0.0       # ER loss
    k0_minus: float = 0.0         # final pool loss (decay/clearance)
    M_tot: float = 1_000.0        # total motors
    include_A_as_bound: bool = True  # count A as motor-bound occupancy

# ---------- 3) Cohort ODE (vectorized across bins) ----------
def build_rhs(k_bins: np.ndarray, N_per_bin: np.ndarray, params: ModelParams, A0_source: np.ndarray = None):
    """
    Returns an RHS function f(t, y) for solve_ivp.
    State layout y = [ER(0..B-1), m1(0..B-1), ..., m_alpha(0..B-1), A(0..B-1)]
    All variables represent *counts* of peptides (not concentrations).
    """
    B = len(k_bins)
    alpha = params.alpha
    k_step = params.k_step
    k_m_plus = params.k_m_plus
    k_er_plus = params.k_er_plus
    k_er_minus = params.k_er_minus
    k0_minus = params.k0_minus
    include_A = params.include_A_as_bound

    if A0_source is None:
        A0_source = np.ones(B)
    A0_source = np.asarray(A0_source, dtype=float)

    # indices/slices
    def unpack(y):
        y = y.reshape(-1)
        ER = y[0:B]
        m = np.empty((alpha, B), dtype=float)
        for i in range(alpha):
            m[i] = y[B*(1+i):B*(2+i)]
        A = y[B*(1+alpha):B*(2+alpha)]
        return ER, m, A

    def rhs(t, y):
        ER, m, A = unpack(y)
        # Bound motors = sum over all m-steps (+ A if counted)
        bound = m.sum() + (A.sum() if include_A else 0.0)
        M_free = params.M_tot - bound
        if M_free < 0:  # numerical guard; prevents negative free motors
            M_free = 0.0

        dER = A0_source * k_er_plus - (k_er_minus + k_m_plus * M_free) * ER + (k_bins * m.sum(axis=0))
        dm = np.empty_like(m)
        # m1
        dm[0] = k_m_plus * M_free * ER - (k_step + k_bins) * m[0]
        # m2..m_alpha
        for i in range(1, alpha):
            dm[i] = k_step * m[i-1] - (k_step + k_bins) * m[i]
        dA = k_step * m[-1] - (k0_minus + k_bins) * A  # motor unbinds from A at rate k_bins

        # Pack
        out = np.concatenate([dER] + [dm[i] for i in range(alpha)] + [dA])
        return out

    return rhs

# ---------- 4) Simulation helper ----------
def simulate_cohorts(omega0_logdens: Callable[[np.ndarray], np.ndarray],
                     kmin: float, kmax: float, B: int, N_total: int,
                     params: ModelParams,
                     t_span=(0.0, 200.0), t_eval=None,
                     init_mode: str = "all_in_ER",
                     A0_source_scale: float = 0.0,
                     log_base: float = np.e):
    """
    Build bins, set initial conditions, integrate ODEs, and return (sol, k_bins, N_per_bin).
    init_mode:
      - "all_in_ER": distribute the N_total peptides into ER per bin (no external source).
      - "external_source": start empty and drive ER with A0_source = N_per_bin * A0_source_scale (per-time inflow).
    A0_source_scale has units of 1/time (fraction of the bin's cohort injected per unit time).
    """
    k_bins, w_bins, N_per_bin = make_log_bins(omega0_logdens, kmin, kmax, B, log_base, N_total)
    B = len(k_bins)
    alpha = params.alpha

    # Initial conditions
    if init_mode == "all_in_ER":
        ER0 = N_per_bin.astype(float)
        A0_src = np.zeros(B)
    elif init_mode == "external_source":
        ER0 = np.zeros(B)
        # A0_source proportional to the cohort size (peptides per time)
        A0_src = N_per_bin.astype(float) * float(A0_source_scale)
    else:
        raise ValueError("init_mode must be 'all_in_ER' or 'external_source'.")

    m0 = np.zeros((alpha, B), dtype=float)
    A0 = np.zeros(B, dtype=float)

    y0 = np.concatenate([ER0] + [m0[i] for i in range(alpha)] + [A0])
    rhs = build_rhs(k_bins, N_per_bin, params, A0_src)

    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 501)

    # Choose a stiff solver if rates span orders of magnitude
    sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, method="BDF", atol=1e-8, rtol=1e-6, vectorized=False)
    return sol, k_bins, N_per_bin

# ---------- 5) Example Ω0 (log-normal log-density) ----------
def omega0_lognormal_lndensity(k, mu=np.log(0.1), sigma=1.0):
    """
    Returns log-density Ω0(k) = density per d(ln k) for a log-normal LN(mu, sigma).
    For LN, pdf over k is: f(k)= (1/(k sigma sqrt(2π))) exp(-(ln k - mu)^2/(2 sigma^2))
    Then Ω0(k) = f(k) * k = (1/(sigma sqrt(2π))) * exp(-(ln k - mu)^2/(2 sigma^2))
    """
    k = np.asarray(k, dtype=float)
    z = (np.log(k) - mu) / sigma
    return np.exp(-0.5 * z * z) / (sigma * np.sqrt(2 * np.pi))

# ---------- 6) Quick demo run ----------
if __name__ == "__main__":
    params = ModelParams(
        alpha=2,         # matches your earlier derivation (power 2)
        k_step=0.2,
        k_m_plus=1e-4,
        k_er_plus=0.0,   # closed system by default
        k_er_minus=0.0,
        k0_minus=0.01,
        M_tot=2_000,     # try also smaller values to see depletion
        include_A_as_bound=True
    )

    sol, k_bins, N_per_bin = simulate_cohorts(
        omega0_logdens=omega0_lognormal_lndensity,
        kmin=1e-4, kmax=10.0, B=150, N_total=10_000,
        params=params,
        t_span=(0.0, 500.0),
        init_mode="all_in_ER"  # or "external_source" with k_er_plus>0 and A0_source_scale>0
    )

    # Derived outputs
    B = len(k_bins); alpha = params.alpha
    def unpack(y):
        ER = y[0:B]
        m = np.vstack([y[B*(1+i):B*(2+i)] for i in range(alpha)])
        A = y[B*(1+alpha):B*(2+alpha)]
        return ER, m, A

    ER, m, A = unpack(sol.y[:, -1])  # steady state snapshot
    print(m)
    plt.plot(k_bins, omega0_lognormal_lndensity(k_bins), ls = '--')
    plt.plot(k_bins, m.T)
    plt.plot(k_bins, ER.T)
    plt.plot(k_bins,A.T)
    plt.yscale("log")
    plt.xscale("log")
    plt.show()
    bound = m.sum(axis=0).sum() + (A.sum() if params.include_A_as_bound else 0.0)
    P_free_ss = (params.M_tot - bound) / params.M_tot
    print(f"Steady-state P_free ≈ {P_free_ss:.4f}  (M_free={params.M_tot - bound:.1f} of M_tot={params.M_tot})")

    total_peptides = ER.sum() + m.sum() + A.sum()
    print(f"Total peptides (closed system) ≈ {total_peptides:.1f} (should stay near 10,000 if k_er_minus=k0_minus=0).")
