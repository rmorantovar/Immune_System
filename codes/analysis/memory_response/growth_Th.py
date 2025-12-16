import sys
sys.path.append('../../my_lib/')
from funcs import*
plt.rcParams['text.usetex'] = True
import math
import heapq
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple, Optional

text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'
project = 'memory_response'
subproject = 'growth_Th'

output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
os.makedirs(output_plot, exist_ok=True)

# ==========================
# Model utilities
# ==========================

def rho_A(t: float, rhoA0: float, lamA: float) -> float:
    """
    External antigen concentration/abundance.

    We want rho_A(0) = 0 and exponential growth afterward:
        rho_A(t) = rhoA0 * (exp(lamA*t) - 1)

    Parameters
    ----------
    t : float
        Time.
    rhoA0 : float
        Scale of antigen growth.
    lamA : float
        Exponential growth rate.

    Returns
    -------
    float
        Antigen concentration at time t.
    """
    return rhoA0 * (np.exp(lamA * t) - 1.0)


def hill_saturation(pi: float, K: float, n: float = 1.0) -> float:
    """
    Saturating dependence of help success on presented antigen.
    H(pi) = pi^n / (K^n + pi^n)

    Notes
    -----
    - For n=1, this is Michaelis-Menten.
    - In the small-pi regime pi << K, H(pi) ~ (pi/K)^n.

    Returns in [0, 1].
    """
    if pi <= 0.0:
        return 0.0
    pin = pi**n
    Kn = K**n
    return pin / (Kn + pin)


def p_evolve(t0: float,p0: float,t: float,phi: float,kp: float,deltapi: float,rhoA0: float,lamA: float,) -> float:
    """
    Deterministic evolution of presented antigen pi(t) for a single cell
    between events, under:

        dpi/dt = kp * phi * rho_A(t) - deltapi * pi
        rho_A(t) = rhoA0*(exp(lamA*t) - 1)

    This has an analytic closed form (piecewise for deltapi=0).

    Parameters
    ----------
    t0 : float
        Start time of the interval (event time).
    pi0 : float
        pi(t0), presented antigen at start time.
    t : float
        End time.
    phi : float
        Clone-specific affinity factor.
    kp : float
        Capture/presentation accumulation coefficient.
    deltapi : float
        Turnover/decay rate of presented antigen (>=0).
    rhoA0 : float
        Antigen scale.
    lamA : float
        Antigen exponential growth rate.

    Returns
    -------
    float
        pi(t).
    """
    if t <= t0:
        return p0

    dt = t - t0

    # Special case deltapi = 0: pi integrates kp*phi*rho_A(t)
    if deltapi == 0.0:
        # ∫_{t0}^{t} rhoA0*(e^{lamA*s}-1) ds
        # = rhoA0 * [ (e^{lamA*t}-e^{lamA*t0})/lamA - (t-t0) ]
        term = rhoA0 * ((np.exp(lamA * t) - np.exp(lamA * t0)) / lamA - dt)
        return p0 + kp * phi * term

    # General deltapi > 0: solution of linear ODE with forcing
    # pi(t) = p0 e^{-deltapi dt} + kp*phi ∫_{t0}^{t} e^{-deltapi (t-s)} rho_A(s) ds
    ed = np.exp(-deltapi * dt)

    # Compute the integral analytically:
    # ∫ e^{-deltapi (t-s)} rhoA0*(e^{lamA*s}-1) ds
    #
    # Let’s change variable u = s - t0, dt = t - t0, s = t0 + u
    # Then:
    #   ∫_{0}^{dt} e^{-deltapi (dt-u)} rhoA0*(e^{lamA*(t0+u)}-1) du
    # = rhoA0 * [ e^{lamA t0} * (e^{lamA dt} - e^{-deltapi dt})/(lamA + deltapi)  - (1 - e^{-deltapi dt})/deltapi ]
    #
    A = np.exp(lamA * t0) * (np.exp(lamA * dt) - np.exp(-deltapi * dt)) / (lamA + deltapi)
    B = (1.0 - np.exp(-deltapi * dt)) / deltapi
    integral = rhoA0 * (A - B)

    p_t = p0 * ed + kp * phi * integral
    return max(0.0, p_t)  # guard against tiny negative from floating error


# ==========================
# Event-driven simulation
# ==========================

@dataclass
class CellState:
    """
    Represents one B cell right after an event.

    We simulate only *division completion* events in the global event queue.
    Each cell that exists has a scheduled division completion time, which is:
        t_div = t_help + tau0

    For each cell we store:
    - clone id
    - last_event_time: time when its current p0 is defined (e.g. division time)
    - p0: presented antigen at that last_event_time
    - t_div: next division completion time (scheduled)
    """
    clone_id: int
    last_event_time: float
    p0: float
    t_div: float


def draw_help_time_by_thinning(
    rng: np.random.Generator,
    t0: float,
    p0: float,
    clone_phi: float,
    kp: float,
    deltapi: float,
    rhoA0: float,
    lamA: float,
    kh: float,
    T_const: float,
    K_help: float,
    n_help: float,
    g_func: Callable[[float, float, float], float] = hill_saturation,
    t_max: float = math.inf,
) -> float:
    """
    Draw the time of the next 'help authorization' event for ONE cell
    starting at time t0 with presentation p0.

    Help events occur with time-varying hazard:
        lambda(t) = kh * T_const * H(pi(t))
    where H(pi) in [0,1].

    We use thinning with constant upper bound:
        lambda_max = kh * T_const
    because H(pi) <= 1.

    Algorithm
    ---------
    Propose times from Exp(lambda_max), accept with prob H(pi(t)).
    This is exact for NHPP with rate lambda(t).

    Parameters
    ----------
    t_max : float
        Optional cutoff (if you want to stop drawing beyond t_end).

    Returns
    -------
    float
        Help time t_help (>= t0). If t_max is finite and exceeded, returns inf.
    """
    lam_max = kh * T_const
    if lam_max <= 0.0:
        return math.inf

    t = t0
    while True:
        # propose next point from homogeneous Poisson of rate lam_max
        t += rng.exponential(1.0 / lam_max)

        if t > t_max:
            return math.inf

        # compute presented antigen at time t
        p_t = p_evolve(t0, p0, t, clone_phi, kp, deltapi, rhoA0, lamA)

        # accept with probability H(p_t)
        accept_prob = g_func(p_t, K_help, n_help)
        if rng.random() < accept_prob:
            return t


def simulate_clonal_expansion(
    *,
    rng: np.random.Generator,
    # clones
    clone_phis: List[float],
    initial_cells_per_clone: List[int],
    # antigen
    rhoA0: float,
    lamA: float,
    # presentation dynamics
    kp: float,
    deltapi: float,
    # help dynamics
    kh: float,
    T_const: float,
    K_help: float,
    n_help: float = 1.0,
    # division program
    tau0: float = 6.0,
    uptake_during_division: bool = False,
    # simulation control
    t_end: float = 72.0,
    sample_dt: float = 0.5,
    max_cells: int = 200_000,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate multi-clone expansion with:
    - antigen starts at 0 and grows exponentially
    - pi accumulates deterministically between events
    - help events are NHPP, sampled exactly by thinning
    - upon help, division completes after tau0
    - at division, pi is diluted into daughters

    Returns
    -------
    times : array shape (num_samples,)
    counts : array shape (num_samples, num_clones)
        counts[k,i] = number of cells in clone i at times[k]
    """
    num_clones = len(clone_phis)
    assert len(initial_cells_per_clone) == num_clones

    # Sampling grid
    times = np.arange(0.0, t_end + 1e-12, sample_dt)
    counts = np.zeros((len(times), num_clones), dtype=int)

    # Global event queue of division completions:
    # heap items: (t_div, unique_id, CellState)
    # unique_id avoids comparing CellState objects in heap ties.
    event_heap: List[Tuple[float, int, CellState]] = []
    uid = 0

    # Track current cells per clone
    current_counts = np.array(initial_cells_per_clone, dtype=int)

    # Helper function to schedule a cell's next division
    def schedule_cell(clone_id: int, t_birth: float, p_birth: float) -> Optional[CellState]:
        nonlocal uid

        # Draw help time
        t_help = draw_help_time_by_thinning(
            rng=rng,
            t0=t_birth,
            p0=p_birth,
            clone_phi=clone_phis[clone_id],
            kp=kp,
            deltapi=deltapi,
            rhoA0=rhoA0,
            lamA=lamA,
            kh=kh,
            T_const=T_const,
            K_help=K_help,
            n_help=n_help,
            t_max=t_end,  # don't bother drawing beyond horizon
        )
        if not math.isfinite(t_help):
            return None

        # Determine pi at division completion time (depending on uptake during division)
        t_div = t_help + tau0

        return CellState(
            clone_id=clone_id,
            last_event_time=t_birth,
            p0=p_birth,
            t_div=t_div,
        )

    # Initialize starting population at t=0, with p0 = 0 (no antigen presented initially)
    for i in range(num_clones):
        for _ in range(initial_cells_per_clone[i]):
            st = schedule_cell(i, t_birth=0.0, p_birth=0.0)
            if st is not None:
                heapq.heappush(event_heap, (st.t_div, uid, st))
                uid += 1

    # Simulation main loop
    sample_idx = 0
    t_now = 0.0

    def record_up_to(t_target: float):
        """Record counts at all sampling times <= t_target."""
        nonlocal sample_idx
        while sample_idx < len(times) and times[sample_idx] <= t_target + 1e-12:
            counts[sample_idx, :] = current_counts
            sample_idx += 1

    record_up_to(0.0)

    # Run events in chronological order
    while event_heap:
        t_div, _, cell = heapq.heappop(event_heap)
        if t_div > t_end:
            break

        # Record samples up to this event time
        record_up_to(t_div)

        # Safety: stop if exploding too big
        total_cells = int(current_counts.sum())
        if total_cells >= max_cells:
            # Fill remainder and stop
            record_up_to(t_end)
            break

        clone_id = cell.clone_id

        # Compute pi at help time and at division time.
        # We do it by reconstructing the cell's trajectory from its last_event_time and p0.
        # First, we need t_help = t_div - tau0
        t_help = t_div - tau0

        # pi at help time (end of "waiting for help" stage)
        p_help = p_evolve(
            t0=cell.last_event_time,
            p0=cell.p0,
            t=t_help,
            phi=clone_phis[clone_id],
            kp=kp,
            deltapi=deltapi,
            rhoA0=rhoA0,
            lamA=lamA,
        )

        # Now evolve through division program interval [t_help, t_div]
        if uptake_during_division:
            # Continue same ODE during the division program
            p_div = p_evolve(
                t0=t_help,
                p0=p_help,
                t=t_div,
                phi=clone_phis[clone_id],
                kp=kp,
                deltapi=deltapi,
                rhoA0=rhoA0,
                lamA=lamA,
            )
        else:
            # No new uptake during division program; only turnover
            if deltapi == 0.0:
                p_div = p_help
            else:
                p_div = p_help * math.exp(-deltapi * tau0)

        # Division happens: mother is replaced by 2 daughters
        # Mother count stays same? In a branching view, mother "becomes" two daughters -> +1 net cell.
        # So clone count increases by +1 (from 1 to 2).
        current_counts[clone_id] += 1

        # Daughters inherit diluted presentation
        p_daughter = 0.5 * p_div

        # Schedule next division for each daughter
        # (Two daughters: we need to schedule TWO future division events)
        for _ in range(2):
            st = schedule_cell(clone_id, t_birth=t_div, p_birth=p_daughter)
            if st is not None:
                heapq.heappush(event_heap, (st.t_div, uid, st))
                uid += 1

    # Record remaining sample times if any
    record_up_to(t_end)

    return times, counts


# ==========================
# Example usage
# ==========================

if __name__ == "__main__":
    rng = np.random.default_rng()

    # Four clones with different affinity/probability factors phi
    clone_phis = [10., 1., 0.1, 0.01]
    initial_cells_per_clone = [1, 1, 1, 1]

    # Antigen
    rhoA0 = 1.0
    lamA = 0.08  # per hour

    # Presentation dynamics
    kp = 0.5     # per hour per antigen unit
    deltapi = 0.05    # per hour (turnover)

    # Help dynamics
    kh = 0.8     # per hour per T cell unit
    T_const = 1.0
    K_help = 10.0
    n_help = 1.0

    # Division program
    tau0 = 8.0  # hours

    times, counts = simulate_clonal_expansion(
        rng=rng,
        clone_phis=clone_phis,
        initial_cells_per_clone=initial_cells_per_clone,
        rhoA0=rhoA0,
        lamA=lamA,
        kp=kp,
        deltapi=deltapi,
        kh=kh,
        T_const=T_const,
        K_help=K_help,
        n_help=n_help,
        tau0=tau0,
        uptake_during_division=False,
        t_end=96.0,
        sample_dt=0.5,
        max_cells=200_000,
    )
    
    fig_cells, ax_cells = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

    ax_cells.plot(times, counts)
    ax_cells.plot(times, np.exp((np.log(2)/tau0)*(times-20)), '--', color = 'k', label=r'$\tau_0^{-1}$')
    ax_cells.set_xlabel("Time (hours)")
    ax_cells.set_ylabel("Clone counts")
    ax_cells.set_yscale('log')
    ax_cells.set_ylim(0.9, counts.max()*1.2)
    ax_cells.legend([f"Clone {i} (phi={clone_phis[i]})" for i in range(len(clone_phis))])
    fig_cells.savefig(output_plot + '/cells.pdf')
    # Print final clone sizes
    print("Final clone sizes at t_end:", counts[-1])
