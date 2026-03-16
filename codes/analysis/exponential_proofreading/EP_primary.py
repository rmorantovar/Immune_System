
"""
Simulation with:
- Antigen: NA(t) = exp(lambda_A * t)
- π dynamics (per affinity K):
    dπ/dt = u_on(t) * g(K) - lambda_B * π
  where u_on(t) = k_on_eff * NA(t)
- B cell divisions happen while π > threshold:
    dB/dt = r_div * B   if π > threshold else 0
  (Note: with constant lambda_B, π is *not* explicitly diluted by division; if you
   want division-driven dilution too, say so and we’ll add + r_div*ln2*π term.)
Outputs full time series for t, NA(t), π(t,K), B(t,K) so you can plot.
"""

import sys
sys.path.append('../../library/')
from funcs import*
plt.rcParams['text.usetex'] = True
from dataclasses import dataclass
from typing import Callable, Sequence, Tuple, Dict
from scipy.integrate import solve_ivp

@dataclass
class Params:
    lambda_A: float                 # antigen exponential growth rate
    k_on_eff: float                 # u_on(t) = k_on_eff * NA(t)
    lambda_B: float                 # constant decay/dilution rate for pi
    delta: float                    # decay rate for pi even when not dividing
    pi_threshold: float             # division condition: pi > threshold
    t_span: Tuple[float, float]     # (t0, tf)
    t_eval: np.ndarray              # time points for output

# def dNAdtNaive(t, N, lambda_A, lambda_B):
#     pb = (1+(1e-9/(1e6*60*60*24*np.exp(2.0*(t))/N_Avg)))**(-1)
#     return (lambda_A * (1 - pb) - 2*pb) * N

# def NA(t: float, lambda_A: float) -> float:

#     # Initial condition
#     NA0 = 1.0
#     # Solve the ODE over the time span of your data
#     solNaive = solve_ivp(dNAdtNaive, t_span=p.t_span, y0=[NA0], t_eval=p.t_eval, args=(lambda_A, p.lambda_B), method="RK45", rtol=1e-6, atol=1e-9)
#     # Result
#     NA_real = solNaive.y[0]  # solution N(t) evaluated at t_vals
#     # Interpolate to get value at time t
#     return float(np.interp(t, p.t_eval, NA_real))

# def u_on(t: float, p: Params) -> float:
#     return p.k_on_eff * NA(t, p.lambda_A)


def simulate(
        
    Ks: Sequence[float],
    g: Callable[[float], float],
    p: Params,
    pi0: float = 0.0,
    B0: float = 1.0,
) -> Dict[str, np.ndarray]:
    """
    Returns a dict with:
      - t: (T,)
      - NA: (T,)
      - Ks: (nK,)
      - pi: (T, nK)
      - B:  (T, nK)
      - dividing: (T, nK) boolean mask where pi > threshold
    """
    Ks = np.asarray(Ks, dtype=float)
    n = Ks.size
    gK = np.array([g(K) for K in Ks], dtype=float)


    def rhs(t, y):
        NA = y[0]
        pi = y[1:1+n]
        B  = y[1+n:1+2*n]

        # compute pb from NA (same formula you use in dNAdtNaive)
        pb = (1 + (1e-10/(1e6*60*60*24*np.exp(2.0*t)/N_Avg)))**(-1)  # or whatever dependence you intend
        dNA = (p.lambda_A * (1 - pb) - 3*pb) * NA

        u_on = p.k_on_eff * NA
        dividing = pi > p.pi_threshold
        dpi = np.where(dividing, u_on * gK - p.delta * pi - p.lambda_B*pi, u_on * gK - p.delta * pi)
        dB  = np.where(dividing, np.log(2)*p.lambda_B*B, 0.0)

        return np.concatenate([[dNA], dpi, dB])

    # y0 = np.concatenate([np.full(n, pi0, dtype=float), np.full(n, B0, dtype=float)])
    NA0 = 1.0
    y0 = np.concatenate([[NA0], np.full(n, pi0), np.full(n, B0)])

    sol = solve_ivp(
        rhs,
        t_span=p.t_span,
        y0=y0,
        t_eval=p.t_eval,
        # method="RK45",
        # rtol=1e-6,
        # atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    t = sol.t
    NA_t = sol.y[0, :]          # instead of recomputing via NA(t, ...)
    pi_tK = sol.y[1:1+n, :].T # (T, nK)
    B_tK  = sol.y[1+n:1+2*n, :].T  # (T, nK)

    # NA_t = np.array([NA(t, p.lambda_A) for t in p.t_eval])  # vectorized NA(t)

    return {
        "t": t,
        "NA": NA_t,
        "Ks": Ks,
        "pi": pi_tK,
        "B": B_tK
    }

# -----------------------
# Example usage + plotting
# -----------------------
if __name__ == "__main__":
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/model_dynamics_time'
    os.makedirs(output_plot, exist_ok=True),

    # Example g(K): power law g(K)=K^{-sigma}
    
    def g_power(K: float) -> float:
        sigma = 1
        K_s = 1/((60*7)*3600*24)  # relevant scale for K 
        # K_s = 1e-8
        return (K_s/(K_s+K)) ** (sigma)

    Ks = [8.0e-8, 5.0e-7, 3.0e-6]  # example affinities

    t0, tf = 0.0, 10.5
    t_eval = np.linspace(t0, tf, 101)

    p = Params(
        lambda_A=8.5,
        k_on_eff=1e6*1e6*24*3600/N_Avg,
        lambda_B=8.5*0.4,
        delta=.5,
        pi_threshold=10.0,
        t_span=(t0, tf),
        t_eval=t_eval,
    )

    out = simulate(Ks, g_power, p, pi0=0.0, B0=1.0)

    t = out["t"]
    NA_t = out["NA"]
    pi = out["pi"]
    B  = out["B"]
    Ks = out["Ks"]

    
    # Plot antigen
    fig_antigen, ax_antigen = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_antigen, ax_antigen = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
    ax_antigen.plot(t, NA_t, color = antigen_color, linewidth = 5)
    
    ax_antigen.tick_params(axis='both', which='major', labelsize=28)
    # ax.legend(loc = 2, fontsize = 14, title = r'$K$', title_fontsize = 16)
    ax_antigen.set_yscale('log')
    ax_antigen.set_xlim(left = 0, right = 8)
    fig_antigen.savefig(output_plot + '/antigen_primary_.pdf')
    
    
    my_plot_layout(ax=ax_antigen, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    # ax_antigen.set_xlabel(r"$t$")
    # ax_antigen.set_ylabel(r"$N_A(t)$")
    # ax_antigen.set_title("Antigen")
    ax_antigen.set_ylim(top = 1e11, bottom = 8e-1)
    ax_antigen.set_yscale("log")
    # ax_antigen.set_xticks([])
    fig_antigen.savefig(output_plot + '/antigen_primary.pdf')
    plt.close(fig_antigen)


    # Plot pi trajectories
    my_colors = [my_green_a, my_green_b, my_green_c]
    fig_Pi, ax_Pi =plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.2, 'top': 0.96})
    fig_Pi, ax_Pi =plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
    for i, K in enumerate(out["Ks"]):
        ax_Pi.plot(t, pi[:, i], label=f"K={K}", linewidth = 5, color = my_colors[i])
    # ax_Pi.plot(t, np.exp(6*(t-1))*10, color = 'grey', ls = '--', lw = 2)
    ax_Pi.axhline(p.pi_threshold, linestyle="--", label="threshold", color = "k")
    # my_plot_layout(ax=ax_Pi, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    ax_Pi.tick_params(axis='both', which='major', labelsize=28)
    # ax_Pi.set_xlabel(r"$t$")
    # ax_Pi.set_ylabel(r"$\pi(t,K)$")
    # ax_Pi.set_title("Internalized antigen per cell")
    ax_Pi.set_yscale("log")
    ax_Pi.set_ylim(bottom=0.99, top = 5e3)  # set a reasonable lower limit for log scale
    # ax_Pi.legend()
    # ax_Pi.set_xticks([])
    ax_Pi.set_xlim(left = 0, right = 8)
    fig_Pi.savefig(output_plot + '/pi_primary.pdf')
    plt.close(fig_Pi)


    # Plot B trajectories
    fig_N_b, ax_N_b = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.2, 'top': 0.96})
    fig_N_b, ax_N_b = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
    for i, K in enumerate(out["Ks"]):
        ax_N_b.plot(t, B[:, i], label=f"K={K}", linewidth = 5, color = my_colors[i])
    # ax_N_b.plot(t, np.exp(6*0.4*(t-1)), color = 'grey', ls = '--', lw = 2)
    # my_plot_layout(ax=ax_N_b, yscale = 'log', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    # ax_N_b.set_xlabel("t")
    # ax_N_b.set_ylabel(r"$B(t,K)$")
    # ax_N_b.set_title("B cell population (divides while pi > threshold)")
    ax_N_b.tick_params(axis='both', which='major', labelsize=28)
    ax_N_b.set_yscale("log")
    ax_N_b.set_ylim(bottom = 0.99, top=1e5)  # set a reasonable lower limit for log scale
    ax_N_b.set_xlim(left = 0, right = 8)
    # ax_N_b.legend()
    fig_N_b.savefig(output_plot + '/B_primary.pdf')
    plt.close(fig_N_b)