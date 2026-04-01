"""
Mean-field simulation of T-cell-mediated B-cell clonal expansion.

Simulates the coupled system:
    dN_A/dt = S_A(t) - delta_A * N_A
    dpi/dt  = k_on * psi(DG) * N_A - delta_pi * pi - lambda_B * pi
    dN_B/dt = (lambda_B - delta_B) * N_B
    dN_T/dt = S_T(t) - delta_T * N_T

with:
    lambda_B = (1/h_T + 1/b_0)^{-1}
    h_T = gamma * pi^h/(pi^h + Theta^h) * N_T_free
    N_T_free = N_T / (1 + D)
    D = tau_eng * gamma * sum_i Omega_i * pi^h/(pi^h + Theta^h) * N_B_i * dDG

State variables:
    N_A(t)          : scalar, antigen concentration
    pi(t, DG_i)     : per-cell pMHC for clone i
    N_B(t, DG_i)    : clone size for clone i
    N_T(t)          : scalar, total T cells
"""

import os
import sys
sys.path.append('../../library/')
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass
from mean_field_lib import*


if __name__ == '__main__':
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_dynamics'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    p = Parameters(
        N_A0=1.0,
        lambda_A=6.0,
        delta_A=0.1,
        k_on=1e0*1e6*1e6*24*3600/N_Avg,
        delta_pi=0.1,
        Theta=100.0,
        sigma=1.0,
        beta_star=2.5,
        N_T0=1e4,
        delta_T=0.0,
        gamma=1.0,
        tau_eng=0.1,
        b_0=1.5,
        delta_B=0.0,
        DG_min=0.0,
        DG_max=6.0,
        M=10,
        Omega_0=1.0,
    )
    T = 12
    res = run_simulation(p=p, t_span=(0, T), t_eval=np.linspace(0, T, 1000))

    fig = plot_results(res)
    fig.savefig(os.path.join(output_plot, 'ep_meanfield_results_no_Tcell_limitation_memory.pdf'), dpi=150, bbox_inches='tight')
    print(f"Saved: {os.path.join(output_plot, 'ep_meanfield_results_no_Tcell_limitation_memory.pdf')}")

    # Print some diagnostics
    print(f"\nFront velocity (theory): v = lambda_A / sigma = {p.lambda_A / p.sigma:.3f}")
    Gamma = p.gamma * p.N_T0 * p.k_on * p.N_A0 / (p.delta_pi * p.Theta)
    print(f"Gamma = {Gamma:.4f}")
    print(f"Final demand D = {res['D'][-1]:.3f}")
    print(f"Final N_T_free = {res['N_T_free'][-1]:.3f}")
    print(f"Max clone size = {res['N_B'][:, -1].max():.2e}")