"""
Mean-field simulation of T-cell-mediated B-cell clonal expansion.

Simulates the coupled system:
    dN_A/dt = S_A(t) - delta_A * N_A
    dpi/dt  = k_on * psi(DG) * N_A - delta_pi * pi - lambda_B * pi
    dN_B/dt = (lambda_B - delta_B) * N_B
    dN_T/dt = S_T(t) - delta_T * N_T

with:
    lambda_B = (1/h_T + 1/b0)^{-1}
    h_T = h0 * pi^h/(pi^h + Theta^h) * N_T_free
    N_T_free = N_T / (1 + D)
    D = tau_eng * h0 * sum_i Omega_i * pi^h/(pi^h + Theta^h) * N_B_i * dDG

State variables:
    N_A(t)          : scalar, antigen concentration
    pi(t, DG_i)     : per-cell pMHC for clone i
    N_B(t, DG_i)    : clone size for clone i
    N_T(t)          : scalar, total T cells
"""

import os
import sys
sys.path.append('../../library/')
from lib_mf import*


if __name__ == '__main__':
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_dynamics/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    p = Parameters(N_A0=1.0, lambda_A=6.0, delta_A=0.01,
                   k_on=1e0*1e6*1e6*24*3600/N_Avg, delta_pi=24., Theta=10.0, 
                   hill=3.0, sigma=1.0, beta_star=2.5, 
                   N_T0=1e10, delta_T=0.0, h0=10, 
                   tau_eng=0.5, b0=2.0, delta_B=0.0, 
                   DG_min=0.0, DG_max=5.0,M=200, 
                   Omega_0=1.0, T_lim = 0, memory = True
    )
    # --- Grid mode ---
    print("Running grid mode...")
    res_grid = run_simulation(p=p, t_span=(0, 10), mode='grid')
    fig = plot_diagnostics(res_grid)
    fig.savefig(output_plot + 'sim_grid.png', dpi=150, bbox_inches='tight')
    print(f"  Grid: M={res_grid['M']} bins")
    print(f"  Final L_act = {compute_L_act(res_grid)[-1]:.1f}")
    print(f"  Final N1 = {compute_N1(res_grid)[-1]:.2e}")
 
    # --- Stochastic mode ---
    print("\nRunning stochastic mode...")
    res_stoch = run_simulation(p=p, t_span=(0, 10), mode='stochastic', DG_max_sim=5.0, seed=42)
    fig = plot_diagnostics(res_stoch)
    fig.savefig(output_plot + 'sim_stochastic.png', dpi=150, bbox_inches='tight')
    print(f"  Stochastic: {res_stoch['M']} individual clones")
    print(f"  Final L_act = {compute_L_act(res_stoch)[-1]:.1f}")
    print(f"  Final N1 = {compute_N1(res_stoch)[-1]:.2e}")
 
    # --- Zipf comparison ---
    fig = plot_zipf_comparison(res_grid, res_stoch)
    fig.savefig(output_plot + 'sim_comparison.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved: {output_plot}sim_grid.png, {output_plot}sim_stochastic.png, {output_plot}sim_comparison.png")