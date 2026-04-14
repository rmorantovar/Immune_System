"""
Study N_B_tot (total activated B cells) and L_act (number of activated clones)
across the D~1 crossover.

Uses functions from ep_meanfield_sim.py.
"""

from tkinter.messagebox import IGNORE

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../../library/')
from lib_mf import*


if __name__ == '__main__':
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_Tcell_jam_early/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, delta_A=0.01,
                k_on=1e2*1e6*1e6*24*3600/N_Avg, delta_pi=24., Theta=10.0,
                hill=3.0, sigma=1.0, beta_star=2.5,
                delta_T=0.0, h0=20,
                tau_eng=0.01, b0=2.0, delta_B=0.0,
                DG_min=0.0, DG_max=4.0, M=20,
                Omega_0=1.0, T_lim = 1 , memory = False
    )
    T = 16
    N_ensemble = 20
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    # N_T_values = [1e3] # For testing
    N_T_values = [1e2, 1e3, 1e4, 1e5, 1e20]

    fig_zipf, ax_zipf = plt.subplots(figsize=(6, 5))

    # Storage for summary
    summary = []

    for N_T in N_T_values:
        print(f"Running simulation for N_T={N_T:.1e}")
        ranks = np.linspace(1, 100, 100)
        sizes = np.zeros_like(ranks)
        for i in range(N_ensemble):
            print(f" {i+1}/{N_ensemble}...")
            p = Parameters(**base, N_T0=N_T, lambda_A = 6.)
            # res = run_simulation(p=p, t_span=(0, T), mode='grid')
            res = run_simulation(p=p, t_span=(0, 10), mode='stochastic')

            ranks_i, sizes_i = compute_zipf(res, time_index=-1)
            sizes += sizes_i[:len(ranks)]

            t = res['t']
            N_B_tot = compute_N_B_tot(res)
            L_act = compute_L_act(res)
            t_D = find_t_D(res)

            summary.append({
                'N_T': N_T,
                'lambda_A': 6.0,
                't_D': t_D,
                'N_B_tot_final': N_B_tot[-1],
                'L_act_final': L_act[-1],
                'D_final': res['D'][-1],
            })

            label = f'$N_T=1e{int(np.log10(N_T))}$'

        ax_zipf.loglog(ranks, sizes/N_ensemble, marker='o', ls = '', ms = 4, markeredgecolor='k', label=label, markeredgewidth=0.5)


    # Formatting
    ax_zipf.loglog(ranks, ranks**(-p.sigma*(p.b0/p.lambda_A + 1)/p.beta_star), linestyle = '--', markersize=2, color = my_red, label='Zipf slope')
    ax_zipf.loglog(ranks, ranks**(-2*p.sigma*(p.b0/p.lambda_A + 1)/p.beta_star), linestyle = '--', markersize=2, color = my_blue, label='Zipf slope')
    ax_zipf.set_xlabel('Rank $k$'); ax_zipf.set_ylabel('$N_k$')
    ax_zipf.set_title('Zipf plot (final)')
    ax_zipf.legend(fontsize=7)
    plt.tight_layout()
    fig_zipf.savefig(os.path.join(output_plot, f'Zipf_Tlim-{int(base["T_lim"])}_Mem-{int(p.memory)}.pdf'), dpi=150, bbox_inches='tight')