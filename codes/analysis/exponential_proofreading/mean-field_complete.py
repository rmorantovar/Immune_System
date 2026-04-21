"""
Study N_B_tot (total activated B cells) and L_act (number of activated clones)
across the D~1 crossover.

Uses functions from ep_meanfield_sim.py.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../../library/')
from lib_mf import*


if __name__ == '__main__':
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_complete/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, delta_A=0.01, lambda_A = 6.,
                k_on=1e2*1e6*1e6*24*3600/N_Avg, delta_pi=24., Theta=1000.0,
                hill=1.0, sigma=1.0, beta_star=2.5, K_T = 1e1,
                delta_T=0.00, h0=0.01,
                tau_eng=0.1, b0=2.0, delta_B=0.00,
                DG_min=0.0, DG_max=3.5, M=20,
                Omega_0=1.0, memory = False
    )
    T = 16
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    N_T_values = [1e2] # For testing
    # N_T_values = [1e2, 1e3, 1e4, 1e5, 1e20]
    
    for T_lim in [0, 1]:
        fig, axes = plt.subplots(4, 1, figsize=(8, 12))
        print(f"-- Running simulation for T_lim={T_lim} --")
        base['T_lim'] = T_lim
        for N_T in N_T_values:
        # for lam_A in [5.0, 6.0]:
            print(f"Running simulation for N_T={N_T:.1e}")
            p = Parameters(**base, N_T0=N_T)
            # res = run_simulation(p=p, t_span=(0, T), mode='grid')
            res = run_simulation_complete(p=p, t_span=(0, T), mode='stochastic')

            t = res['t']
            
            label = f'$N_T=1e{int(np.log10(N_T))}$'

            # (a) N_A
            ax = axes[0]
            ax.semilogy(t, res['N_A'], label=label, linewidth = 2, color = antigen_color)

            # (b) pi
            ax = axes[1]
            ax1 = ax.semilogy(t, res['pi'][0, :], label=label, linewidth = 2)

            # (c) T
            ax = axes[2]
            ax1 = ax.semilogy(t, res['N_To'], label=label, linewidth = 2)
            ax.semilogy(t, res['N_Ta'], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
            ax.semilogy(t, res['N_BT'][0, :], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())

            # (d) B
            ax = axes[3]
            ax1 = ax.semilogy(t, res['N_Bo'][0, :], label=label, linewidth = 2)
            ax.semilogy(t, res['N_Ba'][0, :], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
            ax.semilogy(t, res['N_BT'][0, :], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())


        # Formatting
        # axes[0].set_xlabel('Time')
        axes[0].set_ylabel('$N_A$')
        axes[0].set_xticklabels([])

        # axes[1].set_xlabel('Time')
        axes[1].set_ylabel(r'$\pi$')
        axes[1].set_xticklabels([])

        # axes[2].set_xlabel('Time')
        axes[2].set_ylabel('T cells')
        axes[2].axhline(1.0, color='k', linestyle=':', alpha=0.5)
        axes[2].set_xticklabels([])
        # axes[2].legend(fontsize=10)

        axes[3].axhline(1.0, color='k', linestyle='--', alpha=0.5)
        # axes[3].set_xlabel('Time')
        axes[3].set_ylabel('B cells')
        axes[3].set_xticklabels([])

        
        plt.suptitle('Effect of $N_T$ on activation dynamics', fontsize=14)
        plt.tight_layout()
        fig.savefig(os.path.join(output_plot, f'comparison_Tlim-{int(base["T_lim"])}_Mem-{int(p.memory)}.pdf'), dpi=150, bbox_inches='tight')
