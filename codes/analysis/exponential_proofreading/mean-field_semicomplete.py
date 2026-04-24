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
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_semicomplete/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, delta_A=4.0, lambda_A = 6.,
                k_on=1e2*1e6*1e6*24*3600/N_Avg, delta_pi=0.1, Theta=1000.0,
                hill=2.0, sigma=1.0, beta_star=2.5, K_T = 1e5,
                delta_T=0.00, h0=1000.0, Tcell_growth_factor=2.0,
                tau_eng=0.1, b0=2.0, delta_B=0.00,
                DG_min=0.0, DG_max=3.0, M=20,
                Omega_0=1.0, T_lim = True
    )
    T = 14
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    N_T_values = [1e4] # For testing
    # N_T_values = [1e2, 1e3, 1e4, 1e5, 1e20]
    
    for memory in [0]:
        fig, axes = plt.subplots(5, 1, figsize=(8, 15))
        print(f"-- Running simulation for memory={memory} --")
        base['memory'] = memory
        for N_T in N_T_values:
        # for lam_A in [5.0, 6.0]:
            print(f"Running simulation for N_T={N_T:.1e}")
            p = Parameters(**base, N_T0=N_T)
            # res = run_simulation(p=p, t_span=(0, T), mode='grid')
            res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='stochastic')

            t = res['t']
            
            label = f'$N_T=1e{int(np.log10(N_T))}$'

            # (a) N_A
            ax = axes[0]
            ax.semilogy(t, res['N_A'], label=label, linewidth = 2, color = antigen_color)

            # (b) pi
            ax = axes[1]
            ax1 = ax.semilogy(t, res['pi'][0, :], label=label, linewidth = 2)
            ax.semilogy(t, res['pi'][40, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            ax.semilogy(t, 1e12*np.exp(-(p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # (c) lambda_B
            ax = axes[2]
            ax1 = ax.plot(t, p.b0*res['N_Ba'][0, :]/(res['N_Bo'][0, :] + res['N_Ba'][0, :]), label=label, linewidth = 2)
            ax.plot(t, p.b0*res['N_Ba'][40, :]/(res['N_Bo'][40, :] + res['N_Ba'][40, :]), label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.5)
            ax.plot(t, ((p.k_on * res['pi'][0, :]**p.hill * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1), linestyle='--', color=ax1[0].get_color(), linewidth = 2)
            ax.axhline(p.b0, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # (d) T
            ax = axes[3]
            ax1 = ax.semilogy(t, res['N_To'] + res['N_Ta'] , label=label, linewidth = 3)
            ax.semilogy(t, res['N_To'], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
            # ax.semilogy(t, res['N_Ta'], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())

            # (e) B
            ax = axes[4]
            N_B_total = np.sum(res['N_Bo'] + res['N_Ba'] , axis=0)
            N_B_total = N_B_total - N_B_total[0]
            ax1 = ax.semilogy(t, res['N_Bo'][0, :] + res['N_Ba'][0, :], label=label, linewidth = 3)
            ax.semilogy(t, res['N_Bo'][40, :] + res['N_Ba'][40, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())
            # ax.semilogy(t, res['N_Bo'][0, :], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
            # ax.semilogy(t, res['N_Ba'][0, :], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())
            # ax.semilogy(t, N_B_total, label=label, linewidth = 3, linestyle='-', alpha=0.8, color=ax1[0].get_color())
            ax.semilogy(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

        # Formatting
        # axes[0].set_xlabel('Time')
        axes[0].set_ylabel('$N_A$', fontsize = 16)
        axes[0].set_xticklabels([])
        axes[0].set_ylim(bottom = 1e-1)
        axes[0].tick_params(labelsize=14)

        # axes[1].set_xlabel('Time')
        axes[1].set_ylabel(r'$\pi$', fontsize = 16)
        axes[1].set_xticklabels([])
        axes[1].set_ylim(top = 10*np.max(res['pi'][0, :]))

        # axes[1].set_xlabel('Time')
        axes[2].set_ylabel(r'$\lambda_B$', fontsize = 16)
        axes[2].set_xticklabels([])
        axes[2].tick_params(labelsize=14)

        # axes[2].set_xlabel('Time')
        axes[3].set_ylabel('T cells', fontsize = 16)
        axes[3].axhline(1.0, color='k', linestyle=':', alpha=0.5)
        axes[3].set_xticklabels([])
        axes[3].set_yscale('log')
        axes[3].set_ylim(bottom = 5e-1*p.N_T0)
        # axes[3].legend(fontsize=10)
        axes[3].tick_params(labelsize=14)

        axes[4].axhline(1.0, color='k', linestyle='--', alpha=0.5)
        # axes[4].set_xlabel('Time')
        axes[4].set_ylabel('B cells', fontsize = 16)
        # axes[4].set_xticklabels([])
        axes[4].set_xlabel('Time', fontsize = 16)
        axes[4].set_ylim(bottom = 5e-1, top = np.max(res['N_Bo'][0, :])*10)
        axes[4].tick_params(labelsize=14)
        
        plt.suptitle('Effect of $N_T$ on activation dynamics', fontsize=14)
        plt.tight_layout()
        fig.savefig(os.path.join(output_plot, f'results_Mem-{int(p.memory)}.pdf'), dpi=150, bbox_inches='tight')
