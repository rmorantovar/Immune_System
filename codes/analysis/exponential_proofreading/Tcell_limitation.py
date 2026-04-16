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
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_Tcell_jam_2/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, delta_A=0.01,
                k_on=1e-1*1e6*1e6*24*3600/N_Avg, delta_pi=24., Theta=10000.0,
                hill=3.0, sigma=1.0, beta_star=2.5, K_T = 100000.,
                delta_T=0.0, h0=0.01,
                tau_eng=0.01, b0=2.0, delta_B=0.0,
                DG_min=0.0, DG_max=4.0, M=20,
                Omega_0=1.0, T_lim = 0, memory = False
    )
    T = 16
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    N_T_values = [1e3, 1e6] # For testing
    # N_T_values = [1e2, 1e3, 1e4, 1e5, 1e20]
    

    fig, axes = plt.subplots(6, 1, figsize=(8, 12))

    # Storage for summary
    summary = []

    for N_T in N_T_values:
    # for lam_A in [5.0, 6.0]:
        print(f"Running simulation for N_T={N_T:.1e}")
        p = Parameters(**base, N_T0=N_T, lambda_A = 6.)
        # res = run_simulation(p=p, t_span=(0, T), mode='grid')
        res = run_simulation(p=p, t_span=(0, T), mode='stochastic')

        # print(compute_N_B_tot(res))
        # fig_NT = plot_T_cell_analysis(res)
        # fig_NT.savefig(os.path.join(output_plot, f'meanfield_results_Tlim-{int(p.T_lim)}_Mem-{int(p.memory)}_N_T-1e{int(np.log10(N_T))}.pdf'), dpi=150, bbox_inches='tight')

        # fig = plot_diagnostics(res)
        # fig.savefig(output_plot + f'zipf_Tlim-{int(p.T_lim)}_N_T-1e{int(np.log10(N_T))}.pdf', dpi=150, bbox_inches='tight')

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

        # (a) N_A
        ax = axes[0]
        ax.semilogy(t, res['N_A'], label=label, linewidth = 2, color = antigen_color)

        # (b) pi
        ax = axes[1]
        ax1 = ax.semilogy(t, res['pi'][0, :], label=label, linewidth = 2)
        ax.semilogy(t, res['pi'][-1, :], color=ax1[0].get_color(), alpha = 0.7, label=label, linewidth = 2)
        ax.semilogy(t, res['D'] + 1e-10, color=ax1[0].get_color(), linestyle = '--', linewidth = 2)
        
        # (c) N_T_free / N_T
        ax = axes[2]
        ax.semilogy(t, res['N_T_free'] / N_T, label=label, linewidth = 2)

        # (d) Help rate: h(N_T_free) = h0 * (N_T_free/N_T)^sigma
        ax = axes[3]
        ax.semilogy(t, ((res['lambda_B'][0, :])**(-1) - p.b0**(-1))**(-1), label=label, linewidth = 2)

        # (e) Growth rate: d(ln N_B_tot)/dt
        ax = axes[4]
        ax.plot(t, res['lambda_B'][0, :]/p.b0, label=label, linewidth = 2)
        
        # (f) N_B_tot
        ax = axes[5]
        ax1 = ax.semilogy(t, res['N_B'][0, :], linestyle = '-', label=label, linewidth = 2)
        ax.semilogy(t, res['N_B'][-1, :], linestyle = '-', color=ax1[0].get_color(), alpha = 0.7, linewidth = 2)


    # Formatting
    # axes[0].set_xlabel('Time')
    axes[0].set_ylabel('$N_A$')
    axes[0].set_xticklabels([])

    # axes[1].set_xlabel('Time')
    axes[1].set_ylabel(r'$\pi$')
    axes[1].legend(fontsize=10)
    axes[1].set_xticklabels([])

    # axes[2].set_xlabel('Time')
    axes[2].set_ylabel('$N_T^o/N_T$')
    axes[2].axhline(1.0, color='k', linestyle=':', alpha=0.5)
    axes[2].set_xticklabels([])

    axes[3].axhline(1.0, color='k', linestyle='--', alpha=0.5)
    # axes[3].set_xlabel('Time')
    axes[3].set_ylabel(r'$h$')
    axes[3].set_xticklabels([])

    axes[4].axhline(1.0, color='k', linestyle='--', alpha=0.5)
    # axes[4].set_xlabel('Time')
    axes[4].set_ylabel(r'$\lambda_B$')
    axes[4].set_xticklabels([])

    axes[5].set_xlabel('Time')
    axes[5].set_ylabel('$N_B^{tot}$ and $L_{act}$' )
    # axes[5].set_yscale('log')
    

    plt.suptitle('Effect of $N_T$ on activation dynamics', fontsize=14)
    plt.tight_layout()
    fig.savefig(os.path.join(output_plot, f'comparison_Tlim-{int(base["T_lim"])}_Mem-{int(p.memory)}.pdf'), dpi=150, bbox_inches='tight')

    # Print summary
    print(f'\n{"N_T":>8} {"t_D":>8} {"N_B_tot":>12} {"L_act":>10} {"D_final":>10}')
    for s in summary:
        print(f'{s["N_T"]:>8} {s["t_D"]:>8.2f} {s["N_B_tot_final"]:>12.1f} '
            f'{s["L_act_final"]:>10.1f} {s["D_final"]:>10.1f}')