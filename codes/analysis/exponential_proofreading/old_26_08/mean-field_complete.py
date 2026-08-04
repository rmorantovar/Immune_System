"""
Complete model: explicit T-B conjugate intermediate state N_BT.

State variables per clone:
    N_Bo  : resting (free) B cells
    N_BT  : T-B conjugates (B cell engaged with a T cell)
    N_Ba  : activated (dividing) B cells
    N_To  : free T cells
    N_Ta  : activated T cells

Extends the semicomplete model by resolving the engagement step explicitly.
k_on drives conjugate formation: dN_BT/dt = k_on*pi^h*N_Bo*N_To/(N_To+K_T) - h0*N_BT
h0 is the conjugate release/activation rate (not a direct help rate).
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../../library/')
from lib_mf import *


if __name__ == '__main__':
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_complete/'
    os.makedirs(output_plot, exist_ok=True)

    base = dict(N_A0=1.0, delta_A=3.0, lambda_A=6.,
                k_on=1e2*1e6*1e6*24*3600/N_Avg, delta_pi=0.0, Theta=1000.0,
                hill=3.0, sigma=1.0, beta_star=2.5, K_T=1e4,
                delta_T=0.00, h0=100.0,
                tau_eng=0.1, b0=2.0, delta_B=0.00,
                DG_min=0.0, DG_max=3.0, M=20,
                Omega_0=1.0, T_lim=True,
    )
    T = 14

    N_T_values = [1e1]  # For testing
    # N_T_values = [1e2, 1e3, 1e4, 1e5, 1e20]

    for memory in [0]:
        fig, axes = plt.subplots(6, 1, figsize=(8, 18))
        print(f"-- Running simulation for memory={memory} --")
        base['memory'] = memory
        for N_T in N_T_values:
            print(f"Running simulation for N_T={N_T:.1e}")
            p = Parameters(**base, N_T0=N_T)
            res = run_simulation_complete(p=p, t_span=(0, T), mode='stochastic')

            print(f"  M = {res['M']} clones")
            t = res['t']
            label = f'$N_T=1e{int(np.log10(N_T))}$'

            N_B_tot_i = res['N_Bo'] + res['N_BT'] + res['N_Ba']

            # (a) N_A
            ax = axes[0]
            ax.semilogy(t, res['N_A'], label=label, linewidth=2, color=antigen_color)

            # (b) pi — show low/mid/high DG clones
            ax = axes[1]
            ax1 = ax.semilogy(t, res['pi'][0, :],  label=label, linewidth=2)
            ax.semilogy(t, 1e10*np.exp(-p.b0*t), linewidth=1, linestyle='--', color='grey', alpha=0.8)

            # (c) lambda_B/b0 = N_Ba / (N_Bo + N_BT + N_Ba) — activated fraction
            ax = axes[2]
            for idx, alpha in zip([0, -1], [1.0, 0.3]):
                denom = np.where(N_B_tot_i[idx, :] > 0, N_B_tot_i[idx, :], 1.0)
                rate  = res['N_Ba'][idx, :] / denom / p.b0
                if idx == 0:
                    ax1 = ax.plot(t, rate, linewidth=2, label=label)
                else:
                    ax.plot(t, rate, linewidth=2, alpha=alpha, color=ax1[0].get_color())
            # quasi-steady-state approximation for lambda_B
            ax.plot(t,
                    ((p.k_on * res['pi'][0, :]**p.hill * res['N_To'] / (res['N_To'] + p.K_T))**(-1) + p.b0**(-1))**(-1) / p.b0,
                    linestyle='--', color=ax1[0].get_color(), linewidth=2)
            ax.axhline(1.0, linewidth=1, linestyle='--', color='grey', alpha=0.8)

            # (d) T-B conjugates N_BT per clone — unique to complete model
            ax = axes[3]
            ax1 = ax.semilogy(t, res['N_BT'][0, :] + 1e-10, label=label, linewidth=2)
            ax.semilogy(t, res['N_BT'][-1, :] + 1e-10, linewidth=2, alpha=0.3, color=ax1[0].get_color())

            # (e) T cells: total (solid), free N_To (dashed), activated N_Ta (dotted)
            ax = axes[4]
            N_T_tot = res['N_To'] + res['N_Ta'] + np.sum(res['N_BT'], axis=0)
            ax1 = ax.semilogy(t, N_T_tot,    label=label, linewidth=3)
            ax.semilogy(t, res['N_To'],      linewidth=2, linestyle='--', color=ax1[0].get_color())
            ax.semilogy(t, res['N_Ta']+1e-10, linewidth=2, linestyle=':',  color=ax1[0].get_color())

            # (f) B cells per clone: total (solid), N_Bo (dashed), N_Ba (dotted), N_BT (black)
            ax = axes[5]
            ax1 = ax.semilogy(t, N_B_tot_i[0, :],  label=label, linewidth=3)
            ax.semilogy(t, res['N_Bo'][0, :],       linewidth=2, linestyle='--', color=ax1[0].get_color())
            ax.semilogy(t, res['N_Ba'][0, :]+1e-10, linewidth=2, linestyle=':',  color=ax1[0].get_color())
            ax.semilogy(t, res['N_BT'][0, :]+1e-10, linewidth=2, linestyle='-',  color='k', alpha=0.6)
            ax.semilogy(t, 1e-2*np.exp(p.b0*t),    linewidth=1, linestyle='--', color='grey', alpha=0.8)

        # ── Formatting ────────────────────────────────────────────────────────
        axes[0].set_ylabel('$N_A$', fontsize=16)
        axes[0].set_xticklabels([])
        axes[0].set_ylim(bottom=1e-1)
        axes[0].tick_params(labelsize=14)

        axes[1].set_ylabel(r'$\pi$', fontsize=16)
        axes[1].set_xticklabels([])
        axes[1].set_ylim(bottom=1e-4, top=10*np.max(res['pi'][0, :]))
        axes[1].tick_params(labelsize=14)

        axes[2].set_ylabel(r'$\lambda_B/b_0$', fontsize=16)
        axes[2].set_xticklabels([])
        axes[2].tick_params(labelsize=14)

        axes[3].set_ylabel('$N_{BT}$', fontsize=16)
        axes[3].set_xticklabels([])
        axes[3].tick_params(labelsize=14)

        axes[4].set_ylabel('T cells', fontsize=16)
        axes[4].axhline(1.0, color='k', linestyle=':', alpha=0.5)
        axes[4].set_xticklabels([])
        axes[4].set_yscale('log')
        axes[4].set_ylim(bottom=1e-2)
        axes[4].tick_params(labelsize=14)

        axes[5].axhline(1.0, color='k', linestyle='--', alpha=0.5)
        axes[5].set_ylabel('B cells', fontsize=16)
        axes[5].set_xlabel('Time', fontsize=16)
        axes[5].set_ylim(bottom=1e-2, top=np.max(res['N_Bo'][0, :]) * 10)
        axes[5].tick_params(labelsize=14)

        plt.suptitle('Complete model: T-B conjugate dynamics', fontsize=14)
        plt.tight_layout()
        fig.savefig(os.path.join(output_plot, f'results_Mem-{int(p.memory)}.pdf'),
                    dpi=150, bbox_inches='tight')
