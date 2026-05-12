"""
Complete model: explicit T-B conjugate intermediate state N_BT.

State variables per clone:
    N_Bo  : resting (free) B cells
    N_Ba  : activated (dividing) B cells
    N_To  : free T cells
    N_Ta  : activated T cells

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../../library/')
from lib_mf import*
plt.rcParams['text.usetex'] = True


if __name__ == '__main__':
    output_plot = '/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_semicomplete/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, delta_A=6.0, lambda_A = 6.,
                k_on=1e2*1e6*1e6*24*3600/N_Avg, delta_pi=0.1, Theta=1000.0,
                hill=2.0, beta_star=2.5, K_T = 1e5,
                delta_T=0.00, h0=1e-2, Tcell_growth_factor=2.0,
                tau_eng=0.1, b0=2.0, delta_B=0.00,
                DG_min=0.0, DG_max=3.0, M=20,
                Omega_0=1.0, T_lim = True
    )
    T = 12
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    N_T_values = [1e3, 1e5] # For testing
    # N_T_values = [1e2, 1e3, 1e4, 1e5, 1e20]
    fig_Z, ax_Z = plt.subplots(figsize=(8, 5))
    fig_N, ax_N = plt.subplots(figsize=(8, 5))
    color_sigma = [my_blue, my_purple, my_gold]
    for i_m, memory in enumerate([0]):
        print(f"-- Running simulation for memory={memory} --")
        base['memory'] = memory
        N_T = N_T_values[i_m]
        for i_sigma, sigma in enumerate([1.0, 3.0, 5.0]):
            figs, axes = plt.subplots(6, 1, figsize=(8, 18))
            print(f"Running simulation for sigma={sigma}")
            base['sigma'] = sigma
            p = Parameters(**base, N_T0=N_T)
            # res = run_simulation(p=p, t_span=(0, T), mode='grid')
            res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='stochastic', seed=0)

            print(res['M'])
            t = res['t']
            
            label = f'$N_T=1e{int(np.log10(N_T))}$'

            # (a) N_A
            axi = axes[0]
            axi.semilogy(t, res['N_A'], label=label, linewidth = 2, color = antigen_color)

            # (b) pi
            axi = axes[1]
            ax1 = axi.semilogy(t, res['pi'][0, :], label=label, linewidth = 2)
            axi.semilogy(t, res['pi'][40, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            axi.semilogy(t, res['pi'][-1, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            axi.semilogy(t, 1e10*np.exp(-(p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # (c) lambda_B
            axi = axes[2]
            ax1 = axi.plot(t, p.b0*res['N_Ba'][0, :]/(res['N_Bo'][0, :] + res['N_Ba'][0, :])/p.b0, label=label, linewidth = 2)
            axi.plot(t, p.b0*res['N_Ba'][40, :]/(res['N_Bo'][40, :] + res['N_Ba'][40, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.5)
            axi.plot(t, p.b0*res['N_Ba'][-1, :]/(res['N_Bo'][-1, :] + res['N_Ba'][-1, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.2)
            axi.plot(t, ((p.k_on * res['pi'][0, :]**p.hill * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
            axi.axhline(p.b0/p.b0, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
    
            # (d) T
            axi = axes[3]
            ax1 = axi.semilogy(t, res['N_To'] + res['N_Ta'] , label=label, linewidth = 3)
            axi.semilogy(t, res['N_To'], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())

            # (e) B
            axi = axes[4]
            N_B_total = np.sum(res['N_Bo'] + res['N_Ba'] , axis=0)
            N_B_total = N_B_total - N_B_total[0]
            ax1 = axi.semilogy(t, N_B_total, label=label, linewidth = 3)
            axi.semilogy(t, res['N_Bo'][0, :] + res['N_Ba'][0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
            axi.semilogy(t, res['N_Bo'][40, :] + res['N_Ba'][40, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            axi.semilogy(t, res['N_Bo'][-1, :] + res['N_Ba'][-1, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            # axi.semilogy(t, res['N_Bo'][0, :], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
            # axi.semilogy(t, res['N_Ba'][0, :], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())
            # axi.semilogy(t, N_B_total, label=label, linewidth = 3, linestyle='-', alpha=0.8, color=ax1[0].get_color())
            axi.semilogy(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # (f) Potency
            axi = axes[5]
            Z_B = (res['N_Bo'] + res['N_Ba'])*np.exp(-res['DG'][:, None])
            Z_B_total = np.sum(Z_B, axis=0)
            ax1 = axi.semilogy(t, Z_B_total, label=label, linewidth = 3)
            axi.semilogy(t, Z_B[0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
            axi.semilogy(t, Z_B[40, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            axi.semilogy(t, Z_B[-1, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            axi.semilogy(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            

            # Clone size shared axis
            ax1 = ax_N.semilogy(t, N_B_total, label=label, linewidth = 3, color=color_sigma[i_sigma])
            # ax_N.semilogy(t, res['N_Bo'][0, :] + res['N_Ba'][0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
            # ax_N.semilogy(t, res['N_Bo'][40, :] + res['N_Ba'][40, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            # ax_N.semilogy(t, res['N_Bo'][-1, :] + res['N_Ba'][-1, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            ax_N.semilogy(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # Potency shared axis
            N_B = res['N_Bo'] + res['N_Ba']
            activated = N_B > 1.5
            Z_B = N_B*np.exp(-res['DG'][:, None])*activated
            Z_B_total = np.sum(Z_B, axis=0)
            if memory:
                color = my_blue
            else:
                color = my_red
            ax1 = ax_Z.semilogy(t, Z_B_total, label=label, linewidth = 2, color=color_sigma[i_sigma])
            # ax_Z.semilogy(t, Z_B[0, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())
            ax_Z.semilogy(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # Formatting
            # axes[0].set_xlabel('Time')
            axes[0].set_ylabel('$N_A$', fontsize = 16)
            axes[0].set_xticklabels([])
            axes[0].set_ylim(bottom = 1e-1)
            axes[0].set_xlim(0, T)
            axes[0].tick_params(labelsize=14)

            # axes[1].set_xlabel('Time')
            axes[1].set_ylabel(r'$\pi$', fontsize = 16)
            axes[1].set_xticklabels([])
            axes[1].set_ylim(bottom = 1e-4, top = 10*np.max(res['pi'][0, :]))
            axes[1].set_xlim(0, T)

            # axes[2].set_xlabel('Time')
            axes[2].set_ylabel(r'$\lambda_B/b_o$', fontsize = 16)
            axes[2].set_xticklabels([])
            axes[2].set_xlim(0, T)
            axes[2].tick_params(labelsize=14)

            # axes[3].set_xlabel('Time')
            axes[3].set_ylabel('T cells', fontsize = 16)
            axes[3].axhline(1.0, color='k', linestyle=':', alpha=0.5)
            axes[3].set_xticklabels([])
            axes[3].set_xlim(0, T)
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
            axes[4].set_xlim(0, T)
            axes[4].set_yscale('log')
            axes[4].tick_params(labelsize=14)

            # axes[5].set_xlabel('Time')
            axes[5].set_ylabel('Potency, $Z$', fontsize = 16)
            # axes[5].set_xticklabels([])
            axes[5].set_xlabel('Time', fontsize = 16)
            axes[5].set_ylim(bottom = 5e-1, top = 1e11)
            axes[5].set_xlim(0, T)
            axes[5].set_yscale('log')
            axes[5].tick_params(labelsize=14)
            
            plt.suptitle('Effect of $N_T$ on activation dynamics', fontsize=14)
            plt.tight_layout()
            figs.savefig(os.path.join(output_plot, f'results_sigma-{p.sigma}_Mem-{int(p.memory)}.pdf'), dpi=150, bbox_inches='tight')


    ax_N.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    # ax_N.set_xlabel('Time')
    ax_N.set_ylabel(r'Yield, $\bar{N}$', fontsize = 16)
    # ax_N.set_xticklabels([])
    ax_N.set_xlabel('Time', fontsize = 16)
    ax_N.set_ylim(bottom = 5e-1, top = 1e7)
    ax_N.set_xlim(0, T)
    ax_N.set_yscale('log')
    ax_N.tick_params(labelsize=14)
    
    plt.tight_layout()
    fig_N.savefig(os.path.join(output_plot, f'yield_summary.pdf'), dpi=150, bbox_inches='tight')

    ax_Z.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    # ax_Z.set_xlabel('Time')
    ax_Z.set_ylabel('Potency, $Z$', fontsize = 16)
    # ax_Z.set_xticklabels([])
    ax_Z.set_xlabel('Time', fontsize = 16)
    ax_Z.set_ylim(bottom = 5e-1, top = 1e7)
    ax_Z.set_xlim(0, T)
    ax_Z.set_yscale('log')
    ax_Z.tick_params(labelsize=14)
    
    plt.tight_layout()
    fig_Z.savefig(os.path.join(output_plot, f'potency_summary.pdf'), dpi=150, bbox_inches='tight')
