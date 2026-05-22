"""
Complete model: non-explicit T-B conjugate intermediate state N_BT.

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

project = 'exponential_proofreading'
model = 'meanfield'
submodel = 'semicomplete'
subproject = 'potency'
subsubproject = 'exploratory'

if __name__ == '__main__':
    output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}/{subproject}/{subsubproject}/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, lambda_A = 6.0, delta_A=4.0, sigma= 1.0,
                k_on=1e0*1e6*1e6*24*3600/N_Avg, delta_pi=0.1,
                hill=1.0, beta_star=2.5, K_T = 1e4,
                delta_T=0.00, Tcell_growth_factor=2.0,
                tau_eng=0.1, b0=2.0, delta_B=0.00,
                DG_min=0.0, DG_max=8.0, M=30,
                Omega_0=1.0, T_lim = True, N_T0 = 1e6
    )
    T = 15
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    fig_Z, ax_Z = plt.subplots(figsize=(8, 5))
    fig_Z_memory, ax_Z_memory = plt.subplots(figsize=(8, 5))
    axins = ax_Z_memory.inset_axes([0.13, 0.65, 0.3, 0.3])
    fig_N, ax_N = plt.subplots(figsize=(8, 5))
    colors_sim = [my_red, my_blue2, my_purple2, my_gold, my_brown, my_blue, my_green, 'tab:orange', my_purple, my_cyan]
    N_ensembles = 1

    h0s = np.logspace(-4, np.log10(base['b0']), 10)
    pi_stars = (base['b0']/h0s)**(1/base['hill'])
    initial_memory_potency = []
    final_primary_potency = []
    print(f"Running simulation")
    for i_h0, h0 in enumerate(h0s):
        print(f"... for h0={h0:.2g}")
        base['h0'] = h0
        pi_star = (base['b0']/base['h0'])**(1/base['hill'])
        label = f'${pi_star:.2g}$'
        
        for i_m, memory in enumerate([0]):
            print(f"... ... for memory={memory}")
            base['memory'] = memory
            figs, axes = plt.subplots(3, 2, figsize=(16, 15))

            p = Parameters(**base)
            
            initial_memory_potencies_pi_star = []
            final_primary_potencies_pi_star = []

            for i_ensemble in range(N_ensembles):
                # res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='stochastic', seed=None)
                res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='grid')
                # print(res['M'])
                t = res['t']
                N_B = res['N_Bo'] + res['N_Ba']
                # N_B[N_B<2.0] = 0
                N_T = res['N_To'] + res['N_Ta']

                final_primary_potencies_pi_star.append(compute_potency(res))#/compute_yield(res))
                DG_memory, N_memory = produce_memory(res)
                Z0_memory = np.sum(N_memory*np.exp(-DG_memory))#/np.sum(N_memory))
                initial_memory_potencies_pi_star.append(Z0_memory)
                

            ax_Z_memory.plot(DG_memory, N_memory, color=colors_sim[i_h0], alpha=0.5, label=f'${h0:.2g}$')
            initial_memory_potency.append(np.mean(initial_memory_potencies_pi_star))
            final_primary_potency.append(np.mean(final_primary_potencies_pi_star))

            # (a) N_A
            axi = axes[0, 0]
            axi.plot(t, res['N_A'], label=label, linewidth = 2, color = antigen_color)

            # (b) pi
            axi = axes[1, 0]
            ax1 = axi.plot(t, res['pi'][0, :], label=label, linewidth = 2)
            axi.plot(t, res['pi'][10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            axi.plot(t, res['pi'][-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            axi.plot(t, 1e8*np.exp(-(p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            axi.axhline(pi_star, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # (c) lambda_B
            axi = axes[2, 0]
            ax1 = axi.plot(t, p.b0*res['N_Ba'][0, :]/(N_B[0, :])/p.b0, label=label, linewidth = 2)
            axi.plot(t, p.b0*res['N_Ba'][10, :]/(N_B[10, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.5)
            axi.plot(t, p.b0*res['N_Ba'][-10, :]/(N_B[-10, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.2)
            # axi.plot(t, ((p.h0 * (res['pi'][0, :]**p.hill / (res['pi'][0, :]**p.hill + p.pi_star**p.hill)) * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
            axi.plot(t, ((p.h0 * res['pi'][0, :]**p.hill * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
            axi.axhline(p.b0/p.b0, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
    
            # (d) T
            axi = axes[0, 1]
            ax1 = axi.semilogy(t, N_T, label=label, linewidth = 3)
            axi.semilogy(t, res['N_To'], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())

            # (e) B
            axi = axes[1, 1]
            # activated = N_B > 1.5
            N_B_total = compute_N_B_tot(res)
            ax1 = axi.semilogy(t, N_B_total, label=label, linewidth = 4)
            axi.semilogy(t, N_B[0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
            axi.semilogy(t, N_B[10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            axi.semilogy(t, N_B[-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            # axi.semilogy(t, res['N_Bo'][0, :], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
            # axi.semilogy(t, res['N_Ba'][0, :], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())
            # axi.semilogy(t, N_B_total, label=label, linewidth = 3, linestyle='-', alpha=0.8, color=ax1[0].get_color())
            axi.semilogy(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            L_act = compute_L_act(res)
            ax1 = axi.semilogy(t, L_act, label=label, linewidth = 4, color = my_green)

            # (f) Potency
            axi = axes[2, 1]
            Z_B_total = compute_potency_t(res)
            ax1 = axi.semilogy(t, Z_B_total, label=label, linewidth = 4)
            Z_B = N_B * res['weights'][:, None] * np.exp(-res['DG'][:, None])
            Z_B[N_B <= 2.0] = np.nan
            axi.plot(t, Z_B[0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
            axi.plot(t, Z_B[10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            axi.plot(t, Z_B[-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            axi.plot(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            
            # Clone size shared axis
            ax1 = ax_N.semilogy(t, N_B_total, linewidth = 3, color=colors_sim[i_h0], label=label)
            # ax_N.semilogy(t, res['N_Bo'][0, :] + res['N_Ba'][0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
            # ax_N.semilogy(t, res['N_Bo'][10, :] + res['N_Ba'][10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            # ax_N.semilogy(t, res['N_Bo'][-1, :] + res['N_Ba'][-1, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())

            # Potency shared axis
            if memory:
                color = my_blue
            else:
                color = my_red
            ax1 = ax_Z.plot(t, Z_B_total, linewidth = 2, color=colors_sim[i_h0], label=label)
            # ax_Z.semilogy(t, Z_B[0, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())

            # Formatting
            # axes[0].set_xlabel('Time')
            axes[0, 0].set_ylabel('$N_A$', fontsize = 16)
            axes[0, 0].set_xticklabels([])
            axes[0, 0].set_ylim(bottom = 1e-1)
            axes[0, 0].set_xlim(0, T)
            axes[0, 0].tick_params(labelsize=14)

            # axes[1].set_xlabel('Time')
            axes[1, 0].set_ylabel(r'$\pi$', fontsize = 16)
            axes[1, 0].set_xticklabels([])
            axes[1, 0].set_ylim(bottom = 1e-2, top = 1e4)
            axes[1, 0].set_xlim(0, T)
            axes[1, 0].set_yscale('log')

            axes[2, 0].set_xlabel('Time', fontsize = 16)
            axes[2, 0].set_ylabel(r'$\lambda_B/b_o$', fontsize = 16)
            # axes[2, 0].set_xticklabels([])
            axes[2, 0].set_xlim(0, T)
            axes[2, 0].tick_params(labelsize=14)

            # axes[3].set_xlabel('Time')
            axes[0, 1].set_ylabel('T cells', fontsize = 16)
            axes[0, 1].axhline(1.0, color='k', linestyle=':', alpha=0.5)
            axes[0, 1].set_xticklabels([])
            axes[0, 1].set_xlim(0, T)
            axes[0, 1].set_yscale('log')
            axes[0, 1].set_ylim(bottom = 0.1*p.N_T0, top = 1e5*p.N_T0)
            # axes[0, 1].legend(fontsize=10)
            axes[0, 1].tick_params(labelsize=14)

            axes[1, 1].axhline(1.0, color='k', linestyle='--', alpha=0.5)
            # axes[1, 1].set_xlabel('Time')
            axes[1, 1].set_ylabel('B cells', fontsize = 16)
            axes[1, 1].set_xticklabels([])
            # axes[1, 1].set_xlabel('Time', fontsize = 16)
            axes[1, 1].set_ylim(bottom = 5e-1, top = 1e10)
            axes[1, 1].set_xlim(0, T)
            axes[1, 1].set_yscale('log')
            axes[1, 1].tick_params(labelsize=14)

            # axes[2, 0].set_xlabel('Time')
            axes[2, 1].set_ylabel('Potency, $Z$', fontsize = 16)
            # axes[2, 1].set_xticklabels([])
            axes[2, 1].set_xlabel('Time', fontsize = 16)
            axes[2, 1].set_ylim(bottom = 5e-1, top = 1e8)
            axes[2, 1].set_xlim(0, T)
            axes[2, 1].set_yscale('log')
            axes[2, 1].tick_params(labelsize=14)
            
            plt.tight_layout()
            figs.savefig(os.path.join(output_plot, f'summary_pi_star-{pi_star:.2g}_Mem-{int(p.memory)}.pdf'), dpi=150, bbox_inches='tight')


    ax_N.plot(t, 1e-2*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
    ax_N.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    # ax_N.set_xlabel('Time')
    ax_N.set_ylabel(r'Yield, $\bar{N}$', fontsize = 16)
    # ax_N.set_xticklabels([])
    ax_N.set_xlabel('Time', fontsize = 16)
    ax_N.set_ylim(bottom = 5e-1, top = 1e10)
    ax_N.set_xlim(0, T)
    ax_N.set_yscale('log')
    ax_N.tick_params(labelsize=14)
    ax_N.legend(fontsize=14)
    plt.tight_layout()
    fig_N.savefig(os.path.join(output_plot, f'yield.pdf'), dpi=150, bbox_inches='tight')

    ax_Z.plot(t, 1e-2*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
    ax_Z.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    # ax_Z.set_xlabel('Time')
    ax_Z.set_ylabel('Potency, $Z$', fontsize = 16)
    # ax_Z.set_xticklabels([])
    ax_Z.set_xlabel('Time', fontsize = 16)
    ax_Z.set_ylim(bottom = 5e-1, top = 1e7)
    ax_Z.set_xlim(0, T)
    ax_Z.set_yscale('log')
    ax_Z.tick_params(labelsize=14)
    ax_Z.legend(fontsize=14)
    plt.tight_layout()
    fig_Z.savefig(os.path.join(output_plot, f'potency.pdf'), dpi=150, bbox_inches='tight')


    axins.plot(pi_stars, initial_memory_potency, marker='o', color='k', alpha = .7)
    # axins.plot(pi_stars, final_primary_potency, marker='^', color='grey', alpha = .7)
    axins.set_xscale('log')
    axins.set_yscale('log')
    axins.set_xlabel(r'$\pi^*$', fontsize=10)
    axins.set_ylabel(r'$\langle z \rangle$', fontsize = 10)
    axins.tick_params(which='both', labelsize=6)
    ax_Z_memory.plot(res['DG'], res['weights'], marker='o', color='grey')
    # ax_Z_memory.set_ylabel(r'$\mathrm{Potency, } Z_0$', fontsize = 16)
    # ax_Z_memory.set_xticklabels([])
    # ax_Z_memory.set_xlabel(r'$h_0$', fontsize = 16)
    # ax_Z_memory.set_ylim(bottom = 5e-1, top = 1e7)
    # ax_Z_memory.set_xlim(0, T)
    # ax_Z_memory.set_xscale('log')
    ax_Z_memory.set_yscale('log') 
    ax_Z_memory.tick_params(labelsize=14)
    ax_Z_memory.legend(title = r'$\pi^*$', title_fontsize = 12, fontsize=10, loc = 4)
    plt.tight_layout()
    fig_Z_memory.savefig(os.path.join(output_plot, f'Z0_memory.pdf'), dpi=150, bbox_inches='tight')
