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
subsubproject = 'h0'

if __name__ == '__main__':
    output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}/{subproject}/{subsubproject}/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, lambda_A = 6.0, delta_A=6.0,
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
    colors_sim = [my_red, my_blue2, my_purple2, my_gold, my_brown, my_blue, my_green, 'tab:orange', my_purple, my_cyan]
    colors_sim = plt.cm.viridis(np.linspace(0, 1, 50))
    N_ensembles = 1
    fig_Z_memory, ax_Z_memory = plt.subplots(figsize=(8, 5))
    fig_Y_memory, ax_Y_memory = plt.subplots(figsize=(8, 5))

    print(f"Running simulation")
    etas = [1.0, 1.5, 2.0, 3.0]
    for eta in etas:
        print(f"... for eta={eta:.1g}")
        fig_Z_t, ax_Z_t = plt.subplots(figsize=(8, 5))
        fig_N_t, ax_N_t = plt.subplots(figsize=(8, 5))
        fig_hist, ax_hist = plt.subplots(figsize=(8, 5))
        fig_Z_Y_total, ax_Z_Y_total = plt.subplots(figsize=(8, 5))
        fig_Z_Y_mean, ax_Z_Y_mean = plt.subplots(figsize=(8, 5))
        
        h0s = np.logspace(np.log10(base['b0']/4e3), np.log10(base['b0']), 50)
        pi_stars = (base['b0']/h0s)**(1/base['hill'])
        initial_memory_potency = []
        final_primary_potency = []
        initial_memory_yield = []
        final_primary_yield = []
        base['eta'] = eta
        for i_h0, h0 in enumerate(h0s):
            # print(f"... for h0={h0:.2g}")
            base['h0'] = h0
            pi_star = (base['b0']/base['h0'])**(1/base['hill'])
            label = f'${pi_star:.2g}$'
            
            for i_m, memory in enumerate([0]):
                # print(f"... ... for memory={memory}")
                base['memory'] = memory
                # figs, axes = plt.subplots(3, 2, figsize=(16, 15))

                p = Parameters(**base)
                
                initial_memory_potencies_pi_star = []
                final_primary_potencies_pi_star = []
                initial_memory_yield_pi_star = []
                final_primary_yield_pi_star = []

                for i_ensemble in range(N_ensembles):
                    # res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='stochastic', seed=None)
                    res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='grid')
                    # print(res['M'])
                    t = res['t']
                    N_B = res['N_Bo'] + res['N_Ba']
                    # N_B[N_B<2.0] = 0
                    N_T = res['N_To'] + res['N_Ta']

                    final_primary_potencies_pi_star.append(compute_potency(res))
                    final_primary_yield_pi_star.append(compute_yield(res))

                    DG_memory, N_memory = produce_memory(res)
                    Z0_memory = np.sum(N_memory*np.exp(-DG_memory))
                    initial_memory_potencies_pi_star.append(Z0_memory)
                    initial_memory_yield_pi_star.append(np.sum(N_memory))
                    

                ax_hist.plot(DG_memory, N_memory, color=colors_sim[i_h0], alpha=0.5, label=label)
                initial_memory_potency.append(np.mean(initial_memory_potencies_pi_star))
                final_primary_potency.append(np.mean(final_primary_potencies_pi_star))
                initial_memory_yield.append(np.mean(initial_memory_yield_pi_star))
                final_primary_yield.append(np.mean(final_primary_yield_pi_star))

                # # (a) N_A
                # axi = axes[0, 0]
                # axi.plot(t, res['N_A'], label=label, linewidth = 2, color = antigen_color)

                # # (b) pi
                # axi = axes[1, 0]
                # ax1 = axi.plot(t, res['pi'][0, :], label=label, linewidth = 2)
                # axi.plot(t, res['pi'][10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
                # axi.plot(t, res['pi'][-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
                # axi.plot(t, 1e8*np.exp(-(p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
                # axi.axhline(pi_star, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

                # # (c) lambda_B
                # axi = axes[2, 0]
                # ax1 = axi.plot(t, p.b0*res['N_Ba'][0, :]/(N_B[0, :])/p.b0, label=label, linewidth = 2)
                # axi.plot(t, p.b0*res['N_Ba'][10, :]/(N_B[10, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.5)
                # axi.plot(t, p.b0*res['N_Ba'][-10, :]/(N_B[-10, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.2)
                # # axi.plot(t, ((p.h0 * (res['pi'][0, :]**p.hill / (res['pi'][0, :]**p.hill + p.pi_star**p.hill)) * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
                # axi.plot(t, ((p.h0 * res['pi'][0, :]**p.hill * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
                # axi.axhline(p.b0/p.b0, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
        
                # # (d) T
                # axi = axes[0, 1]
                # ax1 = axi.semilogy(t, N_T, label=label, linewidth = 3)
                # axi.semilogy(t, res['N_To'], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())

                # # (e) B
                # axi = axes[1, 1]
                # # activated = N_B > 1.5
                N_B_total = compute_N_B_tot(res)
                # ax1 = axi.semilogy(t, N_B_total, label=label, linewidth = 4)
                # axi.semilogy(t, N_B[0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
                # axi.semilogy(t, N_B[10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
                # axi.semilogy(t, N_B[-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
                # # axi.semilogy(t, res['N_Bo'][0, :], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
                # # axi.semilogy(t, res['N_Ba'][0, :], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())
                # # axi.semilogy(t, N_B_total, label=label, linewidth = 3, linestyle='-', alpha=0.8, color=ax1[0].get_color())
                # axi.semilogy(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
                # L_act = compute_L_act(res)
                # ax1 = axi.semilogy(t, L_act, label=label, linewidth = 4, color = my_green)

                # # (f) Potency
                # axi = axes[2, 1]
                Z_B_total = compute_potency_t(res)
                # ax1 = axi.semilogy(t, Z_B_total, label=label, linewidth = 4)
                # Z_B = N_B * res['weights'][:, None] * np.exp(-res['DG'][:, None])
                # Z_B[N_B <= 2.0] = np.nan
                # axi.plot(t, Z_B[0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
                # axi.plot(t, Z_B[10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
                # axi.plot(t, Z_B[-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
                # axi.plot(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
                
                # Clone size shared axis
                ax1 = ax_N_t.semilogy(t, N_B_total, linewidth = 3, color=colors_sim[i_h0], label=label)
                # ax_N_t.semilogy(t, res['N_Bo'][0, :] + res['N_Ba'][0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
                # ax_N_t.semilogy(t, res['N_Bo'][10, :] + res['N_Ba'][10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
                # ax_N_t.semilogy(t, res['N_Bo'][-1, :] + res['N_Ba'][-1, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())

                # Potency shared axis
                if memory:
                    color = my_blue
                else:
                    color = my_red
                ax1 = ax_Z_t.plot(t, Z_B_total, linewidth = 2, color=colors_sim[i_h0], label=label)
                # ax_Z_t.semilogy(t, Z_B[0, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())

                # Formatting
                # axes[0].set_xlabel('Time')
                # axes[0, 0].set_ylabel('$N_A$', fontsize = 16)
                # axes[0, 0].set_xticklabels([])
                # axes[0, 0].set_ylim(bottom = 1e-1)
                # axes[0, 0].set_xlim(0, T)
                # axes[0, 0].tick_params(labelsize=14)

                # # axes[1].set_xlabel('Time')
                # axes[1, 0].set_ylabel(r'$\pi$', fontsize = 16)
                # axes[1, 0].set_xticklabels([])
                # axes[1, 0].set_ylim(bottom = 1e-2, top = 1e4)
                # axes[1, 0].set_xlim(0, T)
                # axes[1, 0].set_yscale('log')

                # axes[2, 0].set_xlabel('Time', fontsize = 16)
                # axes[2, 0].set_ylabel(r'$\lambda_B/b_o$', fontsize = 16)
                # # axes[2, 0].set_xticklabels([])
                # axes[2, 0].set_xlim(0, T)
                # axes[2, 0].tick_params(labelsize=14)

                # # axes[3].set_xlabel('Time')
                # axes[0, 1].set_ylabel('T cells', fontsize = 16)
                # axes[0, 1].axhline(1.0, color='k', linestyle=':', alpha=0.5)
                # axes[0, 1].set_xticklabels([])
                # axes[0, 1].set_xlim(0, T)
                # axes[0, 1].set_yscale('log')
                # axes[0, 1].set_ylim(bottom = 0.1*p.N_T0, top = 1e5*p.N_T0)
                # # axes[0, 1].legend(fontsize=10)
                # axes[0, 1].tick_params(labelsize=14)

                # axes[1, 1].axhline(1.0, color='k', linestyle='--', alpha=0.5)
                # # axes[1, 1].set_xlabel('Time')
                # axes[1, 1].set_ylabel('B cells', fontsize = 16)
                # axes[1, 1].set_xticklabels([])
                # # axes[1, 1].set_xlabel('Time', fontsize = 16)
                # axes[1, 1].set_ylim(bottom = 5e-1, top = 1e10)
                # axes[1, 1].set_xlim(0, T)
                # axes[1, 1].set_yscale('log')
                # axes[1, 1].tick_params(labelsize=14)

                # # axes[2, 0].set_xlabel('Time')
                # axes[2, 1].set_ylabel('Potency, $Z$', fontsize = 16)
                # # axes[2, 1].set_xticklabels([])
                # axes[2, 1].set_xlabel('Time', fontsize = 16)
                # axes[2, 1].set_ylim(bottom = 5e-1, top = 1e8)
                # axes[2, 1].set_xlim(0, T)
                # axes[2, 1].set_yscale('log')
                # axes[2, 1].tick_params(labelsize=14)
                
                # plt.tight_layout()
                # figs.savefig(os.path.join(output_plot, f'summary_pi_star-{pi_star:.2g}_Mem-{int(p.memory)}.pdf'), dpi=150, bbox_inches='tight')

        ax_N_t.plot(t, 1e-2*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
        ax_N_t.axhline(1.0, color='k', linestyle='--', alpha=0.5)
        # ax_N_t.set_xlabel('Time')
        ax_N_t.set_ylabel(r'Yield, $\bar{N}$', fontsize = 16)
        # ax_N_t.set_xticklabels([])
        ax_N_t.set_xlabel('Time', fontsize = 16)
        ax_N_t.set_ylim(bottom = 5e-1, top = 1e10)
        ax_N_t.set_xlim(0, T)
        ax_N_t.set_yscale('log')
        ax_N_t.tick_params(labelsize=14)
        # ax_N_t.legend(fontsize=14)
        fig_N_t.savefig(os.path.join(output_plot, f'yield_eta-{eta:.1f}.pdf'), dpi=150, bbox_inches='tight')
        plt.close(fig_N_t)
        
        ax_Z_t.plot(t, 1e-2*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
        ax_Z_t.axhline(1.0, color='k', linestyle='--', alpha=0.5)
        # ax_Z_t.set_xlabel('Time')
        ax_Z_t.set_ylabel('Potency, $Z$', fontsize = 16)
        # ax_Z_t.set_xticklabels([])
        ax_Z_t.set_xlabel('Time', fontsize = 16)
        ax_Z_t.set_ylim(bottom = 5e-1, top = 1e7)
        ax_Z_t.set_xlim(0, T)
        ax_Z_t.set_yscale('log')
        ax_Z_t.tick_params(labelsize=14)
        # ax_Z_t.legend(fontsize=14)
        fig_Z_t.savefig(os.path.join(output_plot, f'potency_eta-{eta:.1f}.pdf'), dpi=150, bbox_inches='tight')
        plt.close(fig_Z_t)

        ax_hist.plot(res['DG'], res['weights'], marker='o', color='grey')
        ax_hist.set_ylabel(r'$\mathrm{\# cells}$', fontsize = 16)
        ax_hist.set_xlabel(r'$\Delta G$', fontsize = 16)
        # ax_hist.set_xticklabels([])
        # ax_hist.set_ylim(bottom = 5e-1, top = 1e7)
        ax_hist.set_ylim(bottom = 8e-1, top = 1e6)
        # ax_hist.set_xscale('log')
        ax_hist.set_yscale('log') 
        ax_hist.tick_params(labelsize=14)
        # ax_hist.legend(title = r'$\pi^*$', title_fontsize = 12, fontsize=10, loc = 4)
        fig_hist.savefig(os.path.join(output_plot, f'Z0_memory_eta-{eta:.1f}.pdf'), dpi=150, bbox_inches='tight')
        plt.close(fig_hist)

        ax_Z_Y_total.plot(pi_stars, initial_memory_potency, marker='o', color=my_purple, alpha = .7, ms = 8, lw = 2, label=r'$\mathrm{Memory\ potency}$')
        ax_Z_Y_total.plot(pi_stars, initial_memory_yield, marker='^', color=my_purple, alpha = .7, ms = 8, lw = 2, label=r'$\mathrm{Memory\ yield}$')
        ax_Z_Y_total.plot(pi_stars, final_primary_potency, marker='o', color=my_red, alpha = .7, ms = 8, lw = 2, label=r'$\mathrm{Primary\ potency}$')
        ax_Z_Y_total.plot(pi_stars, final_primary_yield, marker='^', color=my_red, alpha = .7, ms = 8, lw = 2, label=r'$\mathrm{Primary\ yield}$')
        ax_Z_Y_total.set_xscale('log')
        ax_Z_Y_total.set_yscale('log')
        ax_Z_Y_total.set_xlabel(r'$\pi^*$', fontsize=14)
        ax_Z_Y_total.set_ylabel(r'$Z\, \mathrm{\ and \ }\, Y$ ', fontsize = 14)
        ax_Z_Y_total.set_ylim(bottom = 2e-1, top = 5e9)
        ax_Z_Y_total.tick_params(which='both', labelsize=14)
        ax_Z_Y_total.legend(fontsize=14)
        fig_Z_Y_total.savefig(os.path.join(output_plot, f'Z_memory_total_eta-{eta:.1f}.pdf'), dpi=150, bbox_inches='tight')
        plt.close(fig_Z_Y_total)

        ax_Z_Y_mean.plot(pi_stars, np.array(initial_memory_potency)/np.array(initial_memory_yield), marker='D', color=my_purple, alpha = .7, ms = 8, lw = 2, label=r'$\mathrm{Memory\ potency}$')
        ax_Z_Y_mean.plot(pi_stars, np.array(final_primary_potency)/np.array(final_primary_yield), marker='D', color=my_red, alpha = .7, ms = 8, lw = 2, label=r'$\mathrm{Primary\ potency}$')
        ax_Z_Y_mean.set_xscale('log')
        ax_Z_Y_mean.set_yscale('log')
        ax_Z_Y_mean.set_xlabel(r'$\pi^*$', fontsize=14)
        ax_Z_Y_mean.set_ylabel(r'$\langle z \rangle$', fontsize = 14)
        ax_Z_Y_mean.set_ylim(bottom = 5e-4, top = 1e1)
        ax_Z_Y_mean.tick_params(which='both', labelsize=14)
        ax_Z_Y_mean.legend(fontsize=14)
        fig_Z_Y_mean.savefig(os.path.join(output_plot, f'Z_memory_mean_eta-{eta:.1f}.pdf'), dpi=150, bbox_inches='tight')
        plt.close(fig_Z_Y_mean)

        ax_Z_memory.plot(pi_stars, initial_memory_potency, marker='', alpha = .7, ms = 10, lw = 2, label = f'${eta:.1f}$')
        ax_Y_memory.plot(pi_stars, initial_memory_yield, marker='', alpha = .7, ms = 10, lw = 2, label = f'${eta:.1f}$')
        

    ax_Z_memory.set_xscale('log')
    ax_Z_memory.set_yscale('log')
    ax_Z_memory.set_xlabel(r'$\pi^*$', fontsize=14)
    ax_Z_memory.set_ylabel(r'$Z$', fontsize = 14)
    ax_Z_memory.set_ylim(bottom = 5e-1, top = 5e3)
    ax_Z_memory.tick_params(which='both', labelsize=14)
    ax_Z_memory.legend(fontsize=12, title = r'$\eta$', title_fontsize = 14)
    fig_Z_memory.savefig(os.path.join(output_plot, f'Z_memory.pdf'), dpi=150, bbox_inches='tight')

    ax_Y_memory.set_xscale('log')
    ax_Y_memory.set_yscale('log')
    ax_Y_memory.set_xlabel(r'$\pi^*$', fontsize=14)
    ax_Y_memory.set_ylabel(r'$Y$', fontsize = 14)
    ax_Y_memory.set_ylim(bottom = 5e-1, top = 2e4)
    ax_Y_memory.tick_params(which='both', labelsize=14)
    ax_Y_memory.legend(fontsize=12, title = r'$\eta$', title_fontsize = 14)
    fig_Y_memory.savefig(os.path.join(output_plot, f'Y_memory.pdf'), dpi=150, bbox_inches='tight')