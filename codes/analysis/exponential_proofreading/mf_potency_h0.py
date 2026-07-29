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

fig_kw = dict(figsize=(8 * 1.62, 8), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
fig_kw2 = dict(figsize=(8 * 1.62, 5), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
if __name__ == '__main__':
    # Default parameters
    base = dict(N_A0=1.0, lambda_A = 6.0, delta_A=3.0, eta= 1.0,
                k_on=1e1*2e5*1e6*24*3600/N_Avg, delta_pi=0.1,
                hill=1.0, beta_star=2.5, K_T = 1e4,
                delta_T=0.00, Tcell_growth_factor=2.0,
                tau_eng=0.1, b0=2.0, delta_B=0.00,
                DG_min=0.0, DG_max=8.0, M=30,
                omega_0=1.0, T_lim = True, N_T0 = 1e6
    )
    T = 12
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    # fig_Z_memory_hist, ax_Z_memory_hist = plt.subplots(**fig_kw)
    # fig_Z_memory_total, ax_Z_memory_total = plt.subplots(**fig_kw)
    # fig_Z_memory_mean, ax_Z_memory_mean = plt.subplots(**fig_kw)
    
    colors_mem = [my_green, my_blue2]
    styles_sim = ['--', '-', ':', '-.', '-', '--', ':', '-.', '-', '--']
    N_ensemble = 1
    N_h0 = 5
    colors_h0 = plt.cm.brg(np.linspace(0, 1, N_h0))
    widths_h0 = [4, 3, 2, 1, 4, 3, 2, 1, 4, 3]
    h0s = np.flip(np.logspace(np.log10(base['b0']/1e4), np.log10(base['b0']), N_h0))
    pi_stars = (base['b0']/h0s)**(1/base['hill'])
    initial_memory_potency = []
    final_primary_potency = []
    initial_memory_yield = []
    final_primary_yield = []


    fig_NA_shared, ax_NA_shared = plt.subplots(**fig_kw)
    fig_DG_shared, ax_DG_shared = plt.subplots(**fig_kw)
    fig_Z_shared, ax_Z_shared = plt.subplots(**fig_kw)
    fig_n0_shared, ax_n0_shared = plt.subplots(**fig_kw)
    print(f"Running simulation")
    for i_h0, h0 in enumerate(h0s):
        print(f"... for h0={h0:.2g}")
        base['h0'] = h0
        pi_star = (base['b0']/base['h0'])**(1/base['hill'])
        label = f'${pi_star:.2g}$'
        

        for i_m, memory in enumerate([0, 1]):
            output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}/{subproject}/{subsubproject}/memory_{memory}/'
            os.makedirs(output_plot, exist_ok=True)
            print(f"... ... for memory={memory}")
            base['memory'] = memory

            p = Parameters(**base)
            
            initial_memory_potencies_pi_star = []
            final_primary_potencies_pi_star = []
            initial_memory_yield_pi_star = []
            final_primary_yield_pi_star = []

            for i_ensemble in range(N_ensemble):
                # res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='stochastic', seed=None)
                if p.memory == 0:
                    res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='grid')
                    expanded_cells, formed_memory = memory_seed_from_primary(res, p=p, n_mem=int(1e4))
                    ax_n0_shared.plot(res['DG'], formed_memory*res['weights'], color=colors_h0[i_h0], alpha=0.8, label=label)
                    ax_n0_shared.plot(res['DG'], expanded_cells, color=colors_h0[i_h0], alpha=0.8, label=label, ls = styles_sim[i_m])
                else:
                    res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='grid', memory_seed=formed_memory)
                # print(res['M'])
                t = res['t']
                N_B = res['N_Bo'] + res['N_Ba']
                # N_B[N_B<2.0] = 0
                N_T = res['N_To'] + res['N_Ta']

                # final_primary_potencies_pi_star.append(compute_potency(res))
                # final_primary_yield_pi_star.append(compute_yield(res))

                DG_memory, N_memory = produce_memory(res)
                # Z0_memory = np.sum(N_memory*np.exp(-DG_memory))
                # initial_memory_potencies_pi_star.append(Z0_memory)
                # initial_memory_yield_pi_star.append(np.sum(N_memory))
                

            # ax_Z_memory_hist.plot(DG_memory, N_memory, color=colors_mem[i_h0], alpha=0.5, label=label)
            # initial_memory_potency.append(np.mean(initial_memory_potencies_pi_star))
            # final_primary_potency.append(np.mean(final_primary_potencies_pi_star))
            # initial_memory_yield.append(np.mean(initial_memory_yield_pi_star))
            # final_primary_yield.append(np.mean(final_primary_yield_pi_star))

            # N_A
            if p.memory == 0:
                ax_NA_shared.plot(t, res['N_A'], label=label, linewidth = 4, color = antigen_color, ls = styles_sim[i_m])
            else:
                ax_NA_shared.plot(t, res['N_A'], label=label, linewidth = 4, color = colors_h0[i_h0], ls = styles_sim[i_m])

            # # pi
            # ax1 = ax_pi.plot(t, res['pi'][0, :], label=label, linewidth = 4, color=colors_mem[i_m])
            # ax_pi.plot(t, res['pi'][10, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())
            # ax_pi.plot(t, res['pi'][-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            # # ax_pi.plot(t, 1e8*np.exp(-(p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            # ax_pi.axhline(pi_star, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # lambda_B
            # ax1 = ax_lamB.plot(t, p.b0*res['N_Ba'][0, :]/(N_B[0, :])/p.b0, label=label, linewidth = 2, color = colors_mem[i_m])
            # ax_lamB.plot(t, p.b0*res['N_Ba'][10, :]/(N_B[10, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.5)
            # ax_lamB.plot(t, p.b0*res['N_Ba'][-10, :]/(N_B[-10, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.2)
            # # ax_lamB.plot(t, ((p.h0 * (res['pi'][0, :]**p.hill / (res['pi'][0, :]**p.hill + p.pi_star**p.hill)) * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
            # ax_lamB.plot(t, ((p.h0 * res['pi'][0, :]**p.hill * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
            # ax_lamB.axhline(p.b0/p.b0, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
    
            # T
            # ax1 = ax_NT.semilogy(t, N_T, label=label, linewidth = 3)
            # ax_NT.semilogy(t, res['N_To'], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())

            # B
            # activated = N_B > 1.5
            N_B_total = compute_N_B_tot(res)
            L_act = compute_L_act(res)
            # ax1 = ax_NB.plot(t, N_B_total, label=label, linewidth = 4, color = ax1[0].get_color())
            # ax1 = ax_NB.plot(t, N_B[0, :], label=label, linewidth = 4, alpha=1, color=colors_mem[i_m])
            # ax_NB.plot(t, N_B[10, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())
            # ax_NB.plot(t, N_B[-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            # ax_NB.plot(t, res['N_Bo'][0, :], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
            # ax_NB.plot(t, res['N_Ba'][0, :], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())
            # ax_NB.plot(t, N_B_total, label=label, linewidth = 3, linestyle='-', alpha=0.8, color=ax1[0].get_color())
            # ax_NB.plot(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # Y
            # ax1 = ax_Y.plot(t, N_B_total, linewidth = 4, color=colors_mem[i_m], label=label)
            # ax_Y.semilogy(t, res['N_Bo'][0, :] + res['N_Ba'][0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
            # ax_Y.semilogy(t, res['N_Bo'][10, :] + res['N_Ba'][10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            # ax_Y.semilogy(t, res['N_Bo'][-1, :] + res['N_Ba'][-1, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())

            # Z
            Z_B_total = compute_potency_t(res)
            Z_B = N_B * res['weights'][:, None] * np.exp(-res['DG'][:, None])
            Z_B[N_B <= 2.0] = np.nan

            if p.memory == 0:
                ax_Z_shared.plot(t[Z_B_total>0], Z_B_total[Z_B_total>0], color=colors_h0[i_h0], ls = styles_sim[i_m], lw = 3)
                # lambda_prime = p.lambda_A*p.eta*p.beta_star*(1-0.5)
                # ax_Z_shared.plot(t[t<4.5], 1e-12*np.exp((p.b0 + lambda_prime)*t[t<4.5]), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
                # ax_Z_shared.plot(t[t<8], 2e-2*np.exp((p.b0)*t[t<8]), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            else:
                ax_Z_shared.plot(t[Z_B_total>0], Z_B_total[Z_B_total>0], color=colors_h0[i_h0], label=label, ls = styles_sim[i_m], lw = 2)
                # ax_Z_shared.plot(t, 5e1*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            
            if p.memory ==0:
                tstar = t[N_B[0, :]>N_B[0, 0]+0.1*N_B[0, 0]][0]
                tpeak = t[res['N_A']==np.max(res['N_A'])][0]
                ton = t[(t<tpeak) & (t>tstar)]
                toff = t[(t>=tpeak)]
                DGpeak = p.lambda_A/p.eta * (tpeak - tstar)
                DG_on = p.lambda_A/p.eta * (ton - tstar)
                DG_off = p.b0 /p.eta * (tpeak - toff) + DGpeak
                DG_off_null = p.delta_A/p.eta * (tpeak - toff) + DGpeak
                ax_DG_shared.plot(ton, DG_on, color=colors_h0[i_h0], ls = 'dotted', label=label)
                ax_DG_shared.plot(toff[DG_off>p.DG_min], DG_off[DG_off>p.DG_min], color=colors_h0[i_h0], ls = 'dashed', lw = 3)
                # ax_DG_shared.plot(toff[DG_off_null>p.DG_min], DG_off_null[DG_off_null>p.DG_min], linewidth = widths_h0[i_h0], color=my_brown, ls = 'dashed')


        output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}/{subproject}/{subsubproject}/'
        os.makedirs(output_plot, exist_ok=True)

    
    # ax_NA_shared.set_ylabel('$N_A$', fontsize = 16)
    ax_NA_shared.set_xticklabels([])
    ax_NA_shared.set_ylim(bottom = 1e0, top = 1e15)
    ax_NA_shared.set_xlim(0, T)
    ax_NA_shared.set_yscale('log')
    ax_NA_shared.tick_params(axis='y', labelsize=30)
    ax_NA_shared.tick_params(axis='x', labelsize=30)
    fig_NA_shared.savefig(os.path.join(output_plot, f'N_A_shared.pdf'), dpi=150)

    ax_n0_shared.set_ylabel(r'$N_0$', fontsize = 16)
    ax_n0_shared.set_xlabel(r'$\Delta G$', fontsize = 16)
    ax_n0_shared.legend(fontsize=14)
    ax_n0_shared.set_yscale('log')
    fig_n0_shared.savefig(os.path.join(output_plot, f'n0_shared.pdf'), dpi=150)


        # Formatting
    ax_DG_shared.axhline(p.DG_min, tstar, toff[DG_off>p.DG_min][-1], linewidth = 1, linestyle='--', color='grey', alpha=0.8)
    # ax_DG_shared.set_xlabel('Time', fontsize = 16)
    # ax_DG_shared.set_ylabel('B cells', fontsize = 16)
    # ax_DG_shared.set_xticklabels([])
    ax_DG_shared.set_ylim(bottom = p.DG_min-0.1)
    ax_DG_shared.set_xlim(0, T)
    # ax_DG_shared.set_yscale('log')
    ax_DG_shared.tick_params(axis='y', labelsize=30)
    ax_DG_shared.tick_params(axis='x', labelsize=30)
    ax_DG_shared.legend(fontsize=14)
    fig_DG_shared.savefig(os.path.join(output_plot, f'DG_shared.pdf'), dpi=150)
    
    # ax_Z_shared.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    # ax_Z_shared.set_xlabel('Time')
    # ax_Z_shared.set_ylabel('Potency, $Z$', fontsize = 16)
    # ax_Z_shared.set_xticklabels([])
    # ax_Z_shared.set_xlabel('Time', fontsize = 16)
    ax_Z_shared.set_ylim(bottom = 5e-1, top = 2e10)
    ax_Z_shared.set_xlim(0, T)
    ax_Z_shared.set_yscale('log')
    ax_Z_shared.tick_params(axis='y', labelsize=30)
    ax_Z_shared.tick_params(axis='x', labelsize=30)
    ax_Z_shared.legend(fontsize=14)
    fig_Z_shared.savefig(os.path.join(output_plot, f'Z_shared.pdf'), dpi=150)

    # ax_Z_memory_hist.plot(res['DG'], res['weights'], marker='o', color='grey')
    # # ax_Z_memory_hist.set_ylabel(r'$\mathrm{Potency, } Z_0$', fontsize = 16)
    # # ax_Z_memory_hist.set_xticklabels([])
    # # ax_Z_memory_hist.set_xlabel(r'$h_0$', fontsize = 16)
    # # ax_Z_memory_hist.set_ylim(bottom = 5e-1, top = 1e7)
    # # ax_Z_memory_hist.set_xlim(0, T)
    # # ax_Z_memory_hist.set_xscale('log')
    # ax_Z_memory_hist.set_yscale('log') 
    # ax_Z_memory_hist.tick_params(labelsize=14)
    # ax_Z_memory_hist.legend(title = r'$\pi^*$', title_fontsize = 12, fontsize=10, loc = 4)
    # fig_Z_memory_hist.savefig(os.path.join(output_plot, f'Z0_memory.pdf'), dpi=150)

    # ax_Z_memory_total.plot(pi_stars, initial_memory_potency, marker='o', color=my_purple, alpha = .7, ms = 10, lw = 2, label=r'$\mathrm{Memory\ potency}$')
    # ax_Z_memory_total.plot(pi_stars, initial_memory_yield, marker='^', color=my_purple, alpha = .7, ms = 10, lw = 2, label=r'$\mathrm{Memory\ yield}$')
    # ax_Z_memory_total.plot(pi_stars, final_primary_potency, marker='o', color=my_green, alpha = .7, ms = 10, lw = 2, label=r'$\mathrm{Primary\ potency}$')
    # ax_Z_memory_total.plot(pi_stars, final_primary_yield, marker='^', color=my_green, alpha = .7, ms = 10, lw = 2, label=r'$\mathrm{Primary\ yield}$')
    # ax_Z_memory_total.set_xscale('log')
    # ax_Z_memory_total.set_yscale('log')
    # ax_Z_memory_total.set_xlabel(r'$\pi^*$', fontsize=14)
    # ax_Z_memory_total.set_ylabel(r'$Z\, \mathrm{\ and \ }\, Y$ ', fontsize = 14)
    # ax_Z_memory_total.tick_params(which='both', labelsize=14)
    # ax_Z_memory_total.legend(fontsize=14)
    # fig_Z_memory_total.savefig(os.path.join(output_plot, f'Z_memory_total.pdf'), dpi=150)

    # ax_Z_memory_mean.plot(pi_stars, np.array(initial_memory_potency)/np.array(initial_memory_yield), marker='D', color=my_purple, alpha = .7, ms = 10, lw = 2, label=r'$\mathrm{Memory\ potency}$')
    # ax_Z_memory_mean.plot(pi_stars, np.array(final_primary_potency)/np.array(final_primary_yield), marker='D', color=my_green, alpha = .7, ms = 10, lw = 2, label=r'$\mathrm{Primary\ potency}$')
    # ax_Z_memory_mean.set_xscale('log')
    # ax_Z_memory_mean.set_yscale('log')
    # ax_Z_memory_mean.set_xlabel(r'$\pi^*$', fontsize=14)
    # ax_Z_memory_mean.set_ylabel(r'$\langle z \rangle$', fontsize = 14)
    # ax_Z_memory_mean.tick_params(which='both', labelsize=14)
    # ax_Z_memory_mean.legend(fontsize=14)
    # fig_Z_memory_mean.savefig(os.path.join(output_plot, f'Z_memory_mean.pdf'), dpi=150)
