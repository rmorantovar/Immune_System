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
subsubproject = 'dynamics'

fig_kw = dict(figsize=(8 * 1.62, 8), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
fig_kw2 = dict(figsize=(8 * 1.62, 5), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
if __name__ == '__main__':
    # Default parameters
    base = dict(N_A0=1.0, lambda_innate = 2.2, threshold_innate = 5e3,
                lambda_A = 6.2, delta_A=3.0, eta= 1.0,
                k_on=1e2*2e5*1e6*24*3600/N_Avg, delta_pi=0.1,
                hill=2.0, beta_star=2.3, K_T = 1e4,
                delta_T=0.00, Tcell_growth_factor=2.0,
                tau_eng=0.1, b0=1.5, delta_B=0.00,
                DG_min=0.0, DG_max=8.0, M=32,
                omega_0=1.0, T_lim = True, N_T0 = 1e6,
                Z_c = 1e3, n_mem = 1e5
    )
    T = 10

    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    # fig_Z_memory_hist, ax_Z_memory_hist = plt.subplots(**fig_kw)
    # fig_Z_memory_total, ax_Z_memory_total = plt.subplots(**fig_kw)
    # fig_Z_memory_mean, ax_Z_memory_mean = plt.subplots(**fig_kw)
    
    colors_sim = [my_green, my_blue, my_purple2, my_gold, my_brown, my_blue, my_green, 'tab:orange', my_purple, my_cyan]
    styles_sim = ['--', '-', ':', '-.', '-', '--', ':', '-.', '-', '--']
    N_ensemble = 1

    h0s = np.array([base['b0']/1000.])
    pi_stars = (base['b0']/h0s)**(1/base['hill'])
    initial_memory_potency = []
    final_primary_potency = []
    initial_memory_yield = []
    final_primary_yield = []

    print(f"Running simulation")
    for i_h0, h0 in enumerate(h0s):
        print(f"... for h0={h0:.2g}")
        base['h0'] = h0
        fig_NA_shared, ax_NA_shared = plt.subplots(**fig_kw2)
        fig_Z_shared, ax_Z_shared = plt.subplots(**fig_kw2)
        fig_pb_shared, ax_pb_shared = plt.subplots(**fig_kw2)

        for i_m, memory in enumerate([0, 1]):
            output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}/{subproject}/{subsubproject}/memory_{memory}/'
            os.makedirs(output_plot, exist_ok=True)
            print(f"... ... for memory={memory}")
            base['memory'] = memory
            if memory == 1:
                base['h0'] = base['b0']
                alpha = p.eta*(1+p.b0/p.lambda_A)
                print('alpha=', alpha)
            pi_star = (base['b0']/base['h0'])**(1/base['hill'])
            label = f'${pi_star:.2g}$'
            p = Parameters(**base)

            figs, axes = plt.subplots(3, 2, figsize=(16, 15))
            fig_NA, ax_NA = plt.subplots(**fig_kw2)
            fig_pi, ax_pi = plt.subplots(**fig_kw2)
            fig_lamB, ax_lamB = plt.subplots(**fig_kw2)
            fig_NB, ax_NB = plt.subplots(**fig_kw2)
            fig_NT, ax_NT = plt.subplots(**fig_kw2)
            fig_DG, ax_DG = plt.subplots(**fig_kw2)
            fig_Y, ax_Y = plt.subplots(**fig_kw2)
            fig_Z, ax_Z = plt.subplots(**fig_kw2)


            for i_ensemble in range(N_ensemble):
                # res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='stochastic', seed=None)
                if p.memory == 0:
                    res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='grid')
                    expanded_cells, formed_memory = memory_seed_from_primary(res, p=p, n_mem=int(p.n_mem))
                else:
                    res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='grid', memory_seed=formed_memory)
                # print(res['M'])
                t = res['t']
                N_B = res['N_Bo'] + res['N_Ba']
                # N_B[N_B<2.0] = 0
                N_T = res['N_To'] + res['N_Ta']

            # N_A
            ax_NA.plot(t, res['N_A'], label=label, linewidth = 4, color = antigen_color)
            if p.memory == 0:
                ax_NA_shared.plot(t, res['N_A'], label=label, linewidth = 4, color = antigen_color, ls = styles_sim[i_m])
            else:
                ax_NA_shared.plot(t, res['N_A'], label=label, linewidth = 4, color = colors_sim[i_m], ls = styles_sim[i_m])

            # pi
            ax1 = ax_pi.plot(t, res['pi'][0, :], label=label, linewidth = 4, color=colors_sim[i_m])
            ax_pi.plot(t, res['pi'][10, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())
            ax_pi.plot(t, res['pi'][-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            # ax_pi.plot(t, 1e8*np.exp(-(p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            if p.memory == 0:
                ax_pi.axhline(pi_star, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            else:
                ax_pi.axhline(pi_star, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # lambda_B
            ax1 = ax_lamB.plot(t, p.b0*res['N_Ba'][0, :]/(N_B[0, :])/p.b0, label=label, linewidth = 2, color = colors_sim[i_m])
            ax_lamB.plot(t, p.b0*res['N_Ba'][10, :]/(N_B[10, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.5)
            ax_lamB.plot(t, p.b0*res['N_Ba'][-10, :]/(N_B[-10, :])/p.b0, label=label, linewidth = 2, color=ax1[0].get_color(), alpha=0.2)
            # ax_lamB.plot(t, ((p.h0 * (res['pi'][0, :]**p.hill / (res['pi'][0, :]**p.hill + p.pi_star**p.hill)) * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
            ax_lamB.plot(t, ((p.h0 * res['pi'][0, :]**p.hill * res['N_To']/(res['N_To'] + p.K_T))**(-1) + (p.b0)**(-1))**(-1)/p.b0, linestyle='--', color=ax1[0].get_color(), linewidth = 2)
            ax_lamB.axhline(p.b0/p.b0, 0, T, linewidth = 1, linestyle='--', color='grey', alpha=0.8)
    
            # T
            ax1 = ax_NT.semilogy(t, N_T, label=label, linewidth = 3)
            ax_NT.semilogy(t, res['N_To'], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())

            # B
            # activated = N_B > 1.5
            N_B_total = compute_N_B_tot(res)
            L_act = compute_L_act(res)
            # ax1 = ax_NB.plot(t, N_B_total, label=label, linewidth = 4, color = ax1[0].get_color())
            ax1 = ax_NB.plot(t, N_B[0, :], label=label, linewidth = 4, alpha=1, color=colors_sim[i_m])
            ax_NB.plot(t, N_B[10, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())
            ax_NB.plot(t, N_B[-10, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())
            # ax_NB.plot(t, res['N_Bo'][0, :], label=label, linewidth = 2, linestyle='--', color=ax1[0].get_color())
            # ax_NB.plot(t, res['N_Ba'][0, :], label=label, linewidth = 2, linestyle=':', color=ax1[0].get_color())
            # ax_NB.plot(t, N_B_total, label=label, linewidth = 3, linestyle='-', alpha=0.8, color=ax1[0].get_color())
            # ax_NB.plot(t, 1e-2*np.exp((p.b0)*t), label=label, linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            # Y
            ax1 = ax_Y.plot(t, N_B_total, linewidth = 4, color=colors_sim[i_m], label=label)
            # ax_Y.semilogy(t, res['N_Bo'][0, :] + res['N_Ba'][0, :], label=label, linewidth = 2, alpha=1, color=ax1[0].get_color())
            # ax_Y.semilogy(t, res['N_Bo'][10, :] + res['N_Ba'][10, :], label=label, linewidth = 2, alpha=0.5, color=ax1[0].get_color())
            # ax_Y.semilogy(t, res['N_Bo'][-1, :] + res['N_Ba'][-1, :], label=label, linewidth = 2, alpha=0.2, color=ax1[0].get_color())

            # Z
            Z_B_total = compute_potency_t(res)
            Z_B = N_B * res['weights'][:, None] * np.exp(-res['DG'][:, None])
            Z_B[N_B <= 2.0] = np.nan

            ax1 = ax_Z.plot(t, Z_B_total, linewidth = 4, color=colors_sim[i_m], label=label)
            ax_Z_shared.plot(t[Z_B_total>0], Z_B_total[Z_B_total>0], linewidth = 4, color=colors_sim[i_m], label=label)
            # ax_Z.semilogy(t, Z_B[0, :], label=label, linewidth = 3, alpha=0.5, color=ax1[0].get_color())
            if p.memory == 0:
                lambda_prime = p.lambda_A*p.eta*p.beta_star*(1-0.5)
                # ax_Z_shared.plot(t[t<4.5], 1e-12*np.exp((p.b0 + lambda_prime)*t[t<4.5]), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
                # ax_Z_shared.plot(t[t<8], 2e-2*np.exp((p.b0)*t[t<8]), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            # else:
                # ax_Z_shared.plot(t, 5e1*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)

            best_DG_mem = np.argmax(N_B[:, 0]>0)
            if p.memory ==0:
                tstar = t[N_B[best_DG_mem, :]>N_B[best_DG_mem, 0]+0.06*N_B[best_DG_mem, 0]][0]
            else:
                tstar = t[N_B[best_DG_mem, :]>N_B[best_DG_mem, 0]+0.01*N_B[best_DG_mem, 0]][0]
            tpeak_pi = t[res['pi'][0, :]==np.max(res['pi'][0, :])][0]
            tpeak_NA = t[res['N_A']==np.max(res['N_A'])][0]
            tpeak = tpeak_NA + (tpeak_pi - tpeak_NA)/2
            print(tpeak)
            ton = t[(t<tpeak) & (t>tstar)]
            toff = t[(t>=tpeak)]
            # pioff = res['pi'][0, :][(t>=tpeak)]
            # tend = toff[pioff<pi_star][0]
            DGpeak = 0.78*p.lambda_A/p.eta * (tpeak - tstar)
            DG_on = 0.78*p.lambda_A/p.eta * (ton - tstar)
            DG_off = p.b0 /p.eta * (tpeak - toff) + DGpeak
            # DG_off2 = p.b0/p.eta * (tend - toff)
            DG_off_null = p.delta_A/p.eta * (tpeak - toff) + DGpeak
            ax_DG.plot(ton, DG_on, linewidth = 3, color=colors_sim[i_m], ls = 'dotted')
            ax_DG.plot(toff[DG_off>p.DG_min], DG_off[DG_off>p.DG_min], linewidth = 3, color=colors_sim[i_m], ls = 'dashed')
            # ax_DG.plot(toff[DG_off>p.DG_min], DG_off2[DG_off>p.DG_min], linewidth = 3, color=colors_sim[i_m], ls = 'dashed')
            ax_DG.plot(toff[DG_off_null>p.DG_min], DG_off_null[DG_off_null>p.DG_min], linewidth = 3, color=my_brown, ls = 'dashed')

            if memory == 0:
                innate_proxy = np.exp(p.lambda_innate * t)          # external, no feedback
                innate = (1.0 + p.threshold_innate / (innate_proxy + 1e-30))**(-1)
                ax_pb_shared.plot(t, innate, label=label, linewidth = 4, color=colors_sim[i_m], ls = '--')

            # active = (res['pi']> pi_star)                      # pi_c = b0/h0, see below
            # clone_mass = (res['N_Bo'] + res['N_Ba']) * active
            # Ab_proxy = np.sum(res['weights'][:, None] * clone_mass)
            # pb = (1.0 + p.Z_c / (Ab_proxy + 1e-30))**(-1)
            # ax_pb_shared.plot(t, pb, label=label, linewidth = 2, color=colors_sim[i_m])
            pb = (1.0 + p.Z_c / (Z_B_total + 1e-30))**(-1)
            ax_pb_shared.plot(t, pb, label=label, linewidth = 4, color=colors_sim[i_m])

            
            # Formatting
            # ax_NA.set_ylabel('$N_A$', fontsize = 16)
            ax_NA.set_xticklabels([])
            ax_NA.set_ylim(bottom = 1e0, top = 1e11)
            ax_NA.set_xlim(0, T)
            ax_NA.set_yscale('log')
            ax_NA.tick_params(axis='y', labelsize=30)
            ax_NA.tick_params(axis='x', labelsize=30)
            fig_NA.savefig(os.path.join(output_plot, f'N_A.pdf'), dpi=150)

            # ax_pi.set_ylabel(r'$\pi$', fontsize = 16)
            ax_pi.set_xticklabels([])
            ax_pi.set_ylim(bottom = 5e-1, top = 1e4)
            ax_pi.set_xlim(0, T)
            ax_pi.set_yscale('log')
            ax_pi.tick_params(axis='y', labelsize=30)
            ax_pi.tick_params(axis='x', labelsize=30)
            fig_pi.savefig(os.path.join(output_plot, f'pi.pdf'), dpi=150)

            # ax_lamB.set_ylabel(r'$\pi$', fontsize = 16)
            ax_lamB.set_xticklabels([])
            ax_lamB.set_ylim(bottom = -0.1, top = 1.1)
            ax_lamB.set_xlim(0, T)
            ax_lamB.tick_params(axis='y', labelsize=30)
            ax_lamB.tick_params(axis='x', labelsize=30)
            fig_lamB.savefig(os.path.join(output_plot, f'lamB.pdf'), dpi=150)

            ax_NB.axhline(1.0, color='k', linestyle='--', alpha=0.5)
            # ax_NB.set_xlabel('Time', fontsize = 16)
            # ax_NB.set_ylabel('B cells', fontsize = 16)
            ax_NB.set_xticklabels([])
            ax_NB.set_ylim(bottom = 5e-1, top = 1e6)
            ax_NB.set_xlim(0, T)
            ax_NB.set_yscale('log')
            ax_NB.tick_params(axis='y', labelsize=30)
            ax_NB.tick_params(axis='x', labelsize=30)
            fig_NB.savefig(os.path.join(output_plot, f'N_B.pdf'), dpi=150)

            ax_DG.axhline(p.DG_min, tstar, toff[DG_off>p.DG_min][-1], linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            # ax_DG.set_xlabel('Time', fontsize = 16)
            # ax_DG.set_ylabel('B cells', fontsize = 16)
            # ax_DG.set_xticklabels([])
            ax_DG.set_ylim(bottom = p.DG_min-0.1, top = 6)
            ax_DG.set_xlim(0, T)
            # ax_DG.set_yscale('log')
            ax_DG.tick_params(axis='y', labelsize=30)
            ax_DG.tick_params(axis='x', labelsize=30)
            fig_DG.savefig(os.path.join(output_plot, f'DG.pdf'), dpi=150)

            ax_Y.plot(t, 1e-2*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            ax_Y.axhline(1.0, color='k', linestyle='--', alpha=0.5)
            # ax_Y.set_xlabel('Time')
            # ax_Y.set_ylabel(r'Yield, $\bar{N}$', fontsize = 16)
            # ax_Y.set_xticklabels([])
            ax_Y.set_xlabel('Time', fontsize = 16)
            # ax_Y.set_ylim(bottom = 5e-1, top = 1e10)
            ax_Y.set_xlim(0, T)
            ax_Y.set_yscale('log')
            ax_Y.tick_params(axis='y', labelsize=30)
            ax_Y.tick_params(axis='x', labelsize=30)
            # ax_Y.legend(fontsize=14)
            fig_Y.savefig(os.path.join(output_plot, f'Y.pdf'), dpi=150)

            # if p.memory == 0:
                # lambda_prime = p.lambda_A*p.eta*p.beta_star*(1-0.5)
                # ax_Z.plot(t, 1e-17*np.exp((p.b0 + lambda_prime)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
                # ax_Z.plot(t, 2e-2*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            # else:
                # ax_Z.plot(t, 2e0*np.exp((p.b0)*t), linewidth = 1, linestyle='--', color='grey', alpha=0.8)
            ax_Z.axhline(1.0, linewidth = 1, linestyle='--', color='k', alpha=1.0)
            ax_Z.axhline(p.Z_c, linewidth = 1, linestyle='--', color='k', alpha=1.0)
            # ax_Z.set_xlabel('Time')
            # ax_Z.set_ylabel('Potency, $Z$', fontsize = 16)
            # ax_Z.set_xticklabels([])
            ax_Z.set_xlabel('Time', fontsize = 16)
            ax_Z.set_ylim(bottom = 5e-3, top = 1e7)
            ax_Z.set_xlim(0, T)
            ax_Z.set_yscale('log')
            ax_Z.tick_params(axis='y', labelsize=30)
            ax_Z.tick_params(axis='x', labelsize=30)
            # ax_Z.legend(fontsize=14)
            fig_Z.savefig(os.path.join(output_plot, f'Z.pdf'), dpi=150)

        output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}/{subproject}/{subsubproject}/'
        os.makedirs(output_plot, exist_ok=True)

        # Formatting
        # ax_NA_shared.set_ylabel('$N_A$', fontsize = 16)
        ax_NA_shared.set_xticklabels([])
        ax_NA_shared.set_ylim(bottom = 1e0, top = 1e11)
        ax_NA_shared.set_xlim(0+0.5, T-3)
        ax_NA_shared.set_yscale('log')
        ax_NA_shared.tick_params(axis='y', labelsize=30)
        ax_NA_shared.tick_params(axis='x', labelsize=30)
        fig_NA_shared.savefig(os.path.join(output_plot, f'N_A_shared.pdf'), dpi=150)
        
        ax_Z_shared.axhline(p.Z_c, linewidth = 1, linestyle='--', color='k', alpha=1.0)
        # ax_Z_shared.set_xlabel('Time')
        # ax_Z_shared.set_ylabel('Potency, $Z$', fontsize = 16)
        # ax_Z_shared.set_xticklabels([])
        # ax_Z_shared.set_xlabel('Time', fontsize = 16)
        ax_Z_shared.set_ylim(bottom = 5e-1, top = 1e6)
        ax_Z_shared.set_xlim(0+0.5, T-3)
        ax_Z_shared.set_yscale('log')
        ax_Z_shared.tick_params(axis='y', labelsize=30)
        ax_Z_shared.tick_params(axis='x', labelsize=30)
        # ax_Z_shared.legend(fontsize=14)
        fig_Z_shared.savefig(os.path.join(output_plot, f'Z_shared.pdf'), dpi=150)

        # ax_pb_shared.axhline(p.Z_c, linewidth = 1, linestyle='--', color='k', alpha=1.0)
        # ax_pb_shared.set_xlabel('Time')
        # ax_pb_shared.set_ylabel('Potency, $Z$', fontsize = 16)
        # ax_pb_shared.set_xticklabels([])
        # ax_pb_shared.set_xlabel('Time', fontsize = 16)
        # ax_pb_shared.set_ylim(bottom = 5e-1, top = 1e6)
        ax_pb_shared.set_xlim(0+0.5, T-3)
        # ax_pb_shared.set_yscale('log')
        ax_pb_shared.tick_params(axis='y', labelsize=30)
        ax_pb_shared.tick_params(axis='x', labelsize=30)
        # ax_pb_shared.legend(fontsize=14)
        fig_pb_shared.savefig(os.path.join(output_plot, f'pb_shared.pdf'), dpi=150)
