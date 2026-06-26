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
from funcs import *

project = 'exponential_proofreading'
model = 'meanfield'
subproject = 'zipf'
subsubproject = 'h0'
fig_kw = dict(figsize=(8 * 1.62, 8), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})


if __name__ == '__main__':
    fig, ax = plt.subplots(**fig_kw)
    for submodel in ['semicomplete', 'null']:
        output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}/{subproject}/{subsubproject}/'
        os.makedirs(output_plot, exist_ok=True)
        # Default parameters
        base = dict(N_A0=1.0, lambda_A = 6.0, delta_A=6.0,
                    k_on=1e0*1e6*1e6*24*3600/N_Avg, delta_pi=0.1,
                    hill=1.0, beta_star=2.3, K_T = 1e4,
                    delta_T=0.00, Tcell_growth_factor=2.0,
                    tau_eng=0.1, b0=2.0, delta_B=0.00,
                    DG_min=0.0, DG_max=8.0, M=40,
                    omega_0=1.0, T_lim = True, N_T0 = 1e6
        )
        T = 15
        # print(compute_N_B_tot(res))
        # ============================================================
        # Scan N_T: move t_D relative to dynamics
        # ============================================================
        colors_sim = [my_red, my_blue2, my_purple2, my_gold, my_brown, my_blue, my_green, 'tab:orange', my_purple, my_cyan]
        colors_sim = plt.cm.viridis(np.linspace(0, 1, 50))
        N_ensemble = 1
        fig_zipf, ax_zipf = plt.subplots(**fig_kw)

        print(f"Running simulation")
        if submodel == 'semicomplete':
            memories = [0, 1]
            my_colors = ['limegreen', my_blue2]
            my_markers = ['*', 'o']
            my_ms = [18, 12]
        else:
            memories = [0]
            my_colors = [my_brown]
            my_markers = ['^']
            my_ms = [12]

        etas = [1.0]
        for eta in etas:
            print(f"... for eta={eta:.1g}")
            
            h0s = np.array([10*base['b0']])
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
                
                for i_m, memory in enumerate(memories):
                    print(f"... ... for memory={memory}")
                    base['memory'] = memory

                    p = Parameters(**base)
                    
                    initial_memory_potencies_pi_star = []
                    final_primary_potencies_pi_star = []
                    initial_memory_yield_pi_star = []
                    final_primary_yield_pi_star = []
                    ranks = np.linspace(1, 100, 100)
                    sizes = np.zeros_like(ranks)
                    for i_ensemble in range(N_ensemble):
                        if submodel == 'semicomplete':
                            res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='grid')
                        if submodel == 'null':
                            res = run_simulation_null(p=p, t_span=(0, T), mode='grid')

                        ranks_i, sizes_i = compute_zipf(res, time_index=-1)
                        sizes += sizes_i[:len(ranks)]
            

                    ax_zipf.plot(ranks, sizes/N_ensemble, marker=my_markers[i_m], ls = '', ms = my_ms[i_m], alpha = 1.0, markerfacecolor='None',  color = my_colors[i_m])
                    ax.plot(ranks, sizes/N_ensemble, marker=my_markers[i_m], ls = '', ms = my_ms[i_m], alpha = 1.0, markerfacecolor='None', color = my_colors[i_m])
                
        zeta = p.eta*(p.b0/p.lambda_A + p.b0/np.min([p.b0, p.delta_A]))/p.beta_star
        zeta_null = p.eta*(p.b0/p.lambda_A + p.b0/p.delta_A)/p.beta_star

        if submodel == 'semicomplete':
            ax_zipf.loglog(ranks, ranks**(-zeta), lw = 3, linestyle = '--', markersize=2, label=f'{zeta:.2g}', color = 'limegreen')
            ax_zipf.loglog(ranks, ranks**(-2*zeta), lw = 3, linestyle = '--', markersize=2, label=f'{zeta:.2g}', color = my_blue)
        else:
            ax_zipf.loglog(ranks, ranks**(-zeta), lw = 3, linestyle = '--', markersize=2, label=f'{zeta:.2g}', color = 'limegreen')
            ax_zipf.loglog(ranks, ranks**(-zeta_null), lw = 3, linestyle = '--', markersize=2, label=f'{zeta_null:.2g}', color = my_brown)


        # Formatting
        
        # ax_zipf.loglog(ranks, ranks**(-2*p.eta*(p.b0/p.lambda_A + 1)/p.beta_star), linestyle = '--', markersize=2, color = my_blue, label='Zipf slope')
        my_plot_layout(ax=ax, yscale='log', xscale='log', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
        # ax_zipf.set_xlabel(r'Rank $k$', fontsize = 16)
        # ax_zipf.set_ylabel(r'$N_k$', fontsize = 16)
        ax_zipf.set_ylim(bottom=2e-2, top=1.1)
        ax_zipf.set_xlim(left = 0.9, right=5e1)
        ax_zipf.tick_params(labelsize=30)
        ax_zipf.legend(fontsize=30, title=r'$\zeta$', title_fontsize=30)
        fig_zipf.savefig(os.path.join(output_plot, f'Zipf.pdf'), dpi=150)
    
    
    ax.loglog(ranks, ranks**(-zeta_null), lw = 3, linestyle = '-', markersize=2, label=f'{zeta_null:.2g}', color = my_brown)
    ax.loglog(ranks, ranks**(-zeta), lw = 3, linestyle = '-', markersize=2, label=f'{zeta:.2g}', color = 'limegreen')
    ax.loglog(ranks, ranks**(-2*zeta), lw = 3, linestyle = '-', markersize=2, label=f'{2*zeta:.2g}', color = my_blue)

    # ax.loglog(ranks, ranks**(-2*p.eta*(p.b0/p.lambda_A + 1)/p.beta_star), linestyle = '--', markersize=2, color = my_blue, label='Zipf slope')
    my_plot_layout(ax=ax, yscale='log', xscale='log', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    # ax.set_xlabel(r'Rank $k$', fontsize = 16)
    # ax.set_ylabel(r'$N_k$', fontsize = 16)
    ax.set_ylim(bottom=2e-2, top=1.1)
    ax.set_xlim(left = 0.9, right=5e1)
    ax.tick_params(which = 'both', labelsize=30)
    ax.legend(fontsize=30, title=r'$\zeta$', title_fontsize=30)
    output_plot_combined = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/combined/{subproject}/{subsubproject}/'
    os.makedirs(output_plot_combined, exist_ok=True)
    fig.savefig(os.path.join(output_plot_combined, f'Zipf.pdf'), dpi=150)