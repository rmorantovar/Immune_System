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

project = 'exponential_proofreading'
model = 'meanfield'
subproject = 'zipf'
subsubproject = 'h0'

for submodel in ['semicomplete', 'null']:
    if __name__ == '__main__':
        output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}/{subproject}/{subsubproject}/'
        os.makedirs(output_plot, exist_ok=True)
        # Default parameters
        base = dict(N_A0=1.0, lambda_A = 6.0, delta_A=6.0,
                    k_on=1e0*1e6*1e6*24*3600/N_Avg, delta_pi=0.1,
                    hill=1.0, beta_star=2.5, K_T = 1e4,
                    delta_T=0.00, Tcell_growth_factor=2.0,
                    tau_eng=0.1, b0=2.0, delta_B=0.00,
                    DG_min=0.0, DG_max=8.0, M=40,
                    Omega_0=1.0, T_lim = True, N_T0 = 1e6
        )
        T = 15
        # print(compute_N_B_tot(res))
        # ============================================================
        # Scan N_T: move t_D relative to dynamics
        # ============================================================
        colors_sim = [my_red, my_blue2, my_purple2, my_gold, my_brown, my_blue, my_green, 'tab:orange', my_purple, my_cyan]
        colors_sim = plt.cm.viridis(np.linspace(0, 1, 50))
        N_ensemble = 1
        fig_zipf, ax_zipf = plt.subplots(figsize=(8, 5))

        print(f"Running simulation")
        etas = [1.0]
        for eta in etas:
            print(f"... for eta={eta:.1g}")
            
            h0s = np.array([base['b0'], 10*base['b0']])
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
                    ranks = np.linspace(1, 100, 100)
                    sizes = np.zeros_like(ranks)
                    for i_ensemble in range(N_ensemble):
                        if submodel == 'semicomplete':
                            res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='grid')
                        if submodel == 'null':
                            res = run_simulation_null(p=p, t_span=(0, T), mode='grid')

                        ranks_i, sizes_i = compute_zipf(res, time_index=-1)
                        sizes += sizes_i[:len(ranks)]
            

                    ax_zipf.loglog(ranks, sizes/N_ensemble, marker='o', ls = '', ms = 4, alpha = 0.5, markeredgecolor= 'k', markeredgewidth=0.5)
                
            zeta = p.eta*(p.b0/p.lambda_A + p.b0/np.min([p.b0, p.delta_A]))/p.beta_star
            ax_zipf.loglog(ranks, ranks**(-zeta), linestyle = '--', markersize=2, label=f'{zeta:.2g}', color = 'k')

            zeta_null = p.eta*(p.b0/p.lambda_A + p.b0/p.delta_A)/p.beta_star
            ax_zipf.loglog(ranks, ranks**(-zeta_null), linestyle = ':', markersize=2, label=f'{zeta_null:.2g}', color = 'grey')
        # Formatting
        
        # ax_zipf.loglog(ranks, ranks**(-2*p.eta*(p.b0/p.lambda_A + 1)/p.beta_star), linestyle = '--', markersize=2, color = my_blue, label='Zipf slope')
        ax_zipf.set_xlabel(r'Rank $k$', fontsize = 16)
        ax_zipf.set_ylabel(r'$N_k$', fontsize = 16)
        ax_zipf.set_ylim(bottom = 6e-2, top = 1.1e0)
        ax_zipf.set_xlim(left = 0.8e0, right = 1.1e2)
        ax_zipf.tick_params(labelsize=14)
        ax_zipf.legend(fontsize=12, title=r'$\zeta$', title_fontsize=14)
        fig_zipf.savefig(os.path.join(output_plot, f'Zipf.pdf'), dpi=150, bbox_inches='tight')