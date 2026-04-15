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


if __name__ == '__main__':
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_Tcell_jam_2/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, delta_A=0.01, lambda_A = 6.,
                k_on=1e2*1e6*1e6*24*3600/N_Avg, delta_pi=24., Theta=10000.0,
                hill=3.0, sigma=1.0, beta_star=2.5, K_T = 100000.,
                delta_T=0.0, h0=0.002, N_T0=1e6,
                tau_eng=0.01, b0=2.0, delta_B=0.0,
                DG_min=0.0, DG_max=2.5, M=20,
                Omega_0=1.0, T_lim = 0
    )
    T = 16
    N_ensemble = 1000
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    fig_potency, ax_potency = plt.subplots(figsize=(6, 5))

    # Storage for summary
    summary = []

    for memory in [False, True]:
        print(f"Running simulation for memory={memory}")
        potencies = np.zeros(N_ensemble)
        for i in range(N_ensemble):
            print(f" {i+1}/{N_ensemble}...")
            p = Parameters(**base,  memory=memory)
            # res = run_simulation(p=p, t_span=(0, T), mode='grid')
            res = run_simulation(p=p, t_span=(0, T), mode='stochastic')

            potency = compute_potency(res, time_index=-1)
            potencies[i] = potency

        label = memory
        log_potencies = (np.log(potencies) - np.mean(np.log(potencies)))/np.std(np.log(potencies)) # normalize and center
        ax_potency.hist(log_potencies, bins=20, alpha=0.5, label=label, density=True)

    ax_potency.plot(np.linspace(-4, 4, 100), stats.norm.pdf(np.linspace(-3, 3, 100)), 'k--', label='Normal(0,1)')


    # Formatting
    ax_potency.set_xlabel('Potency'); ax_potency.set_ylabel('Frequency')
    ax_potency.set_title('Potency distribution (final)')
    ax_potency.legend(fontsize=7)
    ax_potency.set_yscale('log')
    plt.tight_layout()
    fig_potency.savefig(os.path.join(output_plot, f'Potency_Tlim-{int(base["T_lim"])}_Mem-{int(p.memory)}.pdf'), dpi=150, bbox_inches='tight')