"""
Study N_B_tot (total activated B cells) and L_act (number of activated clones)
across the D~1 crossover.

Uses functions from ep_meanfield_sim.py.
"""

from tkinter.messagebox import IGNORE
from turtle import mode

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
subsubproject = 'b0'

if __name__ == '__main__':
    output_plot = f'/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/{project}/{model}/{submodel}//{subproject}/{subsubproject}/'
    os.makedirs(output_plot, exist_ok=True)
    # Default parameters
    base = dict(N_A0=1.0, delta_A=4.0, lambda_A = 6.,
                k_on=1e2*1e6*1e6*24*3600/N_Avg, delta_pi=0.1, Theta=1000.0,
                hill=1.0, h0=1e-0, beta_star=2.5, K_T = 1e5,
                delta_T=0.00, Tcell_growth_factor=2.0,
                tau_eng=0.1, sigma=1.0, delta_B=0.00,
                DG_min=0.0, DG_max=3.0, M=20,
                Omega_0=1.0, T_lim = True, N_T0=1e4,
                memory=False
    )
    T = 12
    N_ensemble = 1
    # print(compute_N_B_tot(res))
    # ============================================================
    # Scan N_T: move t_D relative to dynamics
    # ============================================================

    fig_potency, ax_potency = plt.subplots(figsize=(6, 5))

    # Storage for summary
    summary = []
    b0s = np.linspace(1, 8, 10)
    potencies = []
    Byields = []
    for b0 in b0s:
        print(f"Running simulation for b0={b0}...")
        potencies_i = np.zeros(N_ensemble)
        Byields_i = np.zeros(N_ensemble)
        for i in range(N_ensemble):
            # print(f" {i+1}/{N_ensemble}...")
            p = Parameters(**base, b0=b0)
            # res = run_simulation(p=p, t_span=(0, T), mode='grid')
            res = run_simulation_semicomplete(p=p, t_span=(0, T), mode='stochastic')

            Byield = compute_yield(res)
            potency = compute_potency(res, time_index=-1)
            Byields_i[i] = Byield
            potencies_i[i] = potency

        Byields.append(np.mean(Byields_i))    
        potencies.append(np.mean(potencies_i))
        label = f'$b_0={b0}$'

    ax_potency.plot(b0s, potencies, label=r'$\mathrm{Potency}$', color='tab:blue', lw = 2)
    ax_potency.plot(b0s, Byields, label=r'$\mathrm{Yield}$', color='tab:red', lw = 2)

    # Formatting
    ax_potency.set_xlabel(r'$b_0$')
    # ax_potency.set_ylabel('Frequency')
    ax_potency.legend(fontsize=7)
    ax_potency.set_yscale('log')
    ax_potency.set_xscale('linear')
    ax_potency.tick_params(axis='both', which='major', labelsize=14)
    ax_potency.legend(fontsize=14)
    plt.tight_layout()
    fig_potency.savefig(os.path.join(output_plot, f'potency_b0.pdf'), dpi=150, bbox_inches='tight')