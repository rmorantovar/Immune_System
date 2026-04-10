"""
Study N_B_tot (total activated B cells) and L_act (number of activated clones)
across the D~1 crossover.

Uses functions from ep_meanfield_sim.py.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../../library/')
from lib_mf import*


output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/exponential_proofreading/mean_field_entropy/'
os.makedirs(output_plot, exist_ok=True)

# ============================================================
# Base parameters
# ============================================================

base = dict(
    N_A0=1.0,
    delta_A=0.1,
    k_on=1e0*1e6*1e6*24*3600/N_Avg,
    delta_pi=0.1,
    Theta=100.0,
    sigma=1.0,
    beta_star=2.5,
    delta_T=0.0,
    gamma=100.00,
    tau_eng=0.5,
    b0=2.0,
    delta_B=0.0,
    DG_min=0.0,
    DG_max=6.0,
    M=200,
    Omega_0=1.0,
    T_lim = 1,
    memory = False
)

T = 20.0

# ============================================================
# Scan N_T: move t_D relative to dynamics
# ============================================================

N_T_values = [1e3, 1e4, 1e5, 1e6]

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Storage for summary
summary = []

for N_T in N_T_values:
# for lam_A in [5.0, 6.0]:
    p = Parameters(**base, N_T0=N_T, lambda_A = 6.)
    res = run_simulation(p=p, t_span=(0, T), t_eval=np.linspace(0, T, 1000))

    t = res['t']
    N_B_tot = compute_N_B_tot(res)
    L_act = compute_L_act(res)
    t_D = find_t_D(res)

    summary.append({
        'N_T': N_T,
        'lambda_A': 6.0,
        't_D': t_D,
        'N_B_tot_final': N_B_tot[-1],
        'L_act_final': L_act[-1],
        'D_final': res['D'][-1],
    })

    label = f'$N_T={N_T}$'

    # (a) N_B_tot
    ax = axes[0, 0]
    ax.semilogy(t, N_B_tot, label=label)

    # (b) L_act
    ax = axes[0, 1]
    ax.semilogy(t, L_act + 0.1, label=label)

    # (c) D(t)
    ax = axes[0, 2]
    ax.semilogy(t, res['D'] + 1e-10, label=label)

    # (d) N_T_free / N_T
    ax = axes[1, 0]
    ax.plot(t, res['N_T_free'] / N_T, label=label)

    # (e) Growth rate: d(ln N_B_tot)/dt
    ax = axes[1, 1]
    dln = np.gradient(np.log(np.maximum(N_B_tot, 1.0)), t)
    ax.plot(t, dln, label=label)

    # (f) N_B_tot / N_T
    ax = axes[1, 2]
    ax.plot(t, N_B_tot / N_T, label=label)


# Formatting
axes[0, 0].set_xlabel('Time')
axes[0, 0].set_ylabel('$N_B^{tot}$')
axes[0, 0].set_title('Total activated B cells')
axes[0, 0].legend(fontsize=7)

axes[0, 1].set_xlabel('Time')
axes[0, 1].set_ylabel('$L_{act}$')
axes[0, 1].set_title('Number of activated clones')
axes[0, 1].legend(fontsize=7)

axes[0, 2].set_xlabel('Time')
axes[0, 2].set_ylabel('$D(t)$')
axes[0, 2].set_title('Demand function')
axes[0, 2].axhline(1.0, color='k', linestyle=':', alpha=0.5)
axes[0, 2].legend(fontsize=7)

axes[1, 0].set_xlabel('Time')
axes[1, 0].set_ylabel('$N_T^{free}/N_T$')
axes[1, 0].set_title('Fraction free T cells')
axes[1, 0].axhline(0.5, color='k', linestyle=':', alpha=0.5)
axes[1, 0].legend(fontsize=7)

axes[1, 1].set_xlabel('Time')
axes[1, 1].set_ylabel('$d(\\ln N_B^{tot})/dt$')
axes[1, 1].set_title('Growth rate')
axes[1, 1].axhline(base['b0'], color='k', linestyle=':', alpha=0.5, label='$b0$')
axes[1, 1].legend(fontsize=7)

axes[1, 2].set_xlabel('Time')
axes[1, 2].set_ylabel('$N_B^{tot} / N_T$')
axes[1, 2].set_title('B cells per T cell')
axes[1, 2].legend(fontsize=7)

plt.suptitle('Effect of $N_T$ on activation dynamics', fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(output_plot, f'crossover_study_Tlim-{int(base["T_lim"])}.pdf'), dpi=150, bbox_inches='tight')
print(f"Saved: crossover_study_Tlim-{int(base['T_lim'])}.pdf")

# Print summary
print(f'\n{"N_T":>8} {"t_D":>8} {"N_B_tot":>12} {"L_act":>10} {"D_final":>10}')
for s in summary:
    print(f'{s["N_T"]:>8} {s["t_D"]:>8.2f} {s["N_B_tot_final"]:>12.1f} '
          f'{s["L_act_final"]:>10.1f} {s["D_final"]:>10.1f}')