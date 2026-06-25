"""
EP repertoire figure: naive density of states Omega0 and its 1st/2nd tilted
(EP) versions, with the EVT edge, per-round activation thresholds a1/a2, and
the *active* population shaded from the edge to each threshold.

Drop your own parameter values in the PARAMETERS block. The curve definitions
mirror the Mathematica expressions:
    n=0:  Omega0(x)
    n=1:  Omega0(x) * exp(-a x) * exp(-1/2 s1^2 a^2 + a mu01)
    n=2:  (1EP) * exp(-a x) * exp(-1/2 s1^2 a^2 + a mu11)
"""

import numpy as np
from scipy.stats import differential_entropy
from scipy.stats import norm
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import os
import sys
sys.path.append('../../library/')
from funcs import *


# ── Paths & constants ─────────────────────────────────────────────────────────
project    = 'exponential_proofreading'
subproject = 'data'
root_dir   = (
    f"/Users/robertomorantovar/Dropbox/Research_files/"
    f"Immune_system/{project}/{subproject}/mesin2020"
)
subproject = 'plots'
output_plot = (
    '/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/'
    f'Immune_System/_Repository/Figures/{project}/{subproject}'
)
os.makedirs(output_plot, exist_ok=True)

for alpha in [1.0, 1.3, 2.0, 3.0]:# EP tilt strength per round 
    # ── Figure setup ──────────────────────────────────────────────────────────────

    fig_kw = dict(figsize=(8 * 1.62, 8), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})

    # ----------------------------------------------------------------------
    # PARAMETERS  (replace with your values)
    # ----------------------------------------------------------------------

    zeta = 0.5     # EP proofreading strength (zeta = alpha/beta*)
    sigma = np.sqrt(7)     # width sigma_E
    beta_star = alpha / zeta   # beta* = 1/(k_B T) in DG units (for plotting only)
    print('beta* =', beta_star)
    mu0   = beta_star * sigma**2    # naive peak position (in DG units)
    # the two normalization offsets from your Mathematica expr:
    mu1  = mu0 - alpha*sigma**2   # offset entering the 2nd tilt prefactor (= mu1)
    mu2  = mu0 - 2*alpha*sigma**2   # offset entering the 2nd tilt prefactor (= mu2)
    print('mu0, mu1, mu2 =', mu0, mu1, mu2)

    edge  = 0.0     # EVT edge  DG*  (no clones below this)
    a1    = 6.0     # 1EP activation threshold (upper DG cutoff)
    a2    = 3.0     # 2EP activation threshold

    L0    = np.exp(0.5 * (mu0/sigma)**2)   # implied repertoire size (for amplitude)
    xmax  = mu0 + 0.5*sigma   # x-axis max for plotting

    varphi = norm.pdf      # φ(x)
    Phi = norm.cdf      # Φ(x)

    # ----------------------------------------------------------------------

    def Omega0(x):
        # Gaussian density of states, normalized so peak amplitude ~ L0/(sqrt(2pi)sigma)
        return (1/(np.sqrt(2*np.pi)*sigma)) * np.exp(-(x-mu0)**2/(2*sigma**2))

    def q1(x):
        return Omega0(x) * np.exp(-alpha*x) * np.exp(-0.5*sigma**2*alpha**2 + alpha*mu0)

    def q2(x):
        return q1(x) * np.exp(-alpha*x) * np.exp(-0.5*sigma**2*alpha**2 + alpha*mu1)

    x = np.linspace(-0.5, xmax, 1200)
    a = np.linspace(0.1, 3*xmax/4, 1200)

    # tab10-style colors (Mathematica ColorData[97] analog)
    c0, c1, c2 = 'gray', my_red, my_blue   # naive, 1EP, 2EP
    colors = [c0, c1, c2]

    fig, ax = plt.subplots(**fig_kw)
    fig2, ax2 = plt.subplots(**fig_kw)
    fig3, ax3 = plt.subplots(**fig_kw)

    # --- curves ---
    ax.plot(x, Omega0(x), color=c0, lw=2.2, label=r"$0$")
    ax.plot(x, q1(x),     color=c1, lw=2.2, label=r"$1$")
    ax.plot(x, q2(x),     color=c2, lw=2.2, label=r"$2$")

    D_max = 0.5*(sigma*beta_star)**2
    phi_min = np.exp(-D_max)
    ax2.plot(x, np.ones_like(x)*phi_min, color='k', ls = '--', lw=3)

    phi_max_primary = np.exp(-zeta**2*D_max)
    phi_1_2_recall = 0.5*np.exp(-(D_max - np.sqrt(2/3.14*2*D_max)))
    ax2.plot(x, np.ones_like(x)*phi_max_primary, color=c1, ls = '--', lw=3, zorder=-10)
    ax2.plot(x, np.ones_like(x)*phi_1_2_recall, color=c2, ls = '--', lw=3, zorder=-10)
    
    print(D_max, zeta**2*D_max, D_max - np.sqrt(2/3.14*2*D_max))

    omega1 = (a - mu1)/sigma
    D1 = zeta**2*D_max + alpha*sigma*varphi(omega1)/Phi(omega1) - np.log(Phi(omega1))
    phi1 = np.exp(-D1)
    ax2.plot(a[phi1>phi_min], phi1[phi1>phi_min], color=c1, lw=4, label=r"$\mathrm{primary}$")
   
    omega2 = (a - mu2)/sigma
    D2 = D_max + np.sqrt(2*D_max)*(varphi(0) - varphi(omega2))/(Phi(omega2)-0.5) - np.log(Phi(omega2)-0.5)
    phi2 = np.exp(-D2)

    zp  = (edge - mu2)/sigma          # lower (edge) standardized bound
    zpp = (a    - mu2)/sigma          # upper (gate) standardized bound
    Z2  = Phi(zpp) - Phi(zp)
    g2  = varphi(zp) - varphi(zpp)
    m1_2 = g2 / Z2
    delta2 = (mu2 - mu0)/sigma        # = -2*alpha*sigma  (note: from mu0, not edge)
    D2 = 0.5*delta2**2 + delta2*m1_2 - np.log(Z2)
    phi2 = np.exp(-D2)
    ax2.plot(a[phi2>phi_min], phi2[phi2>phi_min], color=c2, lw=4, label=r"$\mathrm{recall}$")

    # --- filled active regions: edge -> a_n ---
    m1 = (x >= edge) & (x <= a1)
    m2 = (x >= edge) & (x <= a2)
    ax.fill_between(x[m1], 0, q1(x[m1]), color=c1, alpha=0.22, lw=0)
    ax.fill_between(x[m2], 0, q2(x[m2]), color=c2, alpha=0.30, lw=0)

    # --- edge and threshold markers ---
    ax.axvline(edge, color="0.35", lw=1.2, ls="-")
    ax.axvline(a1,   color=c1, lw=1.0, ls="--", alpha=0.8)
    ax.axvline(a2,   color=c2, lw=1.0, ls="--", alpha=0.8)

    ymax = Omega0(mu0)*1.05
    # ax.text(edge, ymax*0.98, r"$\Delta G^*$", color="0.35", ha="center", va="top", fontsize=12)
    # ax.text(a1, ymax*0.55, r"$a_1$", color=c1, ha="left", va="center", fontsize=12)
    # ax.text(a2, ymax*0.78, r"$a_2$", color=c2, ha="left", va="center", fontsize=12)

    data = pd.read_csv(root_dir + "/scaling_results.csv")
    primary_data = data[data['response'] == 'naive']
    recall_data = data[data['response'] == 'recall']
    
    for row in primary_data.itertuples():
        S_c_row = row.S
        N_sample = row.barN

        D_row = D_max - S_c_row
        phi_row = np.exp(-D_row)
        a_inferred = a[phi1<phi_row][-1]
        mask = a > 0                       # monotone region only
        ainv = np.interp(phi_row, phi1[mask], a[mask])   # phi1 must be increasing on mask      
        ax2.scatter(ainv, phi_row, color=c1, s=150, edgecolors='k', alpha=0.8, zorder = 20)        
        ax3.scatter(a_inferred, np.log(N_sample), color=c1, s=150, edgecolors='k', alpha=0.8)
    
    allN = recall_data['barN'].values
    my_norm = mpl.colors.LogNorm(vmin=allN.min(), vmax=allN.max())   # log: N spans decades
    cmap = plt.cm.managua
    for row in recall_data.itertuples():
        S_c_row = row.S
        N_sample = row.barN

        D_row = D_max - S_c_row
        phi_row = np.exp(-D_row)
        a_inferred = a[phi2<phi_row][-1]
        mask = a > 0                       # monotone region only
        ainv = np.interp(phi_row, phi2[mask], a[mask])   # phi2 must be increasing on mask
        col = cmap(my_norm(N_sample))
        ax2.scatter(ainv, phi_row, color=col, s=150, edgecolors='k', alpha=0.8, zorder = 20)
        ax3.scatter(a_inferred, np.log(N_sample), color=c2, s=150, edgecolors='k', alpha=0.8)

    axin2 = ax2.inset_axes([0.08, 0.6, 0.38, 0.38], zorder=30)  # [left, bottom, width, height] in axes fraction
    axin2.plot(10*x, np.ones_like(x)*phi_min, color='k', ls = '--', lw=3)
    axin2.plot(10*x, np.ones_like(x)*phi_1_2_recall, color=c2, ls = '--', lw=3)
    for row in recall_data.itertuples():
        S_c_row = row.S
        N_sample = row.barN

        D_row = D_max - S_c_row
        phi_row = np.exp(-D_row)
        col = cmap(my_norm(N_sample))
        axin2.scatter(N_sample, phi_row, color=col, s=150, edgecolors='k', alpha=0.8)
    

    # --- cosmetics: left+bottom spines only (matches your Frame->{{T,F},{T,F}}) ---
    my_plot_layout(ax=ax, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    # ax.spines["top"].set_visible(False)
    # ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.6, xmax)
    ax.set_ylim(bottom = 1e-3, top = ymax)
    # ax.set_xlabel(r"Binding energy $\Delta G$", fontsize=13)
    # ax.set_ylabel(r"$\Omega(\Delta G)$", fontsize=13)
    ax.tick_params(axis='both', labelsize=30)

    leg = ax.legend(title=r"$\mathrm{Response}$", loc=0, frameon=False, fontsize=20, title_fontsize=20)

    # fig.tight_layout()
    fig.savefig(output_plot + f"/ep_plot_alpha-{alpha}.pdf", transparent=.5)
    ax.set_yscale('log')
    fig.savefig(output_plot + f"/ep_plot_alpha-{alpha}_log.pdf", transparent=.5)


    # --- cosmetics: left+bottom spines only (matches your Frame->{{T,F},{T,F}}) ---
    my_plot_layout(ax=ax2, yscale='log', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    # ax2.spines["top"].set_visible(False)
    # ax2.spines["right"].set_visible(False)
    ax2.set_xlim(-0.6, 2*xmax/3)
    ax2.set_ylim(phi_min/5, 2)
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', labelsize=30)
    # ax2.set_xlabel(r"Binding energy $\Delta G$", fontsize=13)
    # ax2.set_ylabel(r"$\Omega(\Delta G)$", fontsize=13)
    leg = ax2.legend(title=r"$\mathrm{Response}$", loc=0, frameon=False, fontsize=20, title_fontsize=20)

    axin2.set_yscale('log')
    axin2.tick_params(axis='both', labelsize=24)
    axin2.set_facecolor('white')        # opaque background
    axin2.patch.set_alpha(1.0)          # ensure the patch isn't transparent
    axin2.set_zorder(50)                # whole inset above main-axes lines
    axin2.patch.set_visible(True)  

    # fig.tight_layout()
    fig2.savefig(output_plot + f"/ep_plot2_alpha-{alpha}.pdf")

    # --- cosmetics: left+bottom spines only (matches your Frame->{{T,F},{T,F}}) ---
    my_plot_layout(ax=ax3, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    # ax3.spines["top"].set_visible(False)
    # ax3.spines["right"].set_visible(False)
    ax3.set_xlim(1.5, 3.7)
    # ax3.set_ylim(-0.6, 2*xmax/3)
    # ax3.set_xlabel(r"Binding energy $\Delta G$", fontsize=13)
    # ax3.set_ylabel(r"$\Omega(\Delta G)$", fontsize=13)
    ax3.tick_params(axis='both', labelsize=30)

    leg = ax3.legend(title=r"$\mathrm{Response}$", loc=0, frameon=False, fontsize=20, title_fontsize=20)

    # fig.tight_layout()
    fig3.savefig(output_plot + f"/ep_plot3_alpha-{alpha}.pdf", transparent=.5)
   
    print("saved")