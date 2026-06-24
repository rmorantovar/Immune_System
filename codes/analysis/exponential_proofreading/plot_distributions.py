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
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import os
import sys
sys.path.append('../../library/')
from funcs import *


# ── Paths & constants ─────────────────────────────────────────────────────────
project    = 'exponential_proofreading'
subproject = 'plots'
root_dir   = (
    f"/Users/robertomorantovar/Dropbox/Research_files/"
    f"Immune_system/{project}/{subproject}/mesin2020"
)
output_plot = (
    '/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/'
    f'Immune_System/_Repository/Figures/{project}/{subproject}'
)
os.makedirs(output_plot, exist_ok=True)

# ── Figure setup ──────────────────────────────────────────────────────────────

fig_kw = dict(figsize=(8 * 1.62, 8), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})

# ----------------------------------------------------------------------
# PARAMETERS  (replace with your values)
# ----------------------------------------------------------------------
zeta = 0.5     # EP proofreading strength (zeta = alpha/beta*)
sigma = np.sqrt(7)     # width sigma_E
alpha = 1.3    # EP tilt strength per round  (alpha = zeta * beta*)
beta_star = alpha / zeta   # beta* = 1/(k_B T) in DG units (for plotting only)
print('beta* =', beta_star)
mu0   = beta_star * sigma**2    # naive peak position (in DG units)
# the two normalization offsets from your Mathematica expr:
mu1  = mu0 - alpha*sigma**2   # offset entering the 2nd tilt prefactor (= mu1)

edge  = 0.0     # EVT edge  DG*  (no clones below this)
a1    = 6.0     # 1EP activation threshold (upper DG cutoff)
a2    = 3.0     # 2EP activation threshold

L0    = np.exp(0.5 * (mu0/sigma)**2)   # implied repertoire size (for amplitude)
xmax  = mu0 + 0.5*sigma   # x-axis max for plotting
# ----------------------------------------------------------------------

def Omega0(x):
    # Gaussian density of states, normalized so peak amplitude ~ L0/(sqrt(2pi)sigma)
    return (L0/(np.sqrt(2*np.pi)*sigma)) * np.exp(-(x-mu0)**2/(2*sigma**2))

def q1(x):
    return Omega0(x) * np.exp(-alpha*x) * np.exp(-0.5*sigma**2*alpha**2 + alpha*mu0)

def q2(x):
    return q1(x) * np.exp(-alpha*x) * np.exp(-0.5*sigma**2*alpha**2 + alpha*mu1)

x = np.linspace(-0.5, xmax, 1200)

# tab10-style colors (Mathematica ColorData[97] analog)
c0, c1, c2 = 'gray', my_red, my_blue   # naive, 1EP, 2EP
colors = [c0, c1, c2]

fig, ax = plt.subplots(**fig_kw)
fig2, ax2 = plt.subplots(**fig_kw)

# --- curves ---
ax.plot(x, Omega0(x), color=c0, lw=2.2, label=r"$0$")
ax.plot(x, q1(x),     color=c1, lw=2.2, label=r"$1$")
ax.plot(x, q2(x),     color=c2, lw=2.2, label=r"$2$")

# ax2.plot(x, , color=c0, lw=2.2, label=r"$0$")

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

# --- cosmetics: left+bottom spines only (matches your Frame->{{T,F},{T,F}}) ---
my_plot_layout(ax=ax, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
# ax.spines["top"].set_visible(False)
# ax.spines["right"].set_visible(False)
ax.set_xlim(-0.6, xmax)
ax.set_ylim(0, ymax)
# ax.set_xlabel(r"Binding energy $\Delta G$", fontsize=13)
# ax.set_ylabel(r"$\Omega(\Delta G)$", fontsize=13)
ax.tick_params(axis='both', labelsize=30)

leg = ax.legend(title=r"$\mathrm{EP\, rounds}$", loc=0, frameon=False, fontsize=20, title_fontsize=20)

# fig.tight_layout()
fig.savefig(output_plot + "/ep_plot.pdf", bbox_inches='tight', transparent=.5)
# fig.savefig(output_plot + "/ep_plot.png", dpi=150, bbox_inches="tight")
print("saved")