"""
mf_plotting.py -- Styles and reusable plotting helpers for the mean-field models.

Holds everything that used to be copy-pasted across the run scripts:

    * the colour palette and matplotlib rcParams,
    * small helpers (`new_fig`, `style_log_axis`, `pi_star_label`),
    * per-run panel plotters (`plot_NA`, `plot_pi`, `plot_potency`, ...),
    * the multi-panel diagnostic figures (`plot_results`, `plot_diagnostics`),
      which apply to the 'approx' model (they need res['D'], res['N_T_free']).

The driver composes these; nothing here runs a simulation.
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

from mf_lib import (compute_N_B_tot, compute_L_act, compute_zipf, compute_potency_t)

# ============================================================
# Palette
# ============================================================
my_red = np.array((228,75,41))/256.
# _colors_h0_default_naive = plt.cm.autumn(np.linspace(0, 0.7, 6))
# my_red = _colors_h0_default_naive[3]
my_purple = np.array((125, 64, 119)) / 256.
my_purple2 = np.array((116, 97, 164)) / 256.
my_green = np.array((125,165,38))/256.
# _colors_h0_default_naive = plt.cm.summer(np.linspace(0, 0.7, 6))
# my_green = _colors_h0_default_naive[3]
my_blue = np.array((76,109,166))/256.
# _colors_h0_default_memory = plt.cm.winter(np.linspace(0, 0.7, 6))
# my_blue = _colors_h0_default_memory[3]
my_gold = np.array((215, 139, 45)) / 256.
my_brown = np.array((182, 90, 36)) / 256.
my_blue2 = np.array((80, 141, 188)) / 256.
my_yellow = np.array((246, 181, 56)) / 256.
my_yellow2 = np.array((242, 192, 65)) / 256.
# my_green2 = np.array((158, 248, 72)) / 256.
my_green2 = 'darkgreen'
my_cyan = 'tab:cyan'
my_grey = np.array((128, 128, 128)) / 256.
antigen_color = my_yellow

my_green_a = np.array((159, 206, 99)) / 256.
my_green_b = np.array((79, 173, 91)) / 256.
my_green_c = np.array((94, 129, 63)) / 256.

# Default per-series colours / linestyles (index by loop counter).
colors_sim = [my_red, my_blue, my_purple2, my_gold, my_brown,
              my_blue, my_green, 'tab:orange', my_purple, my_cyan]
colors_mem = [my_green, my_green2]
colors_sweep = [my_brown, my_red, my_purple2, my_green, my_blue, my_cyan]
# colors_sweep = [my_red, my_green, my_cyan]
styles_sim = ['--', '-', ':', '-.', '-', '--', ':', '-.', '-', '--']
styles_mem = ['--', '-']

# Standard figure geometries.
FIG_KW = dict(figsize=(8 * 1.62, 8),
              gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
# FIG_KW = dict(figsize=(8 * 1.62, 10),
            #   gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
FIG_KW_WIDE = dict(figsize=(8 * 1.62, 5),
                   gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})


def use_style(usetex=True):
    """Apply the rcParams used across the project. Call once at startup."""
    plt.rcParams['text.usetex'] = usetex
    plt.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.tab10.colors)


# ============================================================
# Small helpers
# ============================================================
def new_fig(wide=False):
    """Return (fig, ax) with the standard geometry."""
    return plt.subplots(**(FIG_KW_WIDE if wide else FIG_KW))


def pi_star_label(pi_star):
    r"""LaTeX label like ${1.5 \cdot 10^{-3}}$ for a pi_star value."""
    e = np.floor(np.log10(pi_star))
    return r'${%.1f \cdot 10^{%d}}$' % (pi_star / 10 ** e, e)


def sweep_colors(n, cmap='summer', vmax=0.8):
    """
    Colours for an n-value parameter sweep.

    cmap : matplotlib colormap name (e.g. 'summer' for the green scale used in
           the h0 sweep). If None, falls back to the discrete `colors_sim`.
    vmax : upper end of the colormap range (0.8 matches the original h0 script).
    """
    if cmap is None:
        return [colors_sim[i % len(colors_sim)] for i in range(n)]
    return plt.get_cmap(cmap)(np.linspace(0, vmax, max(n, 1)))


def style_log_axis(ax, T, ylim=None, xlim=None, labelsize=44, hide_xticklabels=True):
    """Common log-y axis styling for the time-series panels."""
    # if hide_xticklabels:
    #     ax.set_xticklabels([])
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_xlim(*(xlim if xlim is not None else (0, T)))
    ax.set_yscale('log')
    ax.tick_params(axis='y', which='both', labelsize=labelsize, labelleft=True)
    ax.tick_params(axis='x', labelsize=labelsize, labelbottom = True)


def _experiment_tag(cfg):
    """Top-level folder describing the run: 'sweep_<param>' or 'single_run'."""
    name = cfg['sweep'][0]
    return f'sweep_{name}' if name is not None else 'single_run'


def _outdir(cfg, *parts):
    # Everything lands under outroot/<experiment_tag>/... so the folder itself
    # records whether the figures came from a sweep or a single run.
    path = os.path.join(cfg['outroot'], _experiment_tag(cfg),
                        *[str(x) for x in parts])
    os.makedirs(path, exist_ok=True)
    return path

# ============================================================
# Per-run panels (each takes an existing ax + a result dict)
# ============================================================
def plot_NA(ax, res, color=antigen_color, ls='-', label=None):
    ax.plot(res['t'], res['N_A'], lw=6, color=color, ls=ls, label=label)


def plot_pi(ax, res, color, label=None, rows=(0, 10, -10)):
    """Plot pi(t) for a few representative clones (fading alpha)."""
    if 'pi' not in res:
        return None
    t = res['t']
    line = ax.plot(t, res['pi'][rows[0], :], lw=6, color=color, label=label)
    c = line[0].get_color()
    for r, (lw, a) in zip(rows[1:], [(6, 0.5), (4, 0.2)]):
        ax.plot(t, res['pi'][r, :], lw=lw, alpha=a, color=c)
    return c


def plot_NB(ax, res, color, label=None, rows=(0, 10, -10)):
    """Per-clone total B cells N_Bo+N_Ba for representative clones."""
    t = res['t']
    N_B = res['N_Bo'] + res['N_Ba'] if 'N_Ba' in res else res['N_B']
    line = ax.plot(t, N_B[rows[0], :], lw=6, color=color, label=label)
    c = line[0].get_color()
    for r, (lw, a) in zip(rows[1:], [(6, 0.5), (4, 0.2)]):
        ax.plot(t, N_B[r, :], lw=lw, alpha=a, color=c)
    return c


def plot_yield(ax, res, color, label=None):
    ax.plot(res['t'], compute_N_B_tot(res), lw=4, color=color, label=label)


def plot_potency(ax, res, color, ls='-', label=None, mark_Zc=True, mark_t = False,
                 marker='o', marker_face=None, marker_edge='k'):
    """Potency Z(t); optionally mark the first crossing of Z_c."""
    t = res['t']
    Z = compute_potency_t(res)
    ax.plot(t[Z > 0], Z[Z > 0], lw=5, color=color, ls=ls)
    if mark_Zc:
        Z_c = res['params'].Z_c
        hit = np.where(Z >= Z_c)[0]
        if len(hit):
            ax.scatter(t[hit[0]], Z_c, marker = marker,
                       facecolor=color if marker_face is None else marker_face,
                       edgecolor=marker_edge, linewidth=2, s=250, zorder=100, label=label)
    elif mark_t:
        hit = np.where(t >= mark_t)[0]
        if len(hit):
            ax.scatter(mark_t, Z[hit[0]], marker = marker,
                        facecolor=color if marker_face is None else marker_face,
                        edgecolor=marker_edge, linewidth=2, s=180, zorder=100, label=label)
    return Z


def plot_innate(ax, res, color, ls='--', label=None):
    """External innate-suppression fraction (memory-off runs)."""
    p, t = res['params'], res['t']
    innate_proxy = np.exp(p.lambda_innate * t)
    innate = (1.0 + p.threshold_innate / (innate_proxy + 1e-30)) ** (-1)
    ax.plot(t, innate, lw=4, color=color, ls=ls, label=label)


def plot_pb_from_potency(ax, res, color, ls='-', label=None):
    """Neutralization probability pb(t) computed from the potency Z(t)."""
    Z = compute_potency_t(res)
    pb = (1.0 + res['params'].Z_c / (Z + 1e-30)) ** (-1)
    ax.plot(res['t'], pb, lw=4, color=color, ls=ls, label=label)


# ============================================================
# Analytic affinity-front overlay (was duplicated in two scripts)
# ============================================================
def compute_affinity_front(res, on_slope=0.85):
    """
    Reconstruct the on/off affinity front DG(t) from the simulation timing.

    Returns a dict with (ton, DG_on), (toff, DG_off), (toff, DG_off_null), and
    the key times tstar/tpeak. `on_slope` scales the growth-phase slope
    (0.78 in the dynamics scripts, 1.0 in the h0 script).
    """
    t, p = res['t'], res['params']
    N_B = res['N_Bo'] + res['N_Ba'] if 'N_Ba' in res else res['N_B']
    best = np.argmax(N_B[:, 0] > 0)
    bump = 0.06 if not p.memory else 0.01
    tstar = t[N_B[best, :] > N_B[best, 0] + bump * N_B[best, 0]][0]

    tpeak_NA = t[res['N_A'] == np.max(res['N_A'])][0]
    if 'pi' in res:
        tpeak_pi = t[res['pi'][0, :] == np.max(res['pi'][0, :])][0]
        if p.memory==0:
            tpeak = tpeak_NA # + (tpeak_pi - tpeak_NA) / 2
            # tpeak = tpeak_pi
        else:
            tpeak = tpeak_NA
            # tpeak = tpeak_pi
    else:
        tpeak = tpeak_NA

    ton = t[(t < tpeak) & (t > tstar)]
    toff = t[t >= tpeak]
    DGpeak = on_slope * p.lambda_A / p.eta * (tpeak - tstar)
    DG_on = on_slope * p.lambda_A / p.eta * (ton - tstar)
    DG_off = p.b0 / p.eta * (tpeak - toff) + DGpeak
    DG_off_null = p.delta_A / p.eta * (tpeak - toff) + DGpeak
    return dict(tstar=tstar, tpeak=tpeak, ton=ton, toff=toff,
                DG_on=DG_on, DG_off=DG_off, DG_off_null=DG_off_null, DGpeak=DGpeak)


# ============================================================
# Multi-panel diagnostics (approx model: needs D / N_T_free / lambda_B)
# ============================================================
def plot_diagnostics(res):
    """Six-panel overview for the 'approx' model."""
    t, p, DG = res['t'], res['params'], res['DG']
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    ax = axes[0, 0]
    ax.semilogy(t, res['N_A'])
    ax.set_xlabel('Time'); ax.set_ylabel('$N_A$'); ax.set_title('Antigen')

    ax = axes[0, 1]
    ax.plot(t, res['N_T_free'], label='$N_T^{free}$')
    ax.plot(t, res['N_T'], '--', alpha=0.5, label='$N_T$')
    ax.set_xlabel('Time'); ax.set_ylabel('T cells'); ax.set_title('T cell pool')
    ax.legend()

    ax = axes[0, 2]
    ax.semilogy(t, res['D'] + 1e-10)
    ax.axhline(1, color='k', linestyle=':', alpha=0.5)
    ax.set_xlabel('Time'); ax.set_ylabel('$D(t)$'); ax.set_title('Demand')

    ax = axes[1, 0]
    ax.semilogy(t, compute_N_B_tot(res) + 0.1, label=r'$\bar{N}$')
    ax.semilogy(t, compute_L_act(res) + 0.1, label='$L_{act}$')
    ax.set_xlabel('Time'); ax.set_ylabel('Count'); ax.set_title('B cell response')
    ax.legend()

    ax = axes[1, 1]
    N_B_final = res['N_B'][:, -1]
    if res['mode'] == 'stochastic':
        ax.scatter(DG, N_B_final, s=1, alpha=0.3)
    else:
        ax.semilogy(DG, N_B_final)
    ax.set_yscale('log')
    ax.set_xlabel('$\\Delta G$'); ax.set_ylabel('$N_B$')
    ax.set_title('Clone sizes (final)')

    ax = axes[1, 2]
    ranks, sizes = compute_zipf(res, time_index=-1)
    if len(ranks):
        ax.loglog(ranks, sizes, '.', markersize=2)
        ax.set_xlabel('Rank $k$'); ax.set_ylabel('$N_k$')
        ax.set_title('Zipf plot (final)')

    plt.tight_layout()
    return fig
