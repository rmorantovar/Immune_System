"""Zipf analysis of B-cell clone sizes in Mesin et al. 2020 (Cell 180:92), Table S1.

For each group of mice (experiment x phenotype, see GROUPS):
  - rank curve of clone sizes per mouse and averaged over mice,
  - Zipf exponent zeta of the averaged curve, bootstrapped over mice, and a rate-model alternative,
  - per-mouse quantities: N1 (largest clone), L_act (number of clones), N (cells), entropies.
Then figures comparing groups: scaling relations, entropies and pseudo-energy distributions.
Clones are defined by CLONE_COLS.
"""
import os
import sys
import pickle
from functools import lru_cache
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append('../../../library/')
from funcs import *

# ── Paths & constants ─────────────────────────────────────────────────────────
project    = 'exponential_proofreading'
subproject = 'data'
root_dir   = f'{DRIVE}/Research_files/Immune_system/{project}/{subproject}/mesin2020'
output_plot = (f'{DRIVE}/_Documents/Science/research_projects/'
               f'Immune_System/figures/{project}/{subproject}/mesin2020/zipf/all')
os.makedirs(output_plot, exist_ok=True)

EXCEL_FILE   = root_dir + '/1-s2.0-S0092867419313170-mmc1.xlsx'
CLONE_COLS   = ['V', 'J', 'D']     # or ['CDR3:']
MAX_RANK     = 100                 # ranks kept in the averaged rank curve
MAX_RANK_FIT = 20                  # ranks used to fit zeta
N_ENSEMBLE   = 400                 # bootstrap resamples over mice
ALPHA        = 1.3                 # pseudo-energy scale: E = -ln(n/N) / ALPHA
EBETA_STAR   = 23                  # sets the reference entropy line S_STAR
S_STAR       = np.sqrt(2 / np.pi * EBETA_STAR) - np.log(2)

# Experiments to run (keys of GROUPS): '1', '2', '3a', '3b', '2-3a', '2-3b', '4a', '4b', '4c', '4d'
RUN_EXPERIMENTS = {'1', '2-3a', '4a', '4c'}

rng = np.random.default_rng(0)

fig_kw    = dict(figsize=(8 * 1.62, 8), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
fig_kw_sq = dict(figsize=(8, 8),        gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})


# ── Groups ────────────────────────────────────────────────────────────────────
# Each group pools mice from one or more sources (sheet + row filters; a list value keeps any of its values).
# Experiment / Response / phenotype label the group in all outputs.

PHOTO = dict(sheet='Photoactivation CGG', header=1, mouse_col='Mouse')
FATE  = dict(sheet='Fate-mapping CGG',    header=1, mouse_col='Mouse')
FLU   = dict(sheet='Influenza_IGH',       header=2, mouse_col='Experiment / Mouse')
# Authors' influenza sheets: keep location and the original sort labels lost in Influenza_IGH
FLU_CONFETTI = dict(sheet='AID-Confetti Influenza', header=2, mouse_col='Experiment / Mouse')
FLU_S1PR2    = dict(sheet='S1pr2-Tomato Influenza', header=2, mouse_col='Experiment / Mouse')

def src(base, **filters):
    return dict(base, filters=filters)

GROUPS = {
    # naive GC (Fig. 1D)
    '1':    [dict(sources=[src(PHOTO, Figure=1)],
                  Experiment='CGG-alum', Response='primary', phenotype='GC', color=my_red, marker='*', ms=18)],
    # recall GC + fm (Fig. 4A)
    '2':    [dict(sources=[src(FATE, Figure='4A')],
                  Experiment='2', Response='recall', phenotype='GC + fm', color=my_blue, marker='o')],
    # recall per phenotype (Fig. 4C)
    '3a':   [dict(sources=[src(FATE, Figure='4C-H', Phenotype=ph)],
                  Experiment='3', Response='recall', phenotype=ph, color=c, marker='^')
             for ph, c in [('GC + fm', my_blue), ('PB + fm', my_cyan)]],
    # recall, phenotypes combined (Fig. 4C)
    '3b':   [dict(sources=[src(FATE, Figure='4C-H')],
                  Experiment='3', Response='recall', phenotype='combined', color=my_purple, marker='v')],
    # recall GC + fm, pooled Fig. 4A and 4C
    '2-3a': [dict(sources=[src(FATE, Figure='4A'), src(FATE, Figure='4C-H', Phenotype='GC + fm')],
                  Experiment='CGG-alum', Response='recall', phenotype='GC + fm', color=my_blue, marker='s')],
    # recall all phenotypes, pooled Fig. 4A and 4C
    '2-3b': [dict(sources=[src(FATE, Figure='4A'), src(FATE, Figure='4C-H')],
                  Experiment='2-3b', Response='recall', phenotype='GC + pb + fm', color=my_blue2, marker='P')],
    # influenza infection per sort (Fig. 5); other sorts: 'PB + fm', 'M', 'PC'.
    # 'GC' is the total GC of the boosted popliteal LN; 'GC + fm' pools fate-mapped GC cells of the
    # boosted popliteal LN (recall) with HA+ fate-mapped GC cells of the primary mediastinal LN (~2/3).
    '4a':   [dict(sources=[src(FLU, Sort=sort)],
                  Experiment='HA', Response=resp, phenotype=sort, color=c, marker='D')
             for sort, resp, c in [('GC', 'primary', my_red), ('GC + fm', 'recall', my_blue)]],
    # influenza recall, fate-mapped GC cells of the boosted popliteal LN only (day 6 mice; excludes Fig. S5, day 9)
    '4c':   [dict(sources=[src(FLU_CONFETTI, Sort='GC fate-mapped', Figure=[5, 'S4']),
                           src(FLU_S1PR2,    Sort='GC fate-mapped')],
                  Experiment='HA-pLN', Response='recall', phenotype='GC + fm', color=my_blue, marker='h')],
    # influenza recall, fate-mapped plasmablasts of the boosted popliteal LN, alone and pooled with 4c
    '4d':   [dict(sources=[src(FLU_CONFETTI, Sort=sorts, Figure=[5, 'S4']), src(FLU_S1PR2, Sort=sorts)],
                  Experiment='HA-pLN', Response='recall', phenotype=ph, color=c, marker='h')
             for sorts, ph, c in [('PB fate-mapped', 'PB + fm', my_cyan),
                                  (['GC fate-mapped', 'PB fate-mapped'], 'GC + PB + fm', my_purple)]],
    # influenza infection, all fate-mapped cells (Fig. 5)
    '4b':   [dict(sources=[src(FLU, Sort2='fm')],
                  Experiment='4b', Response='recall', phenotype='combined', color=my_blue2, marker='X')],
}


# ── Data ──────────────────────────────────────────────────────────────────────

@lru_cache(maxsize=None)
def read_sheet(sheet, header):
    return pd.read_excel(EXCEL_FILE, sheet_name=sheet, header=header)

def load_counts(sources):
    """Clone sizes, one array per mouse, pooling the mice of all sources."""
    counts = []
    for s in sources:
        data = read_sheet(s['sheet'], s['header'])
        for col, val in s['filters'].items():
            data = data[data[col].isin(val if isinstance(val, list) else [val])]
        sizes = data.groupby([s['mouse_col']] + CLONE_COLS).size()
        counts += [g.to_numpy() for _, g in sizes.groupby(level=0, sort=False)]
    return counts


# ── Models & fits ─────────────────────────────────────────────────────────────

def zipf_logx(logj, z):
    """Zipf: log(N_j/N_1) = -zeta log j."""
    return -z * logj

def rate_logx(logj, z0, c):
    """Rate model: log(N_j/N_1) = log N_1 (j^{-c} - 1), with z0 = c log N_1.
    Reduces to Zipf (-z0 log j) for c -> 0; z0 j^{-c} is the local exponent."""
    return z0 * np.expm1(-c * logj) / c

def fit_both(logj, logx):
    """Fit both models to a normalized rank curve (rank 1 excluded: exact for both)."""
    logj, logx = logj[1:], logx[1:]
    n = len(logj)
    pz, _ = curve_fit(zipf_logx, logj, logx, p0=[0.5])
    rss_z = np.sum((logx - zipf_logx(logj, *pz)) ** 2)
    pr, _ = curve_fit(rate_logx, logj, logx, p0=[pz[0], 0.1], bounds=([0, 1e-6], [10, 3]))
    rss_r = np.sum((logx - rate_logx(logj, *pr)) ** 2)
    aic_z = n * np.log(rss_z / n) + 2 * 1
    aic_r = n * np.log(rss_r / n) + 2 * 2
    return dict(zeta=pz[0], z0=pr[0], c=pr[1], logN1=pr[0] / pr[1], dAIC=aic_r - aic_z)

def fit_rank_curve(x):
    """fit_both on the first MAX_RANK_FIT ranks of a normalized rank curve x_k = N_k/N_1."""
    logj = np.log(np.arange(1, len(x) + 1))[:MAX_RANK_FIT]
    return fit_both(logj, np.log(x)[:MAX_RANK_FIT])

def mean_rank_curve(xs):
    """Average of x_k = N_k/N_1 over mice, for the ranks present in at least 2 mice."""
    n_mice, total = np.zeros(MAX_RANK), np.zeros(MAX_RANK)
    for x in xs:
        x = x[:MAX_RANK]
        n_mice[:len(x)] += 1
        total[:len(x)]  += x / x[0]
    k_max = np.sum(n_mice > 1)
    return total[:k_max] / n_mice[:k_max]

def bootstrap_rank_fits(xs):
    """Rank-curve fits on all mice (full) and on N_ENSEMBLE resamples of the mice (boot)."""
    full = fit_rank_curve(mean_rank_curve(xs))
    boot = [fit_rank_curve(mean_rank_curve([xs[i] for i in rng.integers(len(xs), size=len(xs))]))
            for _ in tqdm(range(N_ENSEMBLE))]
    return full, pd.DataFrame(boot)

def barN(L, Z):
    """Total cells of a continuous Zipf with N1 = L^Z and L clones."""
    return L ** Z / (1 - Z) * (L ** (1 - Z) - 1)

def zipf_p(L, Z):
    p = np.arange(1, int(L) + 1) ** (-Z)
    return p / np.sum(p)


# ── Entropies ─────────────────────────────────────────────────────────────────

def plugin_entropy(c):
    p = c[c > 0] / np.sum(c)
    return -np.sum(p * np.log(p))

def KL_zipf(L, Z):
    """Exact ln(L) - S for a finite Zipf distribution p_k ~ k^-Z, k = 1..L."""
    p = zipf_p(L, Z)
    return np.log(L) + np.sum(p * np.log(p))

def rarefied_entropy(counts, n_grid, n_sub=50):
    """Mean plug-in entropy of n cells subsampled without replacement, for each n in n_grid (nan if n > N)."""
    cells = np.repeat(np.arange(len(counts)), counts)
    valid = n_grid <= len(cells)
    S = np.zeros(np.sum(valid))
    for _ in range(n_sub):
        perm = rng.permutation(cells)
        S += [plugin_entropy(np.bincount(perm[:n])) for n in n_grid[valid]]
    out = np.full(len(n_grid), np.nan)
    out[valid] = S / n_sub
    return out

def zipf_rarefied_entropy(Z, L, n_grid, n_sub=50):
    """Same as rarefied_entropy, for n cells drawn from a Zipf distribution p_k ~ k^-Z, k = 1..L."""
    p = zipf_p(L, Z)
    S = np.zeros(len(n_grid))
    for _ in range(n_sub):
        draws = rng.choice(len(p), np.max(n_grid), p=p)
        S += [plugin_entropy(np.bincount(draws[:n])) for n in n_grid]
    return S / n_sub


# ── Per-mouse quantities ──────────────────────────────────────────────────────

def mouse_stats(counts):
    N = np.sum(counts)
    L = len(counts)
    x = np.sort(counts)[::-1]
    n_fit = min(L, MAX_RANK_FIT)
    (zeta,), _ = curve_fit(zipf_logx, np.log(np.arange(1, n_fit + 1)), np.log(x[:n_fit] / x[0]))
    S = plugin_entropy(counts)
    cell_freqs = np.repeat(counts / N, (1000 * counts / N).astype(int))
    return dict(
        N1=x[0], L_act=L, barN=N,
        S=S, S_corr=S + (L - 1) / (2 * N),                                       # Miller-Madow correction
        lambda2=np.sum(counts * (counts - 1)) / (N * (N - 1)) if N > 1 else np.nan,  # unbiased Simpson
        zeta=zeta,
        counts=counts,
        N=-np.log(counts / N) / ALPHA,           # pseudo-energy per clone (clone-weighted)
        N_cells=-np.log(cell_freqs) / ALPHA,     # pseudo-energy per cell (cell-weighted, ~1000 per mouse)
    )


# ── Plot helpers ──────────────────────────────────────────────────────────────

def save(fig, name):
    fig.savefig(f'{output_plot}/{name}.pdf', dpi=800, transparent=True)

def plot_mouse_top_clones(x, color, name):
    fig, ax = plt.subplots(**fig_kw_sq)
    ax.bar(range(1, len(x[:10]) + 1), x[:10], color=color)
    ax.set_xticks([])
    ax.set_xlim(0.5, 10.5)
    ax.set_yscale('log')
    ax.set_ylim(bottom=0.9, top=50)
    ax.tick_params(axis='both', labelsize=50)
    save(fig, name)
    plt.close(fig)

def plot_zeta_violin(ax, values, values_mice, position, color):
    parts = ax.violinplot([values], positions=[position], showmeans=True, showextrema=False)
    for body in parts['bodies']:
        body.set_facecolor(color)
        body.set_edgecolor('black')
        body.set_alpha(0.5)
    parts['cmeans'].set_color(color)
    ax.scatter(rng.normal(position, 0.04, len(values_mice)), values_mice, color=color, edgecolor='k', s=80, alpha=.5)
    ax.scatter(position, np.mean(values_mice), color=color, edgecolor='k', s=150)

def layout_ranking(ax):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(left=0.9, right=5e1)
    ax.set_ylim(bottom=2e-2, top=1.1)
    ax.tick_params(axis='both', labelsize=30)

def layout_zeta(ax, labels):
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylabel(r'$\zeta$', fontsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.tick_params(axis='x', labelsize=20)


# ── Per-group analysis: rank curves and zeta ──────────────────────────────────

fig_r,    ax_r    = plt.subplots(**fig_kw)
fig_zeta, ax_zeta = plt.subplots(**fig_kw)
layout_ranking(ax_r)
save(fig_r, 'ranking_B_cells_0')

rows       = []   # one per mouse and group
group_fits = []   # one per group
for key, groups in GROUPS.items():
    if key not in RUN_EXPERIMENTS:
        continue
    for g in groups:
        labels = dict(Experiment=g['Experiment'], Response=g['Response'], phenotype=g['phenotype'])
        counts = load_counts(g['sources'])
        xs     = [np.sort(c)[::-1] for c in counts]
        print(f"{key}: {g['Experiment']} / {g['phenotype']}, {len(counts)} mice")

        stats = [mouse_stats(c) for c in counts]
        rows += [dict(labels, **s) for s in stats]
        zetas_mice = np.array([s['zeta'] for s in stats])

        full, boot = bootstrap_rank_fits(xs)
        zeta_mean, zeta_std = boot['zeta'].mean(), boot['zeta'].std()
        print(f'  zeta = {full["zeta"]:.3f}, bootstrap {zeta_mean:.3f} +- {zeta_std:.3f}')
        group_fits.append(dict(
            run=key, **labels, n_mice=len(counts), color=g['color'],
            zeta_hat=full['zeta'], zeta_mean=zeta_mean, zeta_std=zeta_std,
            c_hat=full['c'], c_std=boot['c'].std(), c_hi=boot['c'].quantile(.975),
            frac_c_at_zero=(boot['c'] < 1e-3).mean(),
            dAIC_hat=full['dAIC'], dAIC_lo=boot['dAIC'].quantile(.025), dAIC_hi=boot['dAIC'].quantile(.975),
            mice_zeta_mean=zetas_mice.mean(), mice_zeta_std=zetas_mice.std(),
        ))

        # Rank curves: each mouse, mouse average, Zipf fit and rate-model fit
        color = g['color']
        for x in xs:
            ax_r.step(range(1, len(x) + 1), x / x[0], color=color, alpha=.5, lw=0.5)
        x_avg = mean_rank_curve(xs)
        j = np.arange(1, len(x_avg) + 10)
        ax_r.plot(range(1, len(x_avg) + 1), x_avg, color=color, markerfacecolor='None',
                  ms=g.get('ms', 12), ls='', marker=g['marker'], label=f"${zeta_mean:.2f}$ ; {g['phenotype']}")
        ax_r.plot(j, j ** (-zeta_mean), color=color, alpha=.8, lw=3)
        ax_r.plot(j, np.exp(rate_logx(np.log(j), full['z0'], full['c'])), color=color, lw=2, ls=':')

        plot_zeta_violin(ax_zeta, boot['zeta'], zetas_mice, position=len(group_fits) - 1, color=color)

        for i_m, x in enumerate(xs):
            plot_mouse_top_clones(x, color, f"entropy_mouse_{g['Experiment']}_{g['phenotype']}_{i_m + 1}")

    # Snapshot of the cumulative figures after each experiment
    suffix = key.replace('-', '')
    layout_ranking(ax_r)
    save(fig_r, f'ranking_B_cells_{suffix}')
    layout_zeta(ax_zeta, [f"{gf['Experiment']}\n{gf['phenotype']}" for gf in group_fits])
    save(fig_zeta, f'zetas_{suffix}')

scaling_results = pd.DataFrame(rows)
group_fits      = pd.DataFrame(group_fits)
group_zeta      = {(gf.Experiment, gf.phenotype): gf.zeta_hat for gf in group_fits.itertuples()}

hue_colors = dict(zip(group_fits['phenotype'], group_fits['color']))
hue_colors = {ph: hue_colors[ph] for ph in scaling_results['phenotype'].unique()}   # first group's colour per phenotype
exp_ls     = dict(zip(scaling_results['Experiment'].unique(), ['-', '--', ':', '-.']))

def scatter(ax, x, y, hue='phenotype', palette=None):
    sns.scatterplot(data=scaling_results, x=x, y=y, hue=hue, style='Experiment', ax=ax,
                    s=150, palette=palette or hue_colors, edgecolors='black', alpha=0.8)


# ── Scaling relations ─────────────────────────────────────────────────────────

# N1 vs N
fig, ax = plt.subplots(**fig_kw)
scatter(ax, 'barN', 'N1')
ax.set_xlim(left=30, right=200)
ax.tick_params(axis='both', labelsize=30)
ax.legend(title_fontsize=20, fontsize=20, loc=2)
save(fig, 'scaling_N1_Ntot')

# L_act vs N, with the continuous-Zipf prediction (band: zeta +- 2 std)
fig, ax = plt.subplots(**fig_kw)
L = np.linspace(1, 400, 100)
for gf in group_fits.itertuples():
    ax.plot(barN(L, gf.zeta_mean), L, color=gf.color, linestyle='--')
    ax.fill_betweenx(L, barN(L, gf.zeta_mean - 2 * gf.zeta_std), barN(L, gf.zeta_mean + 2 * gf.zeta_std),
                     color=gf.color, alpha=0.1)
scatter(ax, 'barN', 'L_act')
ax.set_xlim(left=30, right=200)
ax.set_ylim(bottom=0, top=125)
ax.tick_params(axis='both', labelsize=30)
ax.legend(title_fontsize=20, fontsize=20, loc=2)
save(fig, 'scaling_Lact_Ntot')

# L_act vs N1, with L = N1^(1/zeta); shading: above / below L = N1
fig, ax = plt.subplots(**fig_kw)
N1 = np.linspace(0, 70, 100)
ax.fill_between(N1, 0, N1, color=my_blue, alpha=0.2)
ax.fill_between(N1, N1, 200, color=my_red, alpha=0.2)
for gf in group_fits.itertuples():
    ax.plot(N1, N1 ** (1 / gf.zeta_mean), color=gf.color, linestyle='--')
scatter(ax, 'N1', 'L_act', hue='Response', palette=[my_red, my_blue, my_cyan, my_purple][:scaling_results['Response'].nunique()])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xlim(left=0, right=60)
ax.set_ylim(bottom=0, top=125)
ax.tick_params(axis='both', labelsize=30)
ax.legend(title_fontsize=20, fontsize=20, loc=0)
save(fig, 'scaling_Lact_N1')

# N vs its prediction from N1, L_act and the mouse's zeta
scaling_results['barN_prediction'] = [r.N1 * np.sum(np.arange(1, r.L_act + 1) ** (-r.zeta))
                                      for r in scaling_results.itertuples()]
fig, ax = plt.subplots(**fig_kw)
nbar = np.linspace(30, 300, 100)
ax.plot(nbar, nbar, color='k', linestyle='--')
scatter(ax, 'barN_prediction', 'barN')
my_plot_layout(ax=ax, yscale='log', xscale='log', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax.tick_params(axis='both', labelsize=30)
ax.legend(title_fontsize=20, fontsize=20, loc=4)
save(fig, 'Ntot_prediction')


# ── Entropy ───────────────────────────────────────────────────────────────────
x = np.linspace(1, 500, 1000)

# S vs N
fig, ax = plt.subplots(**fig_kw)
scatter(ax, 'barN', 'S_corr')
ax.plot(x, np.full_like(x, S_STAR), color=my_blue, linestyle='--')
my_plot_layout(ax=ax, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax.set_xlabel(r'$N_{\mathrm{sample}}$', fontsize=30)
ax.set_ylabel(r'$S$', fontsize=30)
ax.set_xlim(left=2, right=200)
ax.tick_params(axis='both', labelsize=30)
save(fig, 'entropy_Ntot')

# S vs L_act, with the maximum entropy S = ln(L_act) (all clones equal size)
fig, ax = plt.subplots(**fig_kw)
scatter(ax, 'L_act', 'S_corr')
ax.plot(x, np.full_like(x, S_STAR), color=my_blue, linestyle='--')
ax.plot(x, np.log(x), color='k', linestyle='--')
my_plot_layout(ax=ax, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax.set_xlim(left=1, right=120)
ax.tick_params(axis='both', labelsize=30)
ax.legend(title_fontsize=20, fontsize=20, loc=0)
save(fig, 'entropy_Lact')

# Entropy deficit ln(L_act) - S (KL divergence from uniform) vs L_act,
# with the finite-Zipf prediction for each group's zeta
scaling_results['SlogLact_corr'] = np.log(scaling_results['L_act']) - scaling_results['S_corr']
fig, ax = plt.subplots(**fig_kw)
scatter(ax, 'L_act', 'SlogLact_corr')
L_th = np.arange(2, 121)
for (exp, ph), z in group_zeta.items():
    ax.plot(L_th, [KL_zipf(L, z) for L in L_th], color=hue_colors[ph], linestyle=exp_ls[exp], lw=2)
my_plot_layout(ax=ax, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax.set_xlabel(r'$L_{\mathrm{act}}$', fontsize=30)
ax.set_ylabel(r'$\ln L_{\mathrm{act}} - S$', fontsize=30)
ax.set_xlim(left=1, right=120)
ax.set_ylim(bottom=-0.5, top=1.5)
ax.tick_params(axis='both', labelsize=30)
ax.legend(title_fontsize=20, fontsize=20, loc=2)
save(fig, 'entropy_deficit_Lact')

# Rarefaction: mean entropy of n subsampled cells per group (removes sequencing-depth differences).
# Lines: Zipf prediction with the group's zeta and each mouse's L_act, averaged over the same mice.
MIN_MICE_RARE = 3   # drop n where fewer mice than this have >= n cells
n_grid = np.unique(np.round(np.geomspace(2, scaling_results['barN'].max(), 30)).astype(int))
groups = list(scaling_results.groupby(['Experiment', 'phenotype'], sort=False))

fig, ax = plt.subplots(**fig_kw)
for (exp, ph), grp in groups:
    color  = hue_colors[ph]
    S_data = np.array([rarefied_entropy(c, n_grid) for c in grp['counts']])
    n_mice = np.sum(np.isfinite(S_data), axis=0)
    ok     = n_mice >= MIN_MICE_RARE
    mean   = np.nanmean(S_data[:, ok], axis=0)
    sem    = np.nanstd(S_data[:, ok], axis=0) / np.sqrt(n_mice[ok])
    ax.plot(n_grid[ok], mean, color=color, ls='', marker='o', ms=8, label=f'{exp} ; {ph}')
    ax.fill_between(n_grid[ok], mean - sem, mean + sem, color=color, alpha=0.2)
    S_zipf = np.array([zipf_rarefied_entropy(group_zeta[(exp, ph)], L, n_grid) for L in grp['L_act']])
    S_zipf[~np.isfinite(S_data)] = np.nan   # same mice at each n as the data
    ax.plot(n_grid[ok], np.nanmean(S_zipf[:, ok], axis=0), color=color, ls=exp_ls[exp], lw=2)
my_plot_layout(ax=ax, yscale='linear', xscale='log', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax.set_xlabel(r'$n$ (cells subsampled)', fontsize=30)
ax.set_ylabel(r'$\langle S(n) \rangle$', fontsize=30)
ax.tick_params(axis='both', labelsize=30)
ax.legend(title_fontsize=20, fontsize=20, loc=2)
save(fig, 'entropy_rarefaction')

# Renyi-2 entropy S2 = -ln(lambda), lambda = prob. two cells of the same mouse share a clone.
# Violin: population S2 = -ln(mean lambda over mice), bootstrapped over mice; dots: per-mouse S2.
# Black bar: Zipf prediction with the group's zeta and each mouse's L_act.
fig, ax = plt.subplots(**fig_kw)
for i, ((exp, ph), grp) in enumerate(groups):
    lam   = grp['lambda2'].dropna().to_numpy()
    boots = [-np.log(np.mean(rng.choice(lam, len(lam)))) for _ in range(N_ENSEMBLE)]
    plot_zeta_violin(ax, boots, -np.log(lam[lam > 0]), position=i, color=hue_colors[ph])
    lam_zipf = [np.sum(zipf_p(L, group_zeta[(exp, ph)]) ** 2) for L in grp['L_act']]
    ax.scatter(i, -np.log(np.mean(lam_zipf)), color='k', marker='_', s=1500, lw=3)
ax.set_xticks(range(len(groups)), [f'{exp}\n{ph}' for (exp, ph), _ in groups])
ax.set_ylabel(r'$S_2$', fontsize=30)
ax.tick_params(axis='y', labelsize=30)
ax.tick_params(axis='x', labelsize=20)
save(fig, 'entropy_S2')


# ── Binding-energy densities (data version of manuscript Fig. potency A/B) ────
# EP: n ∝ exp(-alpha_hat dG), so dG = ln(N1/n) / alpha_hat in each mouse (best clone at dG = 0),
# with alpha_hat = ALPHA (primary) or 2 ALPHA (recall). Pooled over the mice of each response:
#   Omega(dG): clones per mouse per unit dG (clone-weighted); predicted log-slope beta* = alpha_hat / zeta
#   q(dG):     pdf of cells, each mouse weighted equally;      predicted log-slope beta* - alpha_hat
# Fits skip the dG = 0 bin (each mouse's best clone, there by construction, as rank 1 in fit_both)
# and stop at the median doubleton energy ln(N1/2) / alpha_hat: beyond it the density is cut by sampling.
ALPHA_HAT   = {'primary': ALPHA, 'recall': 2 * ALPHA}
DG_BIN      = 0.2
resp_colors = {'primary': my_red, 'recall': my_blue}

def fit_log_density(centers, dens, dG_max):
    sel = (dens > 0) & (centers > DG_BIN / 2) & (centers <= dG_max)
    slope, intercept = np.polyfit(centers[sel], np.log(dens[sel]), 1)
    return slope, intercept

fig_clones, ax_clones = plt.subplots(**fig_kw)   # clone-weighted only: same slope beta* predicted for both responses
for resp, grp in scaling_results.groupby('Response', sort=False):
    a_hat   = ALPHA_HAT[resp]
    color   = resp_colors[resp]
    zeta    = group_fits.loc[group_fits['Response'] == resp, 'zeta_hat'].mean()
    beta_pr = a_hat / zeta                       # predicted beta*
    dG      = [np.log(r.N1 / r.counts) / a_hat for r in grp.itertuples()]
    w_cells = [r.counts / r.barN for r in grp.itertuples()]
    dG_fit  = np.median([np.log(r.N1 / 2) / a_hat for r in grp.itertuples()])
    bins    = np.arange(0, np.max(np.concatenate(dG)) + 2 * DG_BIN, DG_BIN) - DG_BIN / 2
    centers = (bins[:-1] + bins[1:]) / 2
    Omega   = np.histogram(np.concatenate(dG), bins)[0] / (len(grp) * DG_BIN)
    q       = np.histogram(np.concatenate(dG), bins, weights=np.concatenate(w_cells))[0] / (len(grp) * DG_BIN)

    # Entropy of q and the check S_c = <ln Omega>_q + S_dG[q] (manuscript Eq. S_c_S_E), on the binned densities
    P, K   = q * DG_BIN, Omega * DG_BIN         # cell fraction and clones per mouse in each bin
    S_dG   = -np.sum(P[P > 0] * np.log(q[P > 0]))
    S_c_dG = np.sum(P[P > 0] * np.log(K[P > 0] / P[P > 0]))
    print(f'{resp}: beta* predicted {beta_pr:.2f}; S_dG[q] = {S_dG:.2f}; '
          f'<ln Omega>_q + S_dG = {S_c_dG:.2f} vs mean clonal entropy S_c = {grp["S"].mean():.2f}')

    fig, ax = plt.subplots(**fig_kw)
    for dens, marker, mfc, slope_pr, name in [(Omega, 'o', 'None', beta_pr, r'\Omega'),
                                              (q, 's', color, beta_pr - a_hat, 'q')]:
        slope, intercept = fit_log_density(centers, dens, dG_fit)
        x_fit = np.linspace(0, dG_fit, 10)
        ax.plot(centers[dens > 0], dens[dens > 0], ls='', marker=marker, ms=10, color=color, mfc=mfc,
                label=fr'${name}$: slope ${slope:.2f}$ (pred. ${slope_pr:.2f}$)')
        ax.plot(x_fit, np.exp(intercept + slope * x_fit), color=color, lw=2)
        x_mid = dG_fit / 2                       # prediction drawn through the fit at mid-range
        ax.plot(x_fit, np.exp(intercept + slope * x_mid + slope_pr * (x_fit - x_mid)), color='k', lw=1.5, ls='--')
        if name == r'\Omega':
            ax_clones.plot(centers[dens > 0], dens[dens > 0], ls='', marker='o', ms=10, color=color,
                           label=fr'{resp}: slope ${slope:.2f}$')
            ax_clones.plot(x_fit, np.exp(intercept + slope * x_fit), color=color, lw=2)
    ax.axvline(dG_fit, color='gray', ls=':')
    my_plot_layout(ax=ax, yscale='log', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax.set_xlabel(r'$\Delta G$', fontsize=30)
    ax.set_ylabel(r'$\Omega(\Delta G),\ q(\Delta G)$', fontsize=30)
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(fontsize=20, loc=0, title=fr'$S_{{\Delta G}}[q] = {S_dG:.2f}$', title_fontsize=20)
    save(fig, f'distribution_dG_{resp}')

my_plot_layout(ax=ax_clones, yscale='log', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax_clones.set_xlabel(r'$\Delta G$', fontsize=30)
ax_clones.set_ylabel(r'$\Omega(\Delta G)$', fontsize=30)
ax_clones.tick_params(axis='both', labelsize=30)
ax_clones.legend(fontsize=20, loc=0)
save(fig_clones, 'distribution_dG_clones')


# ── Save ──────────────────────────────────────────────────────────────────────
with open(root_dir + '/scaling_dict.pkl', 'wb') as f:
    pickle.dump(scaling_results.to_dict('list'), f)
group_fits.drop(columns='color').to_csv(output_plot + '/zeta_stats.csv', index=False)
scaling_results.to_csv(root_dir + '/scaling_results.csv', index=False)
