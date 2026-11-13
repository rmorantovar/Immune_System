"""Zipf analysis of B-cell clone sizes in Mesin et al. 2020 (Cell 180:92), Table S1.

For each group of responses (antigen x response x phenotype, see GROUPS):
  - rank curve of clone sizes per response unit (a mouse, or one lymph node) and averaged over units,
  - Zipf exponent zeta of the averaged curve, bootstrapped over units, and a rate-model alternative,
  - per-unit quantities: N1 (largest clone), L_act (number of clones), N (cells), entropies.
Figures compare groups within each panel (see PANELS): the naive response alone, or the naive response
next to recall GC, recall PB, and recall GC + PB. Clones are defined by CLONE_COLS.
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
from mf_plotting import my_orange, my_blue   # manuscript palette (funcs' my_blue is a different shade)

# ── Paths & constants ─────────────────────────────────────────────────────────
project    = 'exponential_proofreading'
subproject = 'data'
root_dir   = f'{DRIVE}/research_files/Immune_system/{project}/{subproject}/mesin2020'
output_plot = (f'{DRIVE}/documents/Science/research_projects/'
               f'Immune_System/figures/{project}/{subproject}/mesin2020/zipf')
os.makedirs(output_plot, exist_ok=True)

EXCEL_FILE   = root_dir + '/mmc1.xlsx'   # Data S1 as published (the curated 1-s2.0-...-mmc1.xlsx has a corrupted Sequence column)
CLONE_COLS   = ['V', 'J', 'D']     # or ['CDR3:']
MAX_RANK     = 100                 # ranks kept in the averaged rank curve
MAX_RANK_FIT = 20                  # ranks used to fit zeta
N_ENSEMBLE   = 400                 # bootstrap resamples over units
ALPHA        = 1.3                 # pseudo-energy scale: E = -ln(n/N) / ALPHA
EBETA_STAR   = 23                  # sets the reference entropy line S_STAR
S_STAR       = np.sqrt(2 / np.pi * EBETA_STAR) - np.log(2)

rng = np.random.default_rng(0)

fig_kw    = dict(figsize=(8 * 1.62, 8), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
fig_kw_sq = dict(figsize=(8, 8),        gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})


# ── Groups ────────────────────────────────────────────────────────────────────
# Each group pools the response units of one or more sources (sheet + row filters; a list value keeps any
# of its values). A unit is one immune response: rows sharing the source's `unit` columns.

PHOTO = dict(sheet='Photoactivation CGG', header=1, unit=['Mouse'])
FATE  = dict(sheet='Fate-mapping CGG',    header=1, unit=['Figure', 'Mouse'])   # mouse numbers restart per figure
# Influenza: two S1pr2-Tomato mice were boosted in both footpads, so each popliteal LN is a separate unit
FLU_CONFETTI = dict(sheet='AID-Confetti Influenza', header=2, unit=['Experiment / Mouse', 'Location'])
FLU_S1PR2    = dict(sheet='S1pr2-Tomato Influenza', header=2, unit=['Experiment / Mouse', 'Location'])

def src(base, **filters):
    return dict(base, filters=filters)

def flu(sort):
    """Influenza sort of the boosted popliteal LNs, day 6 after the boost (Figs. 5, S4; excludes Fig. S5, day 9)."""
    return [src(FLU_CONFETTI, Sort=sort, Figure=[5, 'S4']), src(FLU_S1PR2, Sort=sort)]

def by_sort(sources, col):
    """Same sources with each sort of a unit as a separate sample (own rank curve), instead of pooling its cells."""
    return [dict(s, unit=s['unit'] + [col]) for s in sources]

GROUPS = {
    # CGG-alum. Naive: primary GC B cells, photoactivated in two GCs per mouse and pooled (Fig. 1).
    # Recall: fate-mapped cells after the boost (Fig. 4A: GC only; Fig. 4C-H: GC and PB of the same mice).
    'CGG naive': dict(sources=[src(PHOTO, Figure=1)],
                      Response='primary', phenotype='GC', color=my_orange),
    'CGG GC':    dict(sources=[src(FATE, Phenotype='GC')],
                      Response='recall', phenotype='GC + fm', color=my_purple),
    'CGG PB':    dict(sources=[src(FATE, Phenotype='PB')],
                      Response='recall', phenotype='PB + fm', color=my_cyan),
    'CGG GC+PB': dict(sources=[src(FATE, Figure='4C-H')],
                      Response='recall', phenotype='GC + PB + fm', color=my_blue),
    # GC and PB samples pooled as separate rank curves (not as cells): the union of 'CGG GC' and 'CGG PB'
    'CGG GC,PB': dict(sources=by_sort([src(FATE)], 'Phenotype'),
                      Response='recall', phenotype='GC, PB + fm', color=my_green),
    # Influenza, HA footpad boost after intranasal infection (Fig. 5). Naive: total GC of the boosted LN
    # (3 of 326 AID-Confetti cells are fate-mapped). Recall: fate-mapped cells of the boosted LN.
    # The primary mediastinal LN ('GC HA+ fate-mapped') is not used.
    'HA naive':  dict(sources=flu('Total GC'),
                      Response='primary', phenotype='GC', color=my_orange),
    'HA GC':     dict(sources=flu('GC fate-mapped'),
                      Response='recall', phenotype='GC + fm', color=my_purple),
    'HA PB':     dict(sources=flu('PB fate-mapped'),
                      Response='recall', phenotype='PB + fm', color=my_cyan),
    'HA GC+PB':  dict(sources=flu(['GC fate-mapped', 'PB fate-mapped']),
                      Response='recall', phenotype='GC + PB + fm', color=my_blue),
    'HA GC,PB':  dict(sources=by_sort(flu(['GC fate-mapped', 'PB fate-mapped']), 'Sort'),
                      Response='recall', phenotype='GC, PB + fm', color=my_green),
}

# Each panel gets its own folder of figures: the naive group alone, next to one recall group, or next to all.
# GC-PB pools the GC and PB cells of each unit; GC-PB_curves pools their rank curves.
PANELS = {f'{exp}_{name}': [f'{exp} naive'] + [f'{exp} {r}' for r in recall]
          for exp in ('CGG', 'HA')
          for name, recall in [('naive', []), ('GC', ['GC']), ('PB', ['PB']), ('GC-PB', ['GC+PB']),
                               ('GC-PB_curves', ['GC,PB']), ('all', ['GC', 'PB', 'GC+PB', 'GC,PB'])]}
# Manuscript: naive and recall (GC + PB cells) of both antigens; colour = response, marker and line style = antigen
PANELS['manuscript']        = ['CGG naive', 'CGG GC+PB', 'HA naive', 'HA GC+PB']
PANELS['manuscript_naive']  = ['CGG naive', 'HA naive']
PANELS['manuscript_recall'] = ['CGG GC+PB', 'HA GC+PB']
RUN_PANELS = list(PANELS)
EXP_LS     = {'CGG': '-', 'HA': '--'}   # line style of model curves, per antigen
EXP_MARKER = {'CGG': 'o', 'HA': 's'}     # marker of data points, per antigen


# ── Data ──────────────────────────────────────────────────────────────────────

@lru_cache(maxsize=None)
def read_sheet(sheet, header):
    return pd.read_excel(EXCEL_FILE, sheet_name=sheet, header=header)

def gene(name):
    """Gene call without allele and functionality ('HJ4*01 F ' -> 'HJ4'): the sheets mix both formats."""
    return name.split('*')[0].strip() if isinstance(name, str) else name

def load_counts(sources):
    """Clone sizes, one array per unit, pooling the units of all sources. Cells missing a CLONE_COLS call are dropped."""
    counts = []
    for s in sources:
        data = read_sheet(s['sheet'], s['header'])
        for col, val in s['filters'].items():
            data = data[data[col].isin(val if isinstance(val, list) else [val])]
        cells = pd.concat([data[s['unit']], data[CLONE_COLS].map(gene)], axis=1)
        sizes = cells.groupby(s['unit'] + CLONE_COLS).size()
        counts += [g.to_numpy() for _, g in sizes.groupby(level=list(range(len(s['unit']))), sort=False)]
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
    """Average of x_k = N_k/N_1 over units, for the ranks present in at least 2 units."""
    n_mice, total = np.zeros(MAX_RANK), np.zeros(MAX_RANK)
    for x in xs:
        x = x[:MAX_RANK]
        n_mice[:len(x)] += 1
        total[:len(x)]  += x / x[0]
    k_max = np.sum(n_mice > 1)
    return total[:k_max] / n_mice[:k_max]

def bootstrap_rank_fits(xs):
    """Rank-curve fits on all units (full) and on N_ENSEMBLE resamples of the units (boot)."""
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


# ── Per-unit quantities ───────────────────────────────────────────────────────

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
        N_cells=-np.log(cell_freqs) / ALPHA,     # pseudo-energy per cell (cell-weighted, ~1000 per unit)
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

def plot_zeta_violin(ax, values, values_mice, position, color, marker='o'):
    parts = ax.violinplot([values], positions=[position], showmeans=True, showextrema=False)
    for body in parts['bodies']:
        body.set_facecolor(color)
        body.set_edgecolor('black')
        body.set_alpha(0.5)
    parts['cmeans'].set_color(color)
    ax.scatter(rng.normal(position, 0.04, len(values_mice)), values_mice, color=color, edgecolor='k', s=80, alpha=.5,
               marker=marker)
    ax.scatter(position, np.mean(values_mice), color=color, edgecolor='k', s=150, marker=marker)

def layout_ranking(ax):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(left=0.9, right=5e1)
    ax.set_ylim(bottom=2e-2, top=1.1)
    ax.tick_params(axis='both', labelsize=30)

def layout_violins(ax, labels, ylabel):
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylabel(ylabel, fontsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.tick_params(axis='x', labelsize=20)


# ── Per-group analysis: rank curves and zeta ──────────────────────────────────

rows       = []   # one per unit and group
group_fits = []   # one per group
rank_fits  = {}   # rank curves and bootstrap fits, per group
os.makedirs(f'{output_plot}/units', exist_ok=True)
for key, g in GROUPS.items():
    if not any(key in PANELS[p] for p in RUN_PANELS):
        continue
    labels = dict(group=key, Experiment=key.split()[0], Response=g['Response'], phenotype=g['phenotype'])
    counts = load_counts(g['sources'])
    xs     = [np.sort(c)[::-1] for c in counts]
    print(f'{key}: {len(counts)} units, cells per unit {[int(np.sum(c)) for c in counts]}')

    stats = [mouse_stats(c) for c in counts]
    rows += [dict(labels, **s) for s in stats]
    zetas_mice = np.array([s['zeta'] for s in stats])

    full, boot = bootstrap_rank_fits(xs)
    zeta_mean, zeta_std = boot['zeta'].mean(), boot['zeta'].std()
    print(f'  zeta = {full["zeta"]:.3f}, bootstrap {zeta_mean:.3f} +- {zeta_std:.3f}')
    group_fits.append(dict(
        **labels, n_mice=len(counts), n_cells=int(sum(np.sum(c) for c in counts)),
        color=g['color'], marker=EXP_MARKER[labels['Experiment']],
        zeta_hat=full['zeta'], zeta_mean=zeta_mean, zeta_std=zeta_std,
        c_hat=full['c'], c_std=boot['c'].std(), c_hi=boot['c'].quantile(.975),
        frac_c_at_zero=(boot['c'] < 1e-3).mean(),
        dAIC_hat=full['dAIC'], dAIC_lo=boot['dAIC'].quantile(.025), dAIC_hi=boot['dAIC'].quantile(.975),
        mice_zeta_mean=zetas_mice.mean(), mice_zeta_std=zetas_mice.std(),
    ))
    rank_fits[key] = dict(xs=xs, full=full, boot=boot, zetas_mice=zetas_mice)

    for i_m, x in enumerate(xs):
        plot_mouse_top_clones(x, g['color'], f"units/entropy_mouse_{key.replace(' ', '_')}_{i_m + 1}")

scaling_results = pd.DataFrame(rows)
group_fits      = pd.DataFrame(group_fits).set_index('group', drop=False)
scaling_results['barN_prediction'] = [r.N1 * np.sum(np.arange(1, r.L_act + 1) ** (-r.zeta))
                                      for r in scaling_results.itertuples()]


# ── Figures of one panel ──────────────────────────────────────────────────────
# EP: n ∝ exp(-alpha_hat dG), so dG = ln(N1/n) / alpha_hat in each unit (best clone at dG = 0),
# with alpha_hat = ALPHA (primary) or 2 ALPHA (recall). Pooled over the units of each group:
#   Omega(dG): clones per unit per unit dG (clone-weighted); predicted log-slope beta* = alpha_hat / zeta
#   q(dG):     pdf of cells, each unit weighted equally;      predicted log-slope beta* - alpha_hat
# Fits skip the dG = 0 bin (each unit's best clone, there by construction, as rank 1 in fit_both)
# and stop at the median doubleton energy ln(N1/2) / alpha_hat: beyond it the density is cut by sampling.
ALPHA_HAT     = {'primary': ALPHA, 'recall': 2 * ALPHA}
DG_BIN        = 0.2
MIN_MICE_RARE = 3   # rarefaction: drop n where fewer units than this have >= n cells

def fit_log_density(centers, dens, dG_max):
    sel = (dens > 0) & (centers > DG_BIN / 2) & (centers <= dG_max)
    slope, intercept = np.polyfit(centers[sel], np.log(dens[sel]), 1)
    return slope, intercept

def make_panel(panel, keys):
    os.makedirs(f'{output_plot}/{panel}', exist_ok=True)
    res   = scaling_results[scaling_results['group'].isin(keys)]
    fits  = group_fits.loc[keys]
    hue   = dict(zip(fits['group'], fits['color']))
    names = [f'{gf.Experiment}\n{gf.phenotype}' for gf in fits.itertuples()]

    def scatter(ax, x, y):
        sns.scatterplot(data=res, x=x, y=y, hue='group', style='Experiment', markers=EXP_MARKER, ax=ax,
                        s=150, palette=hue, edgecolors='black', alpha=0.8)

    # Rank curves: each unit, unit average, Zipf fit and rate-model fit; zeta violins
    fig_r,    ax_r    = plt.subplots(**fig_kw)
    fig_zeta, ax_zeta = plt.subplots(**fig_kw)
    for i, gf in enumerate(fits.itertuples()):
        rf = rank_fits[gf.group]
        for x in rf['xs']:
            ax_r.step(range(1, len(x) + 1), x / x[0], color=gf.color, alpha=.5, lw=0.5)
        x_avg = mean_rank_curve(rf['xs'])
        j = np.arange(1, len(x_avg) + 10)
        ax_r.plot(range(1, len(x_avg) + 1), x_avg, color=gf.color, markerfacecolor='None',
                  ms=12, ls='', marker=gf.marker, label=f'${gf.zeta_mean:.2f}$ ; {gf.Experiment} {gf.phenotype}')
        ax_r.plot(j, j ** (-gf.zeta_mean), color=gf.color, alpha=.8, lw=3, ls=EXP_LS[gf.Experiment])
        ax_r.plot(j, np.exp(rate_logx(np.log(j), rf['full']['z0'], rf['full']['c'])), color=gf.color, lw=2, ls=':')
        plot_zeta_violin(ax_zeta, rf['boot']['zeta'], rf['zetas_mice'], position=i, color=gf.color, marker=gf.marker)
    layout_ranking(ax_r)
    ax_r.legend(fontsize=20, loc=3)
    save(fig_r, f'{panel}/ranking_B_cells')
    layout_violins(ax_zeta, names, r'$\zeta$')
    save(fig_zeta, f'{panel}/zetas')

    # ── Scaling relations ──
    # N1 vs N
    fig, ax = plt.subplots(**fig_kw)
    scatter(ax, 'barN', 'N1')
    ax.set_xlim(left=0, right=1.05 * res['barN'].max())
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(title_fontsize=20, fontsize=20, loc=2)
    save(fig, f'{panel}/scaling_N1_Ntot')

    # L_act vs N, with the continuous-Zipf prediction (band: zeta +- 2 std)
    fig, ax = plt.subplots(**fig_kw)
    L = np.linspace(1, 400, 100)
    for gf in fits.itertuples():
        ax.plot(barN(L, gf.zeta_mean), L, color=gf.color, linestyle=EXP_LS[gf.Experiment])
        ax.fill_betweenx(L, barN(L, gf.zeta_mean - 2 * gf.zeta_std), barN(L, gf.zeta_mean + 2 * gf.zeta_std),
                         color=gf.color, alpha=0.1)
    scatter(ax, 'barN', 'L_act')
    ax.set_xlim(left=0, right=1.05 * res['barN'].max())
    ax.set_ylim(bottom=0, top=1.1 * res['L_act'].max())
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(title_fontsize=20, fontsize=20, loc=2)
    save(fig, f'{panel}/scaling_Lact_Ntot')

    # L_act vs N1, with L = N1^(1/zeta); shading: above / below L = N1
    fig, ax = plt.subplots(**fig_kw)
    N1 = np.linspace(0, 1.1 * res['N1'].max(), 100)
    ax.fill_between(N1, 0, N1, color=my_blue, alpha=0.2)
    ax.fill_between(N1, N1, 200, color=my_orange, alpha=0.2)
    for gf in fits.itertuples():
        ax.plot(N1, N1 ** (1 / gf.zeta_mean), color=gf.color, linestyle=EXP_LS[gf.Experiment])
    scatter(ax, 'N1', 'L_act')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xlim(left=0, right=1.1 * res['N1'].max())
    ax.set_ylim(bottom=0, top=1.1 * res['L_act'].max())
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(title_fontsize=20, fontsize=20, loc=0)
    save(fig, f'{panel}/scaling_Lact_N1')

    # N vs its prediction from N1, L_act and the unit's zeta
    fig, ax = plt.subplots(**fig_kw)
    nbar = np.linspace(res['barN'].min(), res['barN'].max(), 100)
    ax.plot(nbar, nbar, color='k', linestyle='--')
    scatter(ax, 'barN_prediction', 'barN')
    my_plot_layout(ax=ax, yscale='log', xscale='log', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(title_fontsize=20, fontsize=20, loc=4)
    save(fig, f'{panel}/Ntot_prediction')

    # ── Entropy ──
    x = np.linspace(1, 500, 1000)

    # S vs N
    fig, ax = plt.subplots(**fig_kw)
    scatter(ax, 'barN', 'S_corr')
    ax.plot(x, np.full_like(x, S_STAR), color=my_blue, linestyle='--')
    my_plot_layout(ax=ax, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax.set_xlabel(r'$N_{\mathrm{sample}}$', fontsize=30)
    ax.set_ylabel(r'$S$', fontsize=30)
    ax.set_xlim(left=2, right=1.05 * res['barN'].max())
    ax.tick_params(axis='both', labelsize=30)
    save(fig, f'{panel}/entropy_Ntot')

    # S vs L_act, with the maximum entropy S = ln(L_act) (all clones equal size)
    fig, ax = plt.subplots(**fig_kw)
    scatter(ax, 'L_act', 'S_corr')
    ax.plot(x, np.full_like(x, S_STAR), color=my_blue, linestyle='--')
    ax.plot(x, np.log(x), color='k', linestyle='--')
    my_plot_layout(ax=ax, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax.set_xlim(left=1, right=1.1 * res['L_act'].max())
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(title_fontsize=20, fontsize=20, loc=0)
    save(fig, f'{panel}/entropy_Lact')

    # Rarefaction: mean entropy of n subsampled cells per group (removes sequencing-depth differences).
    # Lines: Zipf prediction with the group's zeta and each unit's L_act, averaged over the same units.
    n_grid = np.unique(np.round(np.geomspace(2, res['barN'].max(), 30)).astype(int))
    fig, ax = plt.subplots(**fig_kw)
    for gf in fits.itertuples():
        grp    = res[res['group'] == gf.group]
        S_data = np.array([rarefied_entropy(c, n_grid) for c in grp['counts']])
        n_mice = np.sum(np.isfinite(S_data), axis=0)
        ok     = n_mice >= min(MIN_MICE_RARE, len(grp))
        mean   = np.nanmean(S_data[:, ok], axis=0)
        sem    = np.nanstd(S_data[:, ok], axis=0) / np.sqrt(n_mice[ok])
        ax.plot(n_grid[ok], mean, color=gf.color, ls='', marker=gf.marker, ms=8, label=f'{gf.Experiment} ; {gf.phenotype}')
        ax.fill_between(n_grid[ok], mean - sem, mean + sem, color=gf.color, alpha=0.2)
        S_zipf = np.array([zipf_rarefied_entropy(gf.zeta_hat, L, n_grid) for L in grp['L_act']])
        S_zipf[~np.isfinite(S_data)] = np.nan   # same units at each n as the data
        ax.plot(n_grid[ok], np.nanmean(S_zipf[:, ok], axis=0), color=gf.color, lw=2, ls=EXP_LS[gf.Experiment])
    my_plot_layout(ax=ax, yscale='linear', xscale='log', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax.set_xlabel(r'$n$ (cells subsampled)', fontsize=30)
    ax.set_ylabel(r'$\langle S(n) \rangle$', fontsize=30)
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(title_fontsize=20, fontsize=20, loc=2)
    save(fig, f'{panel}/entropy_rarefaction')

    # Renyi-2 entropy S2 = -ln(lambda), lambda = prob. two cells of the same unit share a clone.
    # Violin: population S2 = -ln(mean lambda over units), bootstrapped over units; dots: per-unit S2.
    # Black bar: Zipf prediction with the group's zeta and each unit's L_act.
    fig, ax = plt.subplots(**fig_kw)
    for i, gf in enumerate(fits.itertuples()):
        grp   = res[res['group'] == gf.group]
        lam   = grp['lambda2'].dropna().to_numpy()
        boots = [-np.log(np.mean(rng.choice(lam, len(lam)))) for _ in range(N_ENSEMBLE)]
        plot_zeta_violin(ax, boots, -np.log(lam[lam > 0]), position=i, color=gf.color, marker=gf.marker)
        lam_zipf = [np.sum(zipf_p(L, gf.zeta_hat) ** 2) for L in grp['L_act']]
        ax.scatter(i, -np.log(np.mean(lam_zipf)), color='k', marker='_', s=1500, lw=3)
    layout_violins(ax, names, r'$S_2$')
    save(fig, f'{panel}/entropy_S2')

    # ── Binding-energy densities (data version of manuscript Fig. potency A/B), one figure per group ──
    fig_clones, ax_clones = plt.subplots(**fig_kw)   # clone-weighted only: same slope beta* predicted for all
    for gf in fits.itertuples():
        grp     = res[res['group'] == gf.group]
        a_hat   = ALPHA_HAT[gf.Response]
        beta_pr = a_hat / gf.zeta_hat               # predicted beta*
        dG      = [np.log(r.N1 / r.counts) / a_hat for r in grp.itertuples()]
        w_cells = [r.counts / r.barN for r in grp.itertuples()]
        dG_fit  = np.median([np.log(r.N1 / 2) / a_hat for r in grp.itertuples()])
        bins    = np.arange(0, np.max(np.concatenate(dG)) + 2 * DG_BIN, DG_BIN) - DG_BIN / 2
        centers = (bins[:-1] + bins[1:]) / 2
        Omega   = np.histogram(np.concatenate(dG), bins)[0] / (len(grp) * DG_BIN)
        q       = np.histogram(np.concatenate(dG), bins, weights=np.concatenate(w_cells))[0] / (len(grp) * DG_BIN)

        # Entropy of q and the check S_c = <ln Omega>_q + S_dG[q] (manuscript Eq. S_c_S_E), on the binned densities
        P, K   = q * DG_BIN, Omega * DG_BIN         # cell fraction and clones per unit in each bin
        S_dG   = -np.sum(P[P > 0] * np.log(q[P > 0]))
        S_c_dG = np.sum(P[P > 0] * np.log(K[P > 0] / P[P > 0]))
        print(f'{panel} / {gf.group}: beta* predicted {beta_pr:.2f}; S_dG[q] = {S_dG:.2f}; '
              f'<ln Omega>_q + S_dG = {S_c_dG:.2f} vs mean clonal entropy S_c = {grp["S"].mean():.2f}')

        fig, ax = plt.subplots(**fig_kw)
        for dens, marker, mfc, slope_pr, name in [(Omega, 'o', 'None', beta_pr, r'\Omega'),
                                                  (q, 's', gf.color, beta_pr - a_hat, 'q')]:
            if np.sum((dens > 0) & (centers > DG_BIN / 2) & (centers <= dG_fit)) < 2:
                continue                                # too few bins to fit
            slope, intercept = fit_log_density(centers, dens, dG_fit)
            x_fit = np.linspace(0, dG_fit, 10)
            ax.plot(centers[dens > 0], dens[dens > 0], ls='', marker=marker, ms=10, color=gf.color, mfc=mfc,
                    label=fr'${name}$: slope ${slope:.2f}$ (pred. ${slope_pr:.2f}$)')
            ax.plot(x_fit, np.exp(intercept + slope * x_fit), color=gf.color, lw=2)
            x_mid = dG_fit / 2                       # prediction drawn through the fit at mid-range
            ax.plot(x_fit, np.exp(intercept + slope * x_mid + slope_pr * (x_fit - x_mid)), color='k', lw=1.5, ls='--')
            if name == r'\Omega':
                ax_clones.plot(centers[dens > 0], dens[dens > 0], ls='', marker=gf.marker, ms=10, color=gf.color,
                               label=fr'{gf.Experiment} {gf.phenotype}: slope ${slope:.2f}$')
                ax_clones.plot(x_fit, np.exp(intercept + slope * x_fit), color=gf.color, lw=2, ls=EXP_LS[gf.Experiment])
        ax.axvline(dG_fit, color='gray', ls=':')
        my_plot_layout(ax=ax, yscale='log', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
        ax.set_xlabel(r'$\Delta G$', fontsize=30)
        ax.set_ylabel(r'$\Omega(\Delta G),\ q(\Delta G)$', fontsize=30)
        ax.tick_params(axis='both', labelsize=30)
        ax.legend(fontsize=20, loc=0, title=fr'$S_{{\Delta G}}[q] = {S_dG:.2f}$', title_fontsize=20)
        save(fig, f"{panel}/distribution_dG_{gf.group.replace(' ', '_').replace('+', '-')}")

    my_plot_layout(ax=ax_clones, yscale='log', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax_clones.set_xlabel(r'$\Delta G$', fontsize=30)
    ax_clones.set_ylabel(r'$\Omega(\Delta G)$', fontsize=30)
    ax_clones.tick_params(axis='both', labelsize=30)
    ax_clones.legend(fontsize=20, loc=0)
    save(fig_clones, f'{panel}/distribution_dG_clones')
    plt.close('all')

for panel in RUN_PANELS:
    make_panel(panel, PANELS[panel])


# ── Save ──────────────────────────────────────────────────────────────────────
with open(root_dir + '/scaling_dict.pkl', 'wb') as f:
    pickle.dump(scaling_results.to_dict('list'), f)
group_fits.drop(columns=['color', 'marker']).to_csv(output_plot + '/zeta_stats.csv', index=False)
scaling_results.to_csv(root_dir + '/scaling_results.csv', index=False)
