"""
share_clones.py

For each pair of mice (within or across experiments), compute:
  - Clone sizes in mouse A vs. mouse B
  - Number of shared / unique-to-A / unique-to-B clones
  - Jaccard index

Clone identity is defined by CLONE_COLS (VDJ gene segments or CDR3 sequence).

Outputs:
  pairwise_overlap_{CLONE_DEF}.csv   — full stats for every pair
  jaccard_hist_{CLONE_DEF}.pdf       — Jaccard distribution by experiment pairing
  jaccard_heatmap_{CLONE_DEF}.pdf    — heatmap of Jaccard indices
  scatter_pairs_{CLONE_DEF}.pdf      — example clone-size scatter plots
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import seaborn as sns
import os
import sys
from itertools import combinations

sys.path.append('../../../library/')
from funcs import *

# ── Paths & constants ─────────────────────────────────────────────────────────
project    = 'exponential_proofreading'
subproject = 'data'
root_dir   = (
    f"/Users/robertomorantovar/Dropbox/Research_files/"
    f"Immune_system/{project}/{subproject}/mesin2020"
)
output_plot = (
    '/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/'
    f'Immune_System/_Repository/Figures/{project}/{subproject}/mesin2020/share_clones'
)
os.makedirs(output_plot, exist_ok=True)

EXCEL_FILE = root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx"

CLONE_DEF  = 'CDR3'                   # 'VDJ' or 'CDR3'
CLONE_COLS = ['V', 'J', 'D'] if CLONE_DEF == 'VDJ' else ['CDR3:']

N_SCATTER_PAIRS = 16                  # example pairs to scatter-plot


# ── Data loading ──────────────────────────────────────────────────────────────

def load_populations(sheet_name, header, clone_cols,
                     mouse_col='Mouse', phenotype_col=None,
                     filter_specs=None, prefix=''):
    """Return {pop_id: clone_count_Series} and {pop_id: metadata_dict}."""
    data = pd.read_excel(EXCEL_FILE, sheet_name=sheet_name, header=header)
    if filter_specs:
        for col, val in filter_specs.items():
            mask = data[col].isin(val) if isinstance(val, list) else data[col] == val
            data = data[mask]

    group_by = [mouse_col] + ([phenotype_col] if phenotype_col else [])
    counts_dict, meta_dict = {}, {}
    for key, grp in data.groupby(group_by):
        parts  = key if isinstance(key, tuple) else (key,)
        pop_id = prefix + '_'.join(str(p) for p in parts)
        counts_dict[pop_id] = grp.groupby(clone_cols).size()
        meta_dict[pop_id]   = {
            'mouse':     str(parts[0]),
            'phenotype': str(parts[1]) if len(parts) > 1 else 'all',
            'experiment': prefix.rstrip('_'),
        }
    return counts_dict, meta_dict


# ── Clone overlap ─────────────────────────────────────────────────────────────

def compute_overlap(counts_a, counts_b):
    """Return DataFrame with clone sizes in population A and B (0 if absent)."""
    idx = counts_a.index.union(counts_b.index)
    return pd.DataFrame({
        'size_a': counts_a.reindex(idx, fill_value=0).values,
        'size_b': counts_b.reindex(idx, fill_value=0).values,
    })


def overlap_stats(df):
    shared   = ((df['size_a'] > 0) & (df['size_b'] > 0)).sum()
    unique_a = ((df['size_a'] > 0) & (df['size_b'] == 0)).sum()
    unique_b = ((df['size_a'] == 0) & (df['size_b'] > 0)).sum()
    total    = shared + unique_a + unique_b
    jaccard  = shared / total if total > 0 else 0.0
    # fraction of cells (not clones) in shared clones
    cells_a  = df['size_a'].sum()
    cells_b  = df['size_b'].sum()
    shared_cells_a = df.loc[(df['size_a'] > 0) & (df['size_b'] > 0), 'size_a'].sum()
    shared_cells_b = df.loc[(df['size_a'] > 0) & (df['size_b'] > 0), 'size_b'].sum()
    frac_cells_a   = shared_cells_a / cells_a if cells_a > 0 else 0.0
    frac_cells_b   = shared_cells_b / cells_b if cells_b > 0 else 0.0
    return dict(
        shared=int(shared), unique_a=int(unique_a), unique_b=int(unique_b),
        jaccard=jaccard,
        frac_shared_cells_a=frac_cells_a, frac_shared_cells_b=frac_cells_b,
    )


# ── Plot helpers ──────────────────────────────────────────────────────────────

def scatter_pair(ax, df, label_a, label_b):
    """Clone-size scatter using log(1+n) so zeros land on the axes."""
    xa = np.log1p(df['size_a'].values)
    xb = np.log1p(df['size_b'].values)
    shared   = (df['size_a'] > 0) & (df['size_b'] > 0)
    unique_a = (df['size_a'] > 0) & (df['size_b'] == 0)
    unique_b = (df['size_a'] == 0) & (df['size_b'] > 0)

    ax.scatter(xa[shared],   xb[shared],   color='steelblue', s=18, alpha=0.6,
               label=f"shared ({shared.sum()})", zorder=3)
    ax.scatter(xa[unique_a], xb[unique_a], color='tomato',    s=12, alpha=0.5,
               marker='|', label=f"only A ({unique_a.sum()})")
    ax.scatter(xa[unique_b], xb[unique_b], color='seagreen',  s=12, alpha=0.5,
               marker='_', label=f"only B ({unique_b.sum()})")

    lim = max(xa.max(), xb.max()) * 1.05 + 0.1
    ax.plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.4)
    ax.set_xlim(-0.1, lim); ax.set_ylim(-0.1, lim)

    ticks_raw = [0, 1, 2, 5, 10, 20, 50, 100]
    ticks_tr  = [np.log1p(t) for t in ticks_raw if np.log1p(t) <= lim]
    ticks_raw = ticks_raw[:len(ticks_tr)]
    ax.set_xticks(ticks_tr, [str(t) for t in ticks_raw], fontsize=8)
    ax.set_yticks(ticks_tr, [str(t) for t in ticks_raw], fontsize=8)
    ax.set_xlabel(label_a, fontsize=9)
    ax.set_ylabel(label_b, fontsize=9)
    ax.legend(fontsize=7, loc='upper left')


# ── Load all experiments ──────────────────────────────────────────────────────

all_counts, all_meta = {}, {}

def _merge(counts_dict, meta_dict):
    all_counts.update(counts_dict)
    all_meta.update(meta_dict)

c, m = load_populations('Photoactivation CGG', header=1, clone_cols=CLONE_COLS,
                        filter_specs={'Figure': 1}, prefix='exp1_')
_merge(c, m)

c, m = load_populations('Fate-mapping CGG', header=1, clone_cols=CLONE_COLS,
                        filter_specs={'Figure': '4A'}, prefix='exp2_')
_merge(c, m)

c, m = load_populations('Fate-mapping CGG', header=1, clone_cols=CLONE_COLS,
                        phenotype_col='Phenotype',
                        filter_specs={'Figure': '4C-H'}, prefix='exp3_')
_merge(c, m)

c, m = load_populations('Influenza_IGH', header=2, clone_cols=CLONE_COLS,
                        mouse_col='Experiment / Mouse', phenotype_col='Sort',
                        filter_specs={'Sort': ['GC', 'GC + fm', 'PB + fm', 'M']},
                        prefix='exp4_')
_merge(c, m)

pop_ids = list(all_counts.keys())
print(f"Total populations: {len(pop_ids)}")


# ── Pairwise overlaps ─────────────────────────────────────────────────────────

records = []
for id_a, id_b in combinations(pop_ids, 2):
    df    = compute_overlap(all_counts[id_a], all_counts[id_b])
    stats = overlap_stats(df)
    ma, mb = all_meta[id_a], all_meta[id_b]
    records.append({
        'pop_a': id_a, 'pop_b': id_b,
        'exp_a': ma['experiment'], 'exp_b': mb['experiment'],
        'phenotype_a': ma['phenotype'], 'phenotype_b': mb['phenotype'],
        'same_exp': ma['experiment'] == mb['experiment'],
        'same_phenotype': ma['phenotype'] == mb['phenotype'],
        **stats,
    })

results = pd.DataFrame(records)
results.to_csv(output_plot + f'/pairwise_overlap_{CLONE_DEF}.csv', index=False)
print(results[['pop_a', 'pop_b', 'shared', 'jaccard']].head(10).to_string())


# ── Jaccard histogram ─────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5),
                         gridspec_kw={'left': .10, 'right': .95, 'bottom': .15, 'top': .92})

for ax, (col, title) in zip(axes, [
    ('same_exp',       'Same vs. different experiment'),
    ('same_phenotype', 'Same vs. different phenotype'),
]):
    bins = np.linspace(0, results['jaccard'].max() + 0.02, 20)
    ax.hist(results.loc[ results[col], 'jaccard'], bins=bins, alpha=0.65,
            color='steelblue', label='Same')
    ax.hist(results.loc[~results[col], 'jaccard'], bins=bins, alpha=0.65,
            color='tomato',    label='Different')
    ax.set_xlabel(r'Jaccard index', fontsize=13)
    ax.set_ylabel(r'Number of pairs', fontsize=13)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=11)

fig.suptitle(f'Clone overlap — {CLONE_DEF} definition', fontsize=13)
fig.savefig(output_plot + f'/jaccard_hist_{CLONE_DEF}.pdf', bbox_inches='tight')


# ── Jaccard heatmap ───────────────────────────────────────────────────────────

jaccard_mat = pd.DataFrame(np.eye(len(pop_ids)), index=pop_ids, columns=pop_ids)
for _, row in results.iterrows():
    jaccard_mat.loc[row['pop_a'], row['pop_b']] = row['jaccard']
    jaccard_mat.loc[row['pop_b'], row['pop_a']] = row['jaccard']

fig_hm, ax_hm = plt.subplots(figsize=(max(8, len(pop_ids) * 0.35 + 2),
                                       max(6, len(pop_ids) * 0.35 + 1)),
                              gridspec_kw={'left': .22, 'right': .97, 'bottom': .22, 'top': .96})
sns.heatmap(jaccard_mat, ax=ax_hm, cmap='Blues', vmin=0, vmax=1,
            xticklabels=True, yticklabels=True,
            cbar_kws={'label': 'Jaccard index', 'shrink': 0.7})
ax_hm.tick_params(axis='both', labelsize=7)
ax_hm.set_title(f'Pairwise Jaccard — {CLONE_DEF}', fontsize=12)
fig_hm.savefig(output_plot + f'/jaccard_heatmap_{CLONE_DEF}.pdf', bbox_inches='tight')


# ── Example scatter plots ─────────────────────────────────────────────────────

n_cols = 4
n_rows = (N_SCATTER_PAIRS + n_cols - 1) // n_cols
fig_sc, axes_sc = plt.subplots(n_rows, n_cols,
                                figsize=(5 * n_cols, 4.5 * n_rows),
                                gridspec_kw={'left': .06, 'right': .97,
                                             'bottom': .08, 'top': .94,
                                             'hspace': 0.45, 'wspace': 0.35})
axes_sc = axes_sc.flatten()

sample = results.nlargest(N_SCATTER_PAIRS, 'shared')

for i, (_, row) in enumerate(sample.iterrows()):
    df = compute_overlap(all_counts[row['pop_a']], all_counts[row['pop_b']])
    la = row['pop_a'].replace('_', r'\_')
    lb = row['pop_b'].replace('_', r'\_')
    scatter_pair(axes_sc[i], df, la, lb)
    axes_sc[i].set_title(
        rf"J$={row['jaccard']:.2f}$ | sh$={row['shared']}$", fontsize=9)

for ax in axes_sc[len(sample):]:
    ax.set_visible(False)

fig_sc.suptitle(f'Clone-size scatter ({CLONE_DEF}) — log$(1+n)$ scale', fontsize=12)
fig_sc.savefig(output_plot + f'/scatter_pairs_{CLONE_DEF}.pdf', bbox_inches='tight')


