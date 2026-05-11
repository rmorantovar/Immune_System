"""
entropy_length.py

How does the Shannon entropy of a mouse's B-cell population depend on
the sequence definition used to identify clones?

VDJ track  — subsets of {V, D, J}: V, D, J, VD, VJ, DJ, VDJ.
             Entropy is computed per mouse, then averaged across mice.

CDR3 track — STUB.
             Variable-length CDR3 sequences require alignment before
             prefix-truncation is biologically meaningful (residue at
             position k should correspond to the same structural position
             across sequences). Recommended options:
               1. IMGT numbering via ANARCI (pip install anarci)
               2. Multiple sequence alignment: BioPython PairwiseAligner
                  or an external MUSCLE/MAFFT binary.
             The function stub below shows the intended interface; fill in
             the alignment step and uncomment to run.
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib.gridspec as gridspec
import pandas as pd
import os
import sys
from itertools import chain, combinations

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
    f'Immune_System/_Repository/Figures/{project}/{subproject}/mesin2020/entropy_length'
)
os.makedirs(output_plot, exist_ok=True)

EXCEL_FILE = root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx"
CDR3_COL   = 'CDR3:'


# ── Helpers ───────────────────────────────────────────────────────────────────

def shannon_entropy(counts):
    counts = np.asarray(counts, dtype=float)
    N = counts.sum()
    if N == 0:
        return np.nan
    p = counts / N
    return float(-np.sum(p[p > 0] * np.log(p[p > 0])))


def per_mouse_entropy(data, clone_cols, mouse_col='Mouse'):
    """Return {mouse_id: entropy} for each mouse in data."""
    result = {}
    for mouse, grp in data.groupby(mouse_col):
        counts = grp.groupby(clone_cols).size().values
        result[str(mouse)] = shannon_entropy(counts)
    return result


def load_raw(sheet_name, header, filter_specs=None):
    data = pd.read_excel(EXCEL_FILE, sheet_name=sheet_name, header=header)
    if filter_specs:
        for col, val in filter_specs.items():
            mask = data[col].isin(val) if isinstance(val, list) else data[col] == val
            data = data[mask]
    return data


def powerset_nonempty(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))


# ── VDJ track ─────────────────────────────────────────────────────────────────

VDJ_GENES   = ['V', 'D', 'J']
VDJ_SUBSETS = list(powerset_nonempty(VDJ_GENES))   # 7 subsets
SUBSET_LABELS = [''.join(s) for s in VDJ_SUBSETS]  # V, D, J, VD, VJ, DJ, VDJ

# One color per experiment
EXP_SPECS = [
    dict(label='Exp 1 (naive GC)',   color=my_red,
         sheet='Photoactivation CGG', header=1, mouse_col='Mouse',
         filter_specs={'Figure': 1}),
    dict(label='Exp 2 (recall GC+fm)', color=my_blue,
         sheet='Fate-mapping CGG', header=1, mouse_col='Mouse',
         filter_specs={'Figure': '4A'}),
    dict(label='Exp 3 (recall)',     color=my_cyan,
         sheet='Fate-mapping CGG', header=1, mouse_col='Mouse',
         filter_specs={'Figure': '4C-H'}),
    dict(label='Exp 4 (influenza GC)', color=my_purple,
         sheet='Influenza_IGH', header=2, mouse_col='Experiment / Mouse',
         filter_specs={'Sort': 'GC'}),
]

fig_vdj, ax_vdj = plt.subplots(
    figsize=(10, 5.5),
    gridspec_kw={'left': .10, 'right': .97, 'bottom': .12, 'top': .93})

x_pos = np.arange(len(SUBSET_LABELS))

for exp in EXP_SPECS:
    data = load_raw(exp['sheet'], exp['header'], exp.get('filter_specs'))
    S_per_subset = []
    for subset in VDJ_SUBSETS:
        em = per_mouse_entropy(data, list(subset), mouse_col=exp['mouse_col'])
        S_per_subset.append(list(em.values()))

    means = [np.nanmean(s) for s in S_per_subset]
    stds  = [np.nanstd(s)  for s in S_per_subset]

    ax_vdj.errorbar(x_pos, means, yerr=stds,
                    fmt='o-', color=exp['color'], capsize=4,
                    ms=8, lw=1.8, label=exp['label'])

    for xi, vals in zip(x_pos, S_per_subset):
        ax_vdj.scatter(
            np.random.normal(xi, 0.07, len(vals)), vals,
            color=exp['color'], alpha=0.35, s=25, zorder=4)

ax_vdj.set_xticks(x_pos, SUBSET_LABELS, fontsize=13)
ax_vdj.set_xlabel(r'Gene combination', fontsize=14)
ax_vdj.set_ylabel(r'$S$ (Shannon entropy)', fontsize=14)
ax_vdj.set_title(r'Entropy vs.\ VDJ gene combination', fontsize=13)
ax_vdj.legend(fontsize=10, loc='upper left')
ax_vdj.tick_params(axis='y', labelsize=12)
fig_vdj.savefig(output_plot + '/entropy_VDJ_subsets.pdf', bbox_inches='tight')


# ── VDJ track: per-mouse heatmap ─────────────────────────────────────────────
# Visualise how entropy changes across subsets for individual mice.

import seaborn as sns

rows = []
for exp in EXP_SPECS:
    data = load_raw(exp['sheet'], exp['header'], exp.get('filter_specs'))
    for subset, label in zip(VDJ_SUBSETS, SUBSET_LABELS):
        em = per_mouse_entropy(data, list(subset), mouse_col=exp['mouse_col'])
        for mouse, s in em.items():
            rows.append({'experiment': exp['label'], 'mouse': mouse,
                         'subset': label, 'S': s})

df_vdj = pd.DataFrame(rows)
df_vdj.to_csv(output_plot + '/entropy_VDJ_subsets.csv', index=False)

pivot = df_vdj.pivot_table(index=['experiment', 'mouse'], columns='subset', values='S')
pivot = pivot[SUBSET_LABELS]   # enforce ordered columns

fig_heat, ax_heat = plt.subplots(
    figsize=(9, max(4, len(pivot) * 0.35 + 1.5)),
    gridspec_kw={'left': .32, 'right': .97, 'bottom': .12, 'top': .94})
sns.heatmap(pivot, ax=ax_heat, cmap='YlOrRd', annot=True, fmt='.2f',
            annot_kws={'size': 7},
            cbar_kws={'label': r'$S$', 'shrink': 0.7})
ax_heat.set_xlabel('Gene combination', fontsize=12)
ax_heat.set_ylabel('')
ax_heat.tick_params(axis='x', labelsize=11)
ax_heat.tick_params(axis='y', labelsize=7)
ax_heat.set_title(r'Per-mouse $S$ across VDJ subsets', fontsize=11)
fig_heat.savefig(output_plot + '/entropy_VDJ_heatmap.pdf', bbox_inches='tight')


# ── CDR3 track (STUB) ────────────────────────────────────────────────────────

def entropy_vs_cdr3_length(data, mouse_col='Mouse', cdr3_col=CDR3_COL):
    """
    STUB — requires sequence alignment before use.

    CDR3 sequences differ in length. Truncating to the first k characters
    conflates structurally non-equivalent positions unless sequences are
    first aligned (anchored at the conserved Cys/Trp ends of CDR3).

    To implement:
      1. Align sequences with ANARCI (IMGT numbering) or
         BioPython PairwiseAligner / MUSCLE.
      2. Replace raw CDR3 strings in `data` with their aligned forms
         (gap-padded to a common length).
      3. Uncomment the loop below and run.

    Intended interface:
      Returns a dict {k: {mouse: entropy}} for k in 1..max_aligned_length.
    """
    # data = data.dropna(subset=[cdr3_col]).copy()
    # data[cdr3_col] = data[cdr3_col].astype(str)
    # # ← insert alignment step here; aligned sequences should all have the
    # #   same length and use '-' for gaps.
    # max_len = data[cdr3_col].str.len().max()
    # result  = {}
    # for k in range(1, max_len + 1):
    #     data['_cdr3_k'] = data[cdr3_col].str[:k]
    #     result[k] = per_mouse_entropy(data, ['_cdr3_k'], mouse_col=mouse_col)
    # data.drop(columns=['_cdr3_k'], inplace=True)
    # return result
    raise NotImplementedError("Implement sequence alignment before calling this function.")


print("\nCDR3 entropy track is a stub — see entropy_vs_cdr3_length() docstring.")

