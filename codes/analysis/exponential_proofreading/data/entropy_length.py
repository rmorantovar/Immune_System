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
import seaborn as sns
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


# ── CDR3 track ───────────────────────────────────────────────────────────────
#
# CDR3 sequences have variable length, so "using k characters" requires
# choosing an anchor end.  We run both:
#   'left'  — first k characters, anchored at the V/Cys end
#   'right' — last  k characters, anchored at the J/Trp end
# Both are approximations; true positional equivalence would need IMGT
# numbering (ANARCI) or MSA (MUSCLE/MAFFT).  The two curves together still
# reveal which end of the CDR3 carries more clone-discriminating information.

def entropy_vs_cdr3_length(data, mouse_col='Mouse', cdr3_col=CDR3_COL, mode='left'):
    """
    Per-mouse Shannon entropy as a function of CDR3 prefix/suffix length k.

    mode: 'left'  — str[:k]  (anchored at V / N-terminal end)
          'right' — str[-k:] (anchored at J / C-terminal end)

    Returns (max_len, {k: {mouse_id: entropy}}).
    """
    data = data.dropna(subset=[cdr3_col]).copy()
    data[cdr3_col] = data[cdr3_col].astype(str).str.strip()
    max_len = int(data[cdr3_col].str.len().max())

    result = {}
    for k in range(1, max_len + 1):
        data['_trunc'] = (data[cdr3_col].str[:k]
                          if mode == 'left'
                          else data[cdr3_col].str[-k:])
        result[k] = per_mouse_entropy(data, ['_trunc'], mouse_col=mouse_col)
    data.drop(columns=['_trunc'], inplace=True)
    return max_len, result


# CDR3 length distribution (informational)
def cdr3_length_hist(data, cdr3_col=CDR3_COL, ax=None, color='gray', label=''):
    lengths = data[cdr3_col].dropna().astype(str).str.strip().str.len()
    if ax is None:
        _, ax = plt.subplots()
    ax.hist(lengths, bins=range(lengths.min(), lengths.max() + 2),
            color=color, alpha=0.5, label=label, density=True)
    return lengths.mode()[0]          # return modal length


# ── CDR3 alignment helpers ───────────────────────────────────────────────────

def center_align_sequences(seqs):
    """
    Insert gaps in the middle of shorter sequences so all reach max_len.
    CDR3 is bounded by conserved Cys (N-term) and Trp/Phe (C-term), making
    the middle the most variable region — padding there preserves both anchors.
    """
    max_len = max(len(s) for s in seqs)
    aligned = []
    for s in seqs:
        n_gaps = max_len - len(s)
        pivot  = (len(s) + 1) // 2          # keep longer half on the left
        aligned.append(s[:pivot] + '-' * n_gaps + s[pivot:])
    return aligned


def _write_fasta(path, id_seq_pairs):
    with open(path, 'w') as f:
        for seq_id, seq in id_seq_pairs:
            f.write(f'>{seq_id}\n{seq}\n')


def _read_fasta(path):
    records, cur_id, cur_seq = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id is not None:
                    records[cur_id] = ''.join(cur_seq)
                cur_id, cur_seq = line[1:].split()[0], []
            else:
                cur_seq.append(line)
    if cur_id is not None:
        records[cur_id] = ''.join(cur_seq)
    return records


def muscle_align_sequences(seqs, muscle_cmd='muscle'):
    """
    Align via MUSCLE external binary. Returns {original_seq: aligned_seq}.
    Tries MUSCLE v5 syntax first (-align / -output), then v3 (-in / -out).
    Raises FileNotFoundError if MUSCLE is not on PATH.
    Install MUSCLE: https://github.com/rcedgar/muscle/releases
    """
    import tempfile
    import subprocess

    unique = list(dict.fromkeys(seqs))      # deduplicate, preserve order
    with tempfile.NamedTemporaryFile(suffix='.fa', mode='w', delete=False) as fin:
        in_path = fin.name
    _write_fasta(in_path, [(f's{i}', s) for i, s in enumerate(unique)])
    out_path = in_path.replace('.fa', '_aln.fa')

    r5 = subprocess.run([muscle_cmd, '-align', in_path, '-output', out_path],
                        capture_output=True)
    if r5.returncode != 0:
        stderr = r5.stderr.decode(errors='replace').strip()
        raise RuntimeError(f"MUSCLE failed (exit {r5.returncode}). stderr:\n{stderr}")

    aln_map = _read_fasta(out_path)
    return {s: aln_map[f's{i}'] for i, s in enumerate(unique)}


_VALID_AA = frozenset('ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx')

def apply_alignment(data, cdr3_col, method='center'):
    """
    Return a copy of data with the CDR3 column replaced by gap-padded aligned
    strings of uniform length.

    method: 'center' — center-alignment (no external dependency)
            'muscle' — MUSCLE MSA   (requires MUSCLE binary on PATH)
    """
    data = data.dropna(subset=[cdr3_col]).copy()
    seqs = data[cdr3_col].astype(str).str.strip()

    # Require non-empty sequences made of standard amino acid characters only
    valid = seqs.apply(lambda s: len(s) > 0 and all(c in _VALID_AA for c in s))
    n_dropped = (~valid).sum()
    if n_dropped > 0:
        bad_chars = set(
            c for s in seqs[~valid] for c in s if c not in _VALID_AA
        )
        print(f"  [{method}] dropping {n_dropped} rows with non-AA chars: "
              f"{bad_chars}")
    data = data[valid]
    seqs = seqs[valid]
    uniq = seqs.unique().tolist()

    if method == 'center':
        mapping = dict(zip(uniq, center_align_sequences(uniq)))
    elif method == 'muscle':
        mapping = muscle_align_sequences(uniq)
    else:
        raise ValueError(f"Unknown alignment method: {method!r}")

    data[cdr3_col] = seqs.map(mapping)
    return data


# ── Run CDR3 track for all experiments ───────────────────────────────────────

fig_cdr3, axes_cdr3 = plt.subplots(
    2, len(EXP_SPECS),
    figsize=(4.5 * len(EXP_SPECS), 9),
    gridspec_kw={'left': .08, 'right': .97, 'bottom': .08,
                 'top': .93, 'hspace': 0.40, 'wspace': 0.30})

fig_len, axes_len = plt.subplots(
    1, len(EXP_SPECS),
    figsize=(4 * len(EXP_SPECS), 3.5),
    gridspec_kw={'left': .08, 'right': .97, 'bottom': .18, 'top': .88, 'wspace': 0.35})

for col, exp in enumerate(EXP_SPECS):
    data_exp = load_raw(exp['sheet'], exp['header'], exp.get('filter_specs'))
    mouse_col = exp['mouse_col']
    mice = data_exp[mouse_col].unique()

    # CDR3 length distribution
    modal_len = cdr3_length_hist(data_exp, ax=axes_len[col],
                                 color=exp['color'], label=exp['label'])
    axes_len[col].axvline(modal_len, color=exp['color'], ls='--', lw=1.5)
    axes_len[col].set_xlabel('CDR3 length', fontsize=11)
    axes_len[col].set_ylabel('Density', fontsize=11)
    axes_len[col].set_title(exp['label'], fontsize=10)

    for row, mode in enumerate(['left', 'right']):
        ax = axes_cdr3[row, col]
        max_len, S_k = entropy_vs_cdr3_length(data_exp, mouse_col=mouse_col, mode=mode)
        k_vals = np.arange(1, max_len + 1)

        # per-mouse curves (thin)
        for mouse in mice:
            s_mouse = [S_k[k].get(str(mouse), np.nan) for k in k_vals]
            ax.plot(k_vals, s_mouse, color=exp['color'], alpha=0.25, lw=0.8)

        # mean ± std (thick)
        means = np.array([np.nanmean(list(S_k[k].values())) for k in k_vals])
        stds  = np.array([np.nanstd( list(S_k[k].values())) for k in k_vals])
        ax.plot(k_vals, means, color=exp['color'], lw=2.2)
        ax.fill_between(k_vals, means - stds, means + stds,
                        color=exp['color'], alpha=0.2)

        # mark modal CDR3 length
        ax.axvline(modal_len, color='k', ls=':', lw=1, alpha=0.5)

        anchor = r'V/Cys end' if mode == 'left' else r'J/Trp end'
        ax.set_xlabel(rf'CDR3 prefix length $k$ ({anchor})', fontsize=10)
        ax.set_ylabel(r'$S$ (Shannon entropy)', fontsize=10)
        if row == 0:
            ax.set_title(exp['label'], fontsize=10)

fig_cdr3.suptitle(r'Entropy vs.\ CDR3 length (left: V-anchored, right: J-anchored)',
                  fontsize=12)
fig_cdr3.savefig(output_plot + '/entropy_CDR3_length.pdf', bbox_inches='tight')

fig_len.suptitle('CDR3 length distribution', fontsize=12)
fig_len.savefig(output_plot + '/CDR3_length_dist.pdf', bbox_inches='tight')


# ── CDR3 track: alignment comparison ─────────────────────────────────────────
# Each panel overlays four curves for one experiment:
#   no-alignment left  (dashed)
#   no-alignment right (dotted)
#   center alignment   (solid)
#   MUSCLE alignment   (dash-dot, skipped gracefully if MUSCLE not on PATH)

ALIGN_METHODS = [
    dict(label=r'No align, left',  align=None,     mode='left',  ls='--', lw=1.8),
    dict(label=r'No align, right', align=None,     mode='right', ls=':',  lw=1.8),
    dict(label=r'Center align',    align='center', mode='left',  ls='-',  lw=2.2),
    dict(label=r'MUSCLE',          align='muscle', mode='left',  ls='-.', lw=2.2),
]

fig_msa, axes_msa = plt.subplots(
    1, len(EXP_SPECS),
    figsize=(4.5 * len(EXP_SPECS), 4.5),
    gridspec_kw={'left': .08, 'right': .97, 'bottom': .14,
                 'top': .88, 'wspace': 0.30})

for col, exp in enumerate(EXP_SPECS):
    ax  = axes_msa[col]
    data_exp  = load_raw(exp['sheet'], exp['header'], exp.get('filter_specs'))
    mouse_col = exp['mouse_col']
    color     = exp['color']

    modal_len = int(
        data_exp[CDR3_COL].dropna().astype(str).str.strip().str.len().mode()[0])
    ax.axvline(modal_len, color='k', ls=':', lw=1, alpha=0.4, label='modal len')

    for meth in ALIGN_METHODS:
        if meth['align'] is not None:
            try:
                data_in = apply_alignment(data_exp, CDR3_COL, method=meth['align'])
            except FileNotFoundError:
                print(f"  [{meth['align']}] skipped — binary not found on PATH")
                continue
            except Exception as e:
                print(f"  [{meth['align']}] skipped — {e}")
                continue
        else:
            data_in = data_exp

        max_len_m, S_k = entropy_vs_cdr3_length(
            data_in, mouse_col=mouse_col, mode=meth['mode'])
        k_vals = np.arange(1, max_len_m + 1)
        means  = np.array([np.nanmean(list(S_k[k].values())) for k in k_vals])
        stds   = np.array([np.nanstd( list(S_k[k].values())) for k in k_vals])

        ax.plot(k_vals, means, color=color,
                ls=meth['ls'], lw=meth['lw'], label=meth['label'])
        ax.fill_between(k_vals, means - stds, means + stds,
                        color=color, alpha=0.12)

    ax.set_xlabel(r'CDR3 length $k$', fontsize=11)
    ax.set_ylabel(r'$S$ (Shannon entropy)', fontsize=11)
    ax.set_title(exp['label'], fontsize=10)
    ax.legend(fontsize=8, loc='lower right')
    ax.tick_params(labelsize=10)

fig_msa.suptitle(r'Entropy vs.\ CDR3 length — alignment comparison', fontsize=12)
fig_msa.savefig(output_plot + '/entropy_CDR3_alignment_comparison.pdf',
                bbox_inches='tight')


