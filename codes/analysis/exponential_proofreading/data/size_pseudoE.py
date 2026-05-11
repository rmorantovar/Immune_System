"""
size_pseudoE.py

1. Align CDR3 sequences via MUSCLE (center-align as fallback).
2. Derive a consensus sequence (column-wise mode of the MSA).
3. Compute Hamming distance of each clone's CDR3 to the consensus.
4. Within each mouse: Spearman r(log clone_size, distance) for alignment
   prefix lengths k = 1 … max_k.

The distance to consensus acts as a pseudo-energy: if affinity maturation
drives selected (large) clones toward a common motif, we expect r < 0.

Outputs:
  corr_vs_length_{ALIGN}.pdf   — mean Spearman r vs. k (mean ± SD across mice)
  size_vs_distance_{ALIGN}.pdf — scatter log(clone_size) vs. distance at full k
  corr_results_{ALIGN}.csv     — per-mouse per-k Spearman r and p-value
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
from scipy import stats
import os, sys, tempfile, subprocess

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
    f'Immune_System/_Repository/Figures/{project}/{subproject}/mesin2020/size_pseudoE'
)
os.makedirs(output_plot, exist_ok=True)

EXCEL_FILE = root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx"
CDR3_COL   = 'CDR3:'
ALIGN      = 'muscle'   # 'muscle' or 'center'

EXP_SPECS = [
    dict(label='Exp 1 (naive GC)',     color=my_red,
         sheet='Photoactivation CGG',  header=1, mouse_col='Mouse',
         filter_specs={'Figure': 1}),
    dict(label='Exp 2 (recall GC+fm)', color=my_blue,
         sheet='Fate-mapping CGG',     header=1, mouse_col='Mouse',
         filter_specs={'Figure': '4A'}),
    dict(label='Exp 3 (recall)',       color=my_cyan,
         sheet='Fate-mapping CGG',     header=1, mouse_col='Mouse',
         filter_specs={'Figure': '4C-H'}),
    dict(label='Exp 4 (influenza GC)', color=my_purple,
         sheet='Influenza_IGH',        header=2, mouse_col='Experiment / Mouse',
         filter_specs={'Sort': 'GC'}),
]

# ── Sequence helpers ──────────────────────────────────────────────────────────

_VALID_AA = frozenset('ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx')


def load_raw(sheet_name, header, filter_specs=None):
    data = pd.read_excel(EXCEL_FILE, sheet_name=sheet_name, header=header)
    if filter_specs:
        for col, val in filter_specs.items():
            mask = data[col].isin(val) if isinstance(val, list) else data[col] == val
            data = data[mask]
    return data


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


def muscle_align(seqs, muscle_cmd='muscle'):
    """Returns {raw_seq: aligned_seq} via MUSCLE v5."""
    unique = list(dict.fromkeys(seqs))
    with tempfile.NamedTemporaryFile(suffix='.fa', delete=False) as fin:
        in_path = fin.name
    out_path = in_path.replace('.fa', '_aln.fa')
    _write_fasta(in_path, [(f's{i}', s) for i, s in enumerate(unique)])
    r = subprocess.run([muscle_cmd, '-align', in_path, '-output', out_path],
                       capture_output=True)
    if r.returncode != 0:
        stderr = r.stderr.decode(errors='replace').strip()
        raise RuntimeError(f"MUSCLE failed (exit {r.returncode}):\n{stderr}")
    aln = _read_fasta(out_path)
    return {s: aln[f's{i}'] for i, s in enumerate(unique)}


def center_align(seqs):
    """Returns {raw_seq: aligned_seq} via gap insertion at the midpoint."""
    unique = list(dict.fromkeys(seqs))
    max_len = max(len(s) for s in unique)
    result  = {}
    for s in unique:
        n_gaps = max_len - len(s)
        pivot  = (len(s) + 1) // 2
        result[s] = s[:pivot] + '-' * n_gaps + s[pivot:]
    return result


def align_sequences(seqs, method='muscle'):
    if method == 'muscle':
        return muscle_align(seqs)
    return center_align(seqs)


def consensus_sequence(aligned_seqs):
    """Column-wise mode, ignoring gap characters."""
    if not aligned_seqs:
        return ''
    L = max(len(s) for s in aligned_seqs)
    cons = []
    for i in range(L):
        col = [s[i] for s in aligned_seqs if i < len(s) and s[i] != '-']
        cons.append(max(set(col), key=col.count) if col else '-')
    return ''.join(cons)


# ── Figures ───────────────────────────────────────────────────────────────────

N = len(EXP_SPECS)

fig_corr, axes_corr = plt.subplots(
    1, N, figsize=(4.5 * N, 4.5),
    gridspec_kw={'left': .09, 'right': .97, 'bottom': .15,
                 'top': .88, 'wspace': 0.35})

fig_sc, axes_sc = plt.subplots(
    1, N, figsize=(4.2 * N, 4.2),
    gridspec_kw={'left': .09, 'right': .97, 'bottom': .15,
                 'top': .88, 'wspace': 0.35})

all_records = []

# ── Main loop ─────────────────────────────────────────────────────────────────

for col, exp in enumerate(EXP_SPECS):
    ax_corr   = axes_corr[col]
    ax_sc     = axes_sc[col]
    color     = exp['color']
    mouse_col = exp['mouse_col']

    # ── Load & clean ──────────────────────────────────────────────────────────
    data = load_raw(exp['sheet'], exp['header'], exp.get('filter_specs'))
    data = data.dropna(subset=[CDR3_COL]).copy()
    data[CDR3_COL] = data[CDR3_COL].astype(str).str.strip()
    valid = data[CDR3_COL].apply(
        lambda s: len(s) > 0 and all(c in _VALID_AA for c in s))
    n_dropped = (~valid).sum()
    if n_dropped:
        print(f"[{exp['label']}] dropping {n_dropped} rows with non-AA chars")
    data = data[valid]

    # ── Clone sizes per (mouse × CDR3) ───────────────────────────────────────
    clone_df = (data.groupby([mouse_col, CDR3_COL])
                    .size()
                    .reset_index(name='clone_size'))

    # ── Align unique CDR3s ────────────────────────────────────────────────────
    unique_seqs = clone_df[CDR3_COL].unique().tolist()
    try:
        aln_map = align_sequences(unique_seqs, method=ALIGN)
        print(f"[{exp['label']}] aligned {len(unique_seqs)} unique CDR3s via {ALIGN}")
    except Exception as e:
        print(f"[{exp['label']}] {ALIGN} failed ({e}), falling back to center")
        aln_map = center_align(unique_seqs)

    clone_df['cdr3_aln'] = clone_df[CDR3_COL].map(aln_map)
    clone_df = clone_df.dropna(subset=['cdr3_aln']).copy()

    # ── Consensus ─────────────────────────────────────────────────────────────
    cons  = consensus_sequence(clone_df['cdr3_aln'].tolist())
    max_k = len(cons)
    print(f"[{exp['label']}] consensus ({max_k} aa): {cons}")

    # ── Position-wise mismatch matrix  (n_clones × max_k) ────────────────────
    aln_seqs = clone_df['cdr3_aln'].values
    mm = np.array([
        [int(s[i] != cons[i]) if i < len(s) else 1 for i in range(max_k)]
        for s in aln_seqs
    ], dtype=np.int8)
    dist_cumsum = mm.cumsum(axis=1)          # dist_cumsum[:, k-1] = Hamming at prefix k

    log_sizes   = np.log(clone_df['clone_size'].values.astype(float))
    mice        = clone_df[mouse_col].unique()
    k_vals      = np.arange(1, max_k + 1)
    r_matrix    = np.full((len(mice), max_k), np.nan)

    # ── Per-mouse Spearman r vs. k ────────────────────────────────────────────
    for mi, mouse in enumerate(mice):
        mask = clone_df[mouse_col].values == mouse
        ls_m = log_sizes[mask]
        dc_m = dist_cumsum[mask, :]

        for k in range(1, max_k + 1):
            d = dc_m[:, k - 1]
            if d.std() < 1e-10 or mask.sum() < 3:
                continue
            r, p = stats.spearmanr(ls_m, d)
            r_matrix[mi, k - 1] = r
            all_records.append({
                'experiment': exp['label'], 'mouse': str(mouse),
                'k': k, 'spearman_r': r, 'p_value': p,
                'n_clones': int(mask.sum()),
            })

    # ── Plot: correlation vs. k ───────────────────────────────────────────────
    mean_r = np.nanmean(r_matrix, axis=0)
    std_r  = np.nanstd( r_matrix, axis=0)

    ax_corr.plot(k_vals, mean_r, color=color, lw=2.2)
    ax_corr.fill_between(k_vals, mean_r - std_r, mean_r + std_r,
                         color=color, alpha=0.20)
    ax_corr.axhline(0, color='k', lw=0.8, ls='--', alpha=0.4)
    ax_corr.set_xlabel(r'Alignment prefix $k$', fontsize=11)
    ax_corr.set_ylabel(r'Spearman $r$', fontsize=11)
    ax_corr.set_title(exp['label'], fontsize=10)
    ax_corr.tick_params(labelsize=9)

    # ── Plot: log(size) vs. distance at full k ────────────────────────────────
    dist_full = dist_cumsum[:, -1]

    ax_sc.scatter(dist_full, log_sizes, color=color, alpha=0.20,
                  s=10, rasterized=True, zorder=2)

    # Binned mean ± SD
    d_min, d_max = int(dist_full.min()), int(dist_full.max())
    if d_max > d_min:
        bins = np.arange(d_min, d_max + 2)
        bin_mean, edges, _ = stats.binned_statistic(
            dist_full, log_sizes, statistic='mean', bins=bins)
        bin_std, _, _  = stats.binned_statistic(
            dist_full, log_sizes, statistic='std',  bins=bins)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax_sc.errorbar(centers, bin_mean, yerr=bin_std,
                       color=color, fmt='o-', capsize=3, ms=5, lw=1.8, zorder=5)

    ax_sc.set_xlabel(r'Hamming distance to consensus', fontsize=11)
    ax_sc.set_ylabel(r'$\log$ clone size', fontsize=11)
    ax_sc.set_title(exp['label'], fontsize=10)
    ax_sc.tick_params(labelsize=9)


# ── Save ──────────────────────────────────────────────────────────────────────

fig_corr.suptitle(
    rf'Spearman $r(\log\,\mathrm{{size}},\,d_{{\mathrm{{cons}}}})$ vs.\ alignment prefix — {ALIGN}',
    fontsize=12)
fig_corr.savefig(output_plot + f'/corr_vs_length_{ALIGN}.pdf', bbox_inches='tight')

fig_sc.suptitle(
    rf'$\log$ clone size vs.\ Hamming distance to CDR3 consensus — {ALIGN}',
    fontsize=12)
fig_sc.savefig(output_plot + f'/size_vs_distance_{ALIGN}.pdf', bbox_inches='tight')

pd.DataFrame(all_records).to_csv(
    output_plot + f'/corr_results_{ALIGN}.csv', index=False)

print("Done.")
