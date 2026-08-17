"""
drift_scan.py -- minimal antigenic-drift crossover scan (v4, faster).

Runs the primary once, forms the memory seed once, then for each drift Delta
rigidly shifts that seed toward higher E (CONSERVING CELL COUNTS) and runs only
the recall. Compares recall potency to the primary at a frozen clock t_c.

Speed changes vs v3:
  * T_SIM = 5.0 (was 10.0) -- still well past t_c ~ 3.5
  * DELTAS coarsened to 0..6 step 1 (refine near the crossover afterwards)
  * max_step cap on recall solves + per-Delta progress print

No edits to run_meanfield.py needed. Run:  python drift_scan.py
"""

import os, sys, time
import numpy as np
import matplotlib.pyplot as plt

# --- point this at wherever mf_lib.py lives ---
sys.path.append(os.path.join(os.path.dirname(__file__), '../../library'))
import mf_lib as mf
import mf_plotting as mfp

# Reuse the exact BASE parameters from the driver so the physics matches.
import run_meanfield as R


# ---- scan knobs ----
DELTAS = np.arange(0.0, 4.01, 0.5)   # coarse first pass; refine near crossover
R_MATURATION = 0.75           # serum fragility factor r (<1); external input
NU_NSTAR_ANALYTIC = None     # set to your analytic nu*N* to print the comparison
T_SIM = 5.0                  # only need to clear t_c (~3.5); 5 gives margin
MODE = 'grid'
MAX_STEP = 0.05              # cap LSODA micro-steps on the stiff recall solves

CONFIG = dict(
    # --- output / figures ---
    figures=['NA_shared', 'Z_shared', 'pb_shared', 'DG_shared', 'n0_shared', 'fig4d', 'P_pi', 'per_run'],
    #   available: 'NA_shared', 'Z_shared', 'pb_shared', 'per_run', 'h0_extras'
    outroot='/Users/robertomorantovar/Library/CloudStorage/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/antigenic_evolution/meanfield/',                  # base directory for saved PDFs
    usetex=True,                       # set True if you have a LaTeX toolchain
    show=False,                         # plt.show() at the end instead of/along saving
)

def shift_seed(seed, delta, p, weights):
    """Rigid drift E -> E + delta with CELL-COUNT conservation.

    seed is per-clone-normalized: seed[i]*weights[i] = real cells at index i.
    We roll the REAL cells (so Z's re-multiplication by weights[i] is honest),
    then convert back to per-clone at the new index.
    """
    if delta <= 0:
        return seed
    dDG = (p.DG_max - p.DG_min) / (p.M - 1)
    k = int(round(delta / dDG))
    if k <= 0:
        return seed
    cells = seed * weights                                  # real cells per bin
    cells_shifted = np.concatenate([np.zeros(k), cells[:-k]])
    with np.errstate(divide='ignore', invalid='ignore'):
        seed_shifted = np.where(weights > 0, cells_shifted / weights, 0.0)
    return seed_shifted


def make_params(memory):
    """Fresh Parameters from the driver's BASE, with the memory flag set.
    Memory phase forces h0 = b0 (the driver's memory_h0_equals_b0 behaviour)."""
    d = dict(R.BASE)
    d['memory'] = memory
    if memory == 1:
        d['h0'] = d['b0']
    return mf.Parameters(**d)


def run_recall(p_rec, seed):
    return mf.run_simulation('semicomplete', p=p_rec, t_span=(0.0, T_SIM),
                             mode=MODE, memory_seed=seed, max_step=MAX_STEP)


def main(cfg=CONFIG):
    t_start = time.time()

    # 1. primary, once (drift-invariant reference)
    p_prim = make_params(memory=0)
    r_prim = mf.run_simulation('semicomplete', p=p_prim, t_span=(0.0, T_SIM),
                               mode=MODE, max_step=MAX_STEP)

    # 2. form the memory seed once
    _, seed0 = mf.memory_seed_from_primary(r_prim, p=p_prim,
                                           n_mem=int(p_prim.n_mem))
    if seed0.sum() <= 0:
        raise RuntimeError("Primary produced no memory seed; check parameters.")
    ln_headstart_sim = np.log(seed0.sum())

    # recall parameters + grid weights (weights needed by shift_seed)
    p_rec = make_params(memory=1)
    ctx_rec = mf.build_context(p_rec, mode=MODE)
    weights_rec = ctx_rec['weights']

    # 3+4. frozen clock from the UNDRIFTED recall; primary reference potency
    print("running undrifted recall (sets frozen clock)...", flush=True)
    r_rec0 = run_recall(p_rec, seed0.copy())
    _, t_c_fixed = mf.find_t_c(r_rec0)
    if not np.isfinite(t_c_fixed):
        raise RuntimeError("Undrifted recall never crosses Z_c; "
                           "lower Z_c or raise T_SIM before scanning.")
    P_prim = mf.compute_potency_P(r_prim, t_c_star=t_c_fixed)[1]
    print(f"  t_c = {t_c_fixed:.3f},  P_primary = {P_prim:.4g}", flush=True)

    # 5. sweep Delta: shift seed (cell-conserving), run recall, potency at t_c
    P_rec = []
    for d in DELTAS:
        t0 = time.time()
        seed_d = shift_seed(seed0.copy(), d, p_rec, weights_rec)
        r_rec = run_recall(p_rec, seed_d)
        Zt = mf.compute_potency_t(r_rec)
        mask = r_rec['t'] > t_c_fixed
        P = Zt[mask][0] if mask.any() else np.nan
        P_rec.append(P)
        print(f"  Delta={d:4.1f}  P_rec={P:.4g}  ({time.time()-t0:.1f}s)",
              flush=True)
    P_rec = np.array(P_rec)

    # 6. crossover Delta* (first downward crossing of the primary line)
    below = np.where(P_rec <= P_prim)[0]
    if len(below) == 0:
        delta_star = np.nan
        print("[warn] recall never drops to primary; extend DELTAS.")
    elif below[0] == 0:
        delta_star = 0.0
        print("[warn] recall <= primary already at Delta=0; check setup.")
    else:
        i = below[0]
        d0, d1 = DELTAS[i - 1], DELTAS[i]
        y0, y1 = P_rec[i - 1], P_rec[i]
        delta_star = d0 + (P_prim - y0) * (d1 - d0) / (y1 - y0)

    # 7. titer conversion
    fold_bare = delta_star / np.log(2) if np.isfinite(delta_star) else np.nan
    fold_serum = fold_bare / R_MATURATION if np.isfinite(fold_bare) else np.nan

    # optional: fit the falling branch to exp(-c*Delta) and print c
    c_fit = np.nan
    if np.isfinite(delta_star):
        i_pk = int(np.nanargmax(P_rec))
        xs, ys = DELTAS[i_pk:], P_rec[i_pk:]
        good = ys > 0
        if good.sum() >= 2:
            c_fit = -np.polyfit(xs[good], np.log(ys[good]), 1)[0]

    # report
    print("\n===== drift crossover =====")
    print(f"lambda_B/lambda_A      = {p_prim.b0/p_prim.lambda_A:.3f}  "
          f"(analytic band 0.15-0.6)")
    print(f"t_c (frozen clock)     = {t_c_fixed:.3f}")
    print(f"P_primary (fixed)      = {P_prim:.4g}")
    print(f"Delta*  (energy units) = {delta_star:.3f}")
    print(f"titer folds (bare)     = {fold_bare:.2f}")
    print(f"titer folds (r={R_MATURATION}) = {fold_serum:.2f}")
    if np.isfinite(c_fit):
        print(f"falling-branch slope c = {c_fit:.3f}  "
              f"(compare to 1+lambda_z/lambda_A)")
    print("\n----- head-start consistency -----")
    print(f"ln(seed.sum())  [sim]  = {ln_headstart_sim:.3f}")
    if NU_NSTAR_ANALYTIC is not None:
        print(f"ln(nu N*) [analytic]   = {np.log(NU_NSTAR_ANALYTIC):.3f}")
        print(f"difference             = "
              f"{ln_headstart_sim - np.log(NU_NSTAR_ANALYTIC):+.3f}")

    print(f"\ntotal wall time: {time.time()-t_start:.1f}s")

    # plot
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.plot(DELTAS, P_rec/P_prim, 'o-', color='C0', lw=2, label=r'recall $P(\Delta)$')
    ax.axhline(1, color='0.4', ls='--', lw=2, label='primary (reference)')
    if np.isfinite(delta_star) and delta_star > 0:
        ax.axvline(delta_star, color='k', ls=':', lw=1.5)
        ax.annotate(rf'$\Delta^*={delta_star:.2f}$',
                    xy=(delta_star, 1), xytext=(6, 8),
                    textcoords='offset points', fontsize=11)
    ax.set_xlabel(r'antigenic drift  $\Delta$  (energy units)')
    ax.set_ylabel(r'potency  $P$  at fixed $t_c$')
    ax.set_yscale('log')
    ax.legend(frameon=False)
    secax = ax.secondary_xaxis('top',
                           functions=(lambda d: d / (R_MATURATION * np.log(2)),
                                      lambda f: f * R_MATURATION * np.log(2)))
    secax.set_xlabel(f'serum titer folds  ($r={R_MATURATION}$)')
    out = os.path.join(cfg['outroot'], 'exploration')
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, 'drift_crossover.png'), dpi=150)
    print(f"saved: {out}")

if __name__ == '__main__':
    main()