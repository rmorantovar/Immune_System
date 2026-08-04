"""
run_meanfield.py -- one config-driven driver for the mean-field models.

Replaces mean-field_null.py, mean-field_semicomplete.py and mf_potency_h0.py.
Everything you used to change by editing a different script is now a field in
CONFIG below:

    * models             : which model(s) to run  -> ['semicomplete'] or
                           ['semicomplete', 'null'] to overlay a null comparison.
    * memory_phases      : [0] primary only, [0, 1] primary + memory.
    * sweep              : (param_name, [values])  -> iterate over any Parameter.
    * figures            : which figures to produce.

Reproducing the three original scripts:

    mean-field_semicomplete.py :
        models=['semicomplete'], sweep=('h0', [b0/1000.]), memory_phases=[0, 1]
    mean-field_null.py :
        models=['semicomplete', 'null'], sweep=('h0', [b0/1000.]),
        memory_phases=[0, 1]
    mf_potency_h0.py :
        models=['semicomplete'], memory_phases=[0, 1],
        sweep=('h0', flip(logspace(log10(b0/1e5), log10(b0/1e0), 6))),
        figures=[..., 'h0_extras']
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../lib'))

import numpy as np
import matplotlib.pyplot as plt

import mf_lib as mf
import mf_plotting as mfp


# ============================================================
# CONFIG  -- edit this block, then run `python run_meanfield.py`
# ============================================================
_b0 = 1.5

BASE = dict(
    N_A0=1.0, lambda_innate=2.2, threshold_innate=5e3,
    lambda_A=6.2, delta_A=3.0, eta=1.0,
    k_on=1e2 * 2e5 * 1e6 * 24 * 3600 / mf.N_Avg, delta_pi=0.1,
    hill=2.0, beta_star=2.3, K_T=1e4,
    delta_T=0.00, Tcell_growth_factor=2.0,
    tau_eng=0.1, b0=_b0, delta_B=0.00, h0=_b0 / 1000.,
    DG_min=0.0, DG_max=8.0, M=32,
    omega_0=1.0, T_lim=True, N_T0=1e6,
    Z_c=1e3, n_mem=1e5,
)

CONFIG = dict(
    # --- what to run ---
    models=['semicomplete'],            # subset of mf.MODELS
    memory_phases=[0, 1],               # 0 = primary, 1 = memory (seeded from 0)
    sweep=('h0', np.flip(np.logspace(np.log10(_b0/1e3), np.log10(_b0/1e2), 2))),                         # None -> no sweep (use BASE params as-is).
    # sweep = None,
                                        # To sweep: ('h0', [v1, v2, ...]) or any
                                        # Parameter name with a list/array of values.
    sweep_cmap='summer',                # colormap for multi-value sweeps
                                        #   ('summer' = your green scale; None = discrete)
    T=10.0,
    mode='grid',                        # 'grid' or 'stochastic'
    base=BASE,

    # --- behaviour toggles ---
    memory_h0_equals_b0=True,           # memory phase forces h0 = b0 (original behaviour)

    # --- output / figures ---
    figures=['NA_shared', 'Z_shared', 'pb_shared', 'per_run'],
    #   available: 'NA_shared', 'Z_shared', 'pb_shared', 'per_run', 'h0_extras'
    outroot='/Users/robertomorantovar/Library/CloudStorage/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/exponential_proofreading/meanfield/figures',                  # base directory for saved PDFs
    usetex=True,                       # set True if you have a LaTeX toolchain
    show=False,                         # plt.show() at the end instead of/along saving
)


# ============================================================
# Parameter construction (fixes the in-place-mutation + p-before-assignment bugs)
# ============================================================
# Defaults filled in for any key omitted from a user config.
CONFIG_DEFAULTS = dict(
    models=['semicomplete'], memory_phases=[0], sweep=None, sweep_cmap='summer',
    T=10.0, mode='grid', base=BASE, memory_h0_equals_b0=True,
    figures=['NA_shared', 'Z_shared', 'pb_shared', 'per_run'],
    outroot='figures', usetex=True, show=False,
)


def normalize_config(cfg):
    """Fill in defaults and list-ify sweep/models/memory_phases.

    Any key omitted from `cfg` falls back to CONFIG_DEFAULTS, so a minimal
    config like dict(models=['null','semicomplete'], memory_phases=[0,1]) works.
    Also lets you pass numpy arrays (e.g. sweep=('h0', np.logspace(...))) without
    hitting `'numpy.ndarray' object has no attribute 'index'`.
    """
    cfg = {**CONFIG_DEFAULTS, **cfg}
    sweep = cfg.get('sweep')
    if sweep is None:                       # no sweep -> single point with BASE
        cfg['sweep'] = (None, [None])
    else:
        name, values = sweep
        cfg['sweep'] = (name, list(values))
    cfg['models'] = list(cfg['models'])
    cfg['memory_phases'] = list(cfg['memory_phases'])
    return cfg


def build_params(cfg, sweep_value, memory):
    """Return a fresh Parameters for one (sweep_value, memory) combination.

    Never mutates cfg['base']; the memory h0 override is applied on a copy.
    """
    d = dict(cfg['base'])
    sweep_name, _ = cfg['sweep']
    if sweep_name is not None:
        d[sweep_name] = sweep_value
    d['memory'] = memory
    if memory == 1 and cfg.get('memory_h0_equals_b0', True):
        d['h0'] = d['b0']
    return mf.Parameters(**d)


def pi_star_of(p):
    return (p.b0 / p.h0) ** (1 / p.hill)

def label_pi_star(res):
    """pi_c label for a run. For memory runs this uses the *naive* run's params
    (stored in res['params_primary']), so the label reports the primary pi_c
    rather than the memory run's forced value."""
    p = res.get('params_primary', res['params'])
    return mfp.pi_star_label(pi_star_of(p))

# ============================================================
# Experiment engine
# ============================================================
def run_experiment(cfg):
    """
    Loop models x sweep x memory, seeding each memory run from its own primary.

    Returns
    -------
    results : dict keyed by (model, sweep_value, memory) -> res dict.
    """
    cfg = normalize_config(cfg)
    sweep_name, sweep_values = cfg['sweep']
    results = {}
    memory_seeds = {}      # (model, sweep_value) -> formed_memory from the primary
    primary_params = {}    # (model, sweep_value) -> Parameters of the naive run

    print("Running experiment")
    for model in cfg['models']:
        print(f"  model = {model}")
        for sval in sweep_values:
            if sweep_name is None:
                print("    (fixed parameters, no sweep)")
            else:
                print(f"    {sweep_name} = {sval:.4g}")
            for memory in cfg['memory_phases']:
                p = build_params(cfg, sval, memory)
                seed = memory_seeds.get((model, sval)) if memory == 1 else None
                if memory == 1 and seed is None:
                    print("      [skip] memory phase requested but no primary "
                          "(add 0 to memory_phases)")
                    continue

                res = mf.run_simulation(model, p=p, t_span=(0.0, cfg['T']),
                                        mode=cfg['mode'], memory_seed=seed)
                # Attach the naive run's parameters so memory plots can recover
                # the primary pi_c (which differs from the memory run's own p).
                res['params_primary'] = primary_params.get((model, sval), p)
                results[(model, sval, memory)] = res

                if memory == 0:
                    primary_params[(model, sval)] = p
                    try:
                        _, formed = mf.memory_seed_from_primary(
                            res, p=p, n_mem=int(p.n_mem))
                        memory_seeds[(model, sval)] = formed
                    except ValueError as e:
                        print(f"      [warn] no memory seed for {model}: {e}")
    return results


# ============================================================
# Plotting
# ============================================================
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


def _sweep_tag(cfg, sval):
    """Folder-name fragment for one sweep value ('fixed' when not sweeping)."""
    name = cfg['sweep'][0]
    return 'fixed' if name is None else f'{name}_{sval:.4g}'


def make_shared_figures(cfg, results):
    """Overlay figures across all runs: antigen, potency, neutralization pb."""
    sweep_name, sweep_values = cfg['sweep']
    T = cfg['T']
    figs = cfg['figures']

    fig_NA, ax_NA = mfp.new_fig(wide=True)
    fig_Z, ax_Z = mfp.new_fig()
    fig_pb, ax_pb = mfp.new_fig(wide=True)

    # colour scale across the sweep (e.g. the green 'summer' scale for an h0 sweep)
    sweep_cols = mfp.sweep_colors(len(sweep_values), cfg.get('sweep_cmap', 'summer'))

    for (model, sval, memory), res in results.items():
        i_s = sweep_values.index(sval)
        i_model = cfg['models'].index(model)
        # multi-value sweep -> colour by sweep value; single value -> by model
        color = (sweep_cols[i_s] if len(sweep_values) > 1
                 else mfp.colors_sim[i_model])
        ls = mfp.styles_sim[memory]
        label = label_pi_star(res)
        if memory == 1:
            print(label)
        # antigen
        na_color = mfp.antigen_color if memory == 0 else color
        mfp.plot_NA(ax_NA, res, color=na_color, ls=ls, label=label)

        # potency + Z_c crossing
        Z = mfp.plot_potency(
            ax_Z, res, color=color, ls=ls, label=(label if memory == 1 else None),
            marker_face=('white' if memory == 0 else color),
            marker_edge=(color if memory == 0 else 'k'))

        # neutralization probability (and the external innate drive for primary)
        if memory == 0:
            mfp.plot_innate(ax_pb, res, color=color, ls='--')
        mfp.plot_pb_from_potency(ax_pb, res, color=color, ls=ls, label=label)

    out = _outdir(cfg, 'shared')

    if 'NA_shared' in figs:
        mfp.style_log_axis(ax_NA, T, ylim=(1e0, 1e11), xlim=(0.5, T - 3))
        fig_NA.savefig(os.path.join(out, 'N_A_shared.pdf'), dpi=200)

    if 'Z_shared' in figs:
        Z_c = cfg['base']['Z_c']
        ax_Z.axhline(Z_c, lw=1, ls='--', color='k')
        
        if sweep_name == 'h0':
            mfp.style_log_axis(ax_Z, T, ylim=(5e-1, 1e5), xlim=(2.0, T - 5.5), hide_xticklabels=False)
            ax_Z.legend(fontsize=22, title=r'$\pi_c$', loc=0, title_fontsize=24)
        else:
            mfp.style_log_axis(ax_Z, T, ylim=(5e-1, 1e5), xlim=(1.5, T - 3))

        fig_Z.savefig(os.path.join(out, 'Z_shared.pdf'), dpi=200, transparent=True)

    if 'pb_shared' in figs:
        ax_pb.set_xlim(0.5, T - 3)
        ax_pb.tick_params(axis='y', labelsize=30)
        ax_pb.tick_params(axis='x', labelsize=30)
        fig_pb.savefig(os.path.join(out, 'pb_shared.pdf'), dpi=200, transparent=True)

    if not cfg['show']:
        plt.close(fig_NA); plt.close(fig_Z); plt.close(fig_pb)


def make_per_run_figures(cfg, results):
    """Standard per-run panels (antigen, pi, N_B, potency)."""
    T = cfg['T']
    for (model, sval, memory), res in results.items():
        color = mfp.colors_sim[cfg['memory_phases'].index(memory)]
        label = label_pi_star(res)
        out = _outdir(cfg, model, _sweep_tag(cfg, sval), f'memory_{memory}')

        fig_NA, ax_NA = mfp.new_fig(wide=True)
        mfp.plot_NA(ax_NA, res, label=label)
        mfp.style_log_axis(ax_NA, T, ylim=(1e0, 1e11))
        fig_NA.savefig(os.path.join(out, 'N_A.pdf'), dpi=150)

        if 'pi' in res:
            fig_pi, ax_pi = mfp.new_fig(wide=True)
            mfp.plot_pi(ax_pi, res, color=color, label=label)
            ax_pi.axhline(pi_star_of(res['params']), lw=1, ls='--',
                          color='grey', alpha=0.8)
            mfp.style_log_axis(ax_pi, T, ylim=(5e-1, 1e4))
            fig_pi.savefig(os.path.join(out, 'pi.pdf'), dpi=150)

        fig_NB, ax_NB = mfp.new_fig(wide=True)
        mfp.plot_NB(ax_NB, res, color=color, label=label)
        ax_NB.axhline(1.0, color='k', ls='--', alpha=0.5)
        mfp.style_log_axis(ax_NB, T, ylim=(5e-1, 1e6))
        fig_NB.savefig(os.path.join(out, 'N_B.pdf'), dpi=150)

        fig_Z, ax_Z = mfp.new_fig(wide=True)
        mfp.plot_potency(ax_Z, res, color=color, label=label)
        ax_Z.axhline(res['params'].Z_c, lw=1, ls='--', color='k')
        mfp.style_log_axis(ax_Z, T, ylim=(5e-3, 1e7))
        ax_Z.set_xlabel('Time', fontsize=16)
        fig_Z.savefig(os.path.join(out, 'Z.pdf'), dpi=150)

        if not cfg['show']:
            plt.close('all')


def make_h0_extra_figures(cfg, results):
    """
    The extra panels from mf_potency_h0.py: the memory-seed distribution n0(DG)
    and the analytic affinity front DG(t). Only meaningful for an h0 sweep.
    """
    sweep_name, sweep_values = cfg['sweep']
    T = cfg['T']
    fig_n0, ax_n0 = mfp.new_fig()
    fig_DG, ax_DG = mfp.new_fig(wide=True)

    colors_h0 = mfp.sweep_colors(len(sweep_values), cfg.get('sweep_cmap', 'summer'))
    for (model, sval, memory), res in results.items():
        i_s = sweep_values.index(sval)
        label = mfp.pi_star_label(pi_star_of(res['params']))

        if memory == 0:
            # affinity front (growth-phase slope 1.0 as in the h0 script)
            fr = mfp.compute_affinity_front(res, on_slope=1.0)
            ax_DG.plot(fr['ton'], fr['DG_on'], color=colors_h0[i_s],
                       ls='dotted', label=label)
            off = fr['DG_off'] > res['params'].DG_min
            ax_DG.plot(fr['toff'][off], fr['DG_off'][off],
                       color=colors_h0[i_s], ls='dashed', lw=3)
            # per-clone memory pool that would be seeded from this primary
            expanded, _ = mf.memory_seed_from_primary(
                res, p=res['params'], n_mem=int(res['params'].n_mem))
            ax_n0.plot(res['DG'], expanded, color=colors_h0[i_s], lw=2, label=label)

    out = _outdir(cfg, 'shared')
    ax_n0.set_yscale('log')
    ax_n0.set_ylim(1, 2.5e9)
    ax_n0.legend(fontsize=16, title=r'$\pi_c$')
    fig_n0.savefig(os.path.join(out, 'n0_shared.pdf'), dpi=200, transparent=True)

    ax_DG.set_ylim(bottom=cfg['base']['DG_min'] - 0.1)
    ax_DG.set_xlim(0, T)
    ax_DG.legend(fontsize=14)
    fig_DG.savefig(os.path.join(out, 'DG_shared.pdf'), dpi=150)

    if not cfg['show']:
        plt.close(fig_n0); plt.close(fig_DG)


# ============================================================
# Entry point
# ============================================================
def main(cfg=CONFIG):
    cfg = normalize_config(cfg)
    mfp.use_style(usetex=cfg['usetex'])
    results = run_experiment(cfg)

    figs = cfg['figures']
    if any(f in figs for f in ('NA_shared', 'Z_shared', 'pb_shared')):
        make_shared_figures(cfg, results)
    if 'per_run' in figs:
        make_per_run_figures(cfg, results)
    if 'h0_extras' in figs:
        make_h0_extra_figures(cfg, results)

    if cfg['show']:
        plt.show()
    print("Done.")
    return results


if __name__ == '__main__':
    main()