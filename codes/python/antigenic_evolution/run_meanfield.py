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
sys.path.append(os.path.join(os.path.dirname(__file__), '../../library'))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

import mf_lib as mf
import mf_plotting as mfp


# ============================================================
# CONFIG  -- edit this block, then run `python run_meanfield.py`
# ============================================================
_b0 = 1.5
# Fig 4D analytic-D parameters (naive N(0, sigma2); edge DG_star = -beta* sigma2)
DG4D_STAR = -16.0        # high-affinity edge (negative)
DG4D_SIG2 = 7.0          # naive variance sigma^2

BASE = dict(
    N_A0=1.0, lambda_innate=2.2, threshold_innate=5e3,
    lambda_A=6.2, delta_A=3.0, eta=1.0,
    k_on=1e2 * 2e5 * 1e6 * 24 * 3600 / mf.N_Avg, delta_pi=0.1,
    hill=2.0, beta_star=2.3, K_T=1e4,
    delta_T=0.00, Tcell_growth_factor=2.0,
    tau_eng=0.1, b0=_b0, delta_B=0.00, h0=_b0 / 1000.,
    DG_min=0.0, DG_max=12.0, M=40,
    omega_0=1.0, T_lim=True, N_T0=1e6,
    Z_c=1.4e3, n_mem=4e4,
)

_alpha_o = BASE['b0'] / BASE['lambda_A'] + BASE['b0'] / BASE['delta_A']   # null EP exponent
FIG4D_MU  = {'primary': -8.0,    # mu = -zeta*|DG*|,  zeta ~ 0.5  (tilted once)
             'recall':  -16.0,   # mu ~ DG*,          zeta ~ 1    (tilted twice)
             'null': -5.16
            # 'null':    -(_alpha_o / BASE['beta_star']) * abs(DG4D_STAR),          # ~ -5.16
            }

CONFIG = dict(
    # --- what to run ---
    models=['semicomplete', 'null'],            # subset of mf.MODELS
    memory_phases=[0, 1],               # 0 = primary, 1 = memory (seeded from 0)
    # sweep=('h0', np.flip(np.logspace(np.log10(_b0/1e5), np.log10(_b0/1e0), 6))),                         # None -> no sweep (use BASE params as-is).
    # sweep=('h0', np.flip([_b0/1e5, _b0/1e3, _b0/1e1])),                         # None -> no sweep (use BASE params as-is).
    sweep = None,
                                        # To sweep: ('h0', [v1, v2, ...]) or any
                                        # Parameter name with a list/array of values.
    null_sval= 1.,      # single value: used to RUN null and to PLOT the null point                               
    sweep_cmap='autumn',                # colormap for multi-value sweeps
                                        #   ('summer' = your green scale; None = discrete)
    T=10.0,
    mode='grid',                        # 'grid' or 'stochastic'
    base=BASE,
    drift_delta=0.0,

    # --- behaviour toggles ---
    memory_h0_equals_b0=True,           # memory phase forces h0 = b0 (original behaviour)

    # --- output / figures ---
    figures=['NA_shared', 'Z_shared', 'pb_shared', 'DG_shared', 'n0_shared', 'fig4d', 'P_pi', 'per_run'],
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
    T=12.0, mode='grid', base=BASE, memory_h0_equals_b0=True,
    figures=['NA_shared', 'Z_shared', 'pb_shared', 'per_run'],
    outroot='figures', usetex=True, show=False, drift_delta=0.0,
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
        if model == 'null' and cfg.get('null_sval') is not None:
            model_vals = [cfg['null_sval']]                       # run null only here
        else:
            model_vals = cfg.get('model_sweep_values', {}).get(model, sweep_values)

        for sval in model_vals:
            if sweep_name is None:
                print("    (fixed parameters, no sweep)")
            else:
                print(f"    {sweep_name} = {sval:.4g}")
            for memory in cfg['memory_phases']:
                p = build_params(cfg, sval, memory)
                seed = memory_seeds.get((model, sval)) if memory == 1 else None
                delta = cfg.get('drift_delta', 0.0)
                if memory == 1 and seed is not None and delta > 0:
                    dDG = (cfg['base']['DG_max'] - cfg['base']['DG_min']) / (cfg['base']['M'] - 1)
                    k = int(round(delta / dDG))
                    if k > 0:
                        seed = np.concatenate([np.zeros(k), seed[:-k]])   # push counts to higher E
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

def _sweep_tag(cfg, sval):
    """Folder-name fragment for one sweep value ('fixed' when not sweeping)."""
    name = cfg['sweep'][0]
    return 'fixed' if name is None else f'{name}_{sval:.4g}'


def make_shared_figures(cfg, results):
    """Overlay figures across all runs: antigen, potency, pb, Fig4D, DG front, n0."""
    sweep_name, sweep_values = cfg['sweep']
    T = cfg['T']
    figs = cfg['figures']

    fig_NA, ax_NA = mfp.new_fig()
    fig_Z,  ax_Z  = mfp.new_fig()
    fig_pb, ax_pb = mfp.new_fig()
    fig_4d, ax_4d = mfp.new_fig()
    fig_DG, ax_DG = mfp.new_fig()
    fig_n0, ax_n0 = mfp.new_fig()
    fig_n02, ax_n02 = mfp.new_fig()
    fig_P_pi, ax_P_pi = mfp.new_fig()

    all_figs = [fig_NA, fig_Z, fig_pb, fig_4d, fig_DG, fig_n0, fig_n02, fig_P_pi]

    # sweep_cols_primary = mfp.sweep_colors(len(sweep_values), 'autumn')   # memory 0
    sweep_cols_primary =  mfp.colors_sweep   # memory 0
    # sweep_cols_memory  = mfp.sweep_colors(len(sweep_values), 'winter')   # memory 1
    sweep_cols_memory  = mfp.colors_sweep   # memory 1

    fig4d_pts = {'primary': [], 'recall': []}
    fig4d_map = {('semicomplete', 0): 'primary', ('semicomplete', 1): 'recall'}

    opt_sval = None
    if len(sweep_values) > 1:
        t_c_star = np.inf
        for sv in sweep_values:
            r = results.get(('semicomplete', sv, 1))
            if r is None:
                continue
            Z = mf.compute_potency_t(r)
            t_c = r['t'][Z >= r['params'].Z_c][0] if np.any(Z >= r['params'].Z_c) else np.inf   
            # _, P = mf.compute_potency_P(r)
            # if not np.isnan(P) and P > _best:
            if not np.isnan(t_c) and t_c < t_c_star:
                t_c_star, opt_sval = t_c, sv

    for (model, sval, memory), res in results.items():
        # colour: by pi_c when sweeping; else brown (null) / green-blue (single run)
        if len(sweep_values) > 1 and sval in sweep_values:
            pi_c = (_b0 / sval) ** (1 / res['params'].hill)  # for debugging: check that the label matches the run's pi_c
            if model == 'null':
                color, ls = mfp.my_grey, mfp.styles_mem[memory]
            else:
                i_s = sweep_values.index(sval)
                color = sweep_cols_primary[i_s] if memory == 0 else sweep_cols_memory[i_s]
                ls = mfp.styles_mem[memory]
        else:
            pi_c = (_b0 / res['params'].h0) ** (1 / res['params'].hill)
            if model == 'null':
                color, ls = mfp.my_grey, mfp.styles_mem[memory]
            else:
                color, ls = mfp.colors_mem[memory],  mfp.styles_mem[memory]
        label = label_pi_star(res)

        na_color = mfp.antigen_color if memory == 0 else color
        mfp.plot_NA(ax_NA, res, color=na_color, ls=ls, label=label)

        if len(sweep_values) > 1 and sval in sweep_values:
            if memory == 0:
                mfp.plot_innate(ax_pb, res, color='grey', ls='--')
                if opt_sval is None or sval == opt_sval or model == 'null':
                    mfp.plot_potency(ax_Z, res, color=color, ls=ls,
                            marker_face=color, marker_edge='k',
                            marker='o', mark_Zc=False, mark_t = t_c_star)
                    t_peak = res['t'][np.argmax(res['N_A'])]
                    ax_Z.axvline(t_peak, lw=2, ls='--', color='grey')

            else:
                mfp.plot_potency(ax_Z, res, color=color, ls=ls,
                    label=label,
                    marker_face=color, marker_edge='k',
                    marker='s', mark_Zc=False, mark_t = t_c_star)
                if opt_sval is None or sval == opt_sval:
                    ax_Z.axvline(t_c_star, lw=2, ls='--', color='k')

        mfp.plot_pb_from_potency(ax_pb, res, color=color, ls=ls, label=label)

        if (model, memory) in fig4d_map:
            name = fig4d_map[(model, memory)]
            D, (_, P)  = mf.compute_specificity(res), mf.compute_potency_P(res, t_c_star=t_c_star)
            D_ana = mf.dkl_analytic(FIG4D_MU[name], DG4D_STAR + mf.DG_mf(res),
                                    DG_star=DG4D_STAR + mf.DG_min(res), sigma2=DG4D_SIG2)
            if not (np.isnan(D) or np.isnan(P)):
                if memory == 0:
                    if sval == opt_sval:
                        fig4d_pts[name].append((D, P, color, D_ana, pi_c))
                else:
                    fig4d_pts[name].append((D, P, color, D_ana, pi_c))
            ax_P_pi.scatter(pi_c, -P, color=color, ls=ls, marker='o', s=64, label=label)

        if memory == 0 and model == 'semicomplete' and sval in sweep_values:
            fr = mfp.compute_affinity_front(res, on_slope=1.0)
            ax_DG.plot(fr['ton'], fr['DG_on'], color=color, ls='--', label=label, lw=3)
            off = fr['DG_off'] > res['params'].DG_min
            ax_DG.plot(fr['toff'][off], fr['DG_off'][off], color=color, ls='--', lw=3)
            expanded, produced = mf.memory_seed_from_primary(res, p=res['params'],
                                                      n_mem=int(res['params'].n_mem))
            ax_n0.plot(res['DG'], expanded, color=color, lw=2, ls = '--')
            ax_n0.plot(res['DG'], produced * res['weights'], color=color, lw=3, label=label)
        if memory == 1 and model == 'semicomplete' and sval in sweep_values:
            fr = mfp.compute_affinity_front(res, on_slope=1.0)
            ax_DG.plot(fr['ton'], fr['DG_on'], color=color, ls=':', lw=3)
            off = fr['DG_off'] > res['params'].DG_min
            ax_DG.plot(fr['toff'][off], fr['DG_off'][off], color=color, ls=':', lw=3)
            expanded, produced = mf.memory_seed_from_primary(res, p=res['params'],
                                                        n_mem=int(res['params'].n_mem))
            ax_n02.plot(res['DG'], expanded, color=color, lw=2, label=label)

    out = mfp._outdir(cfg, 'shared')

    if 'NA_shared' in figs:
        mfp.style_log_axis(ax_NA, T, ylim=(1e0, 1e11), xlim=(0.5, 7.0), hide_xticklabels=False)
        fig_NA.savefig(os.path.join(out, 'N_A_shared.pdf'), dpi=200)

    if 'Z_shared' in figs:
        # ax_Z.axhline(cfg['base']['Z_c'], lw=1, ls='--', color='k')
        if sweep_name == 'h0':
            mfp.style_log_axis(ax_Z, T, ylim=(8e-1, 6e3), xlim=(2.45, 4.5), hide_xticklabels=False)
            model0 = cfg['models'][0]
            # handles = [(Line2D([], [], color=sweep_cols_primary[i], lw=4),
                        # Line2D([], [], color=sweep_cols_memory[i], lw=4))]
            handles = [(Line2D([], [], color=sweep_cols_primary[i], lw=4))
                       for i in range(len(sweep_values))]
            labels = [label_pi_star(results[(model0, sv, 0)]) for sv in sweep_values]
            ax_Z.legend(handles, labels, handler_map={tuple: HandlerTuple(ndivide=2)},
                        handlelength=3, fontsize=22, title=r'$\pi_c$', title_fontsize=24)
        else:
            mfp.style_log_axis(ax_Z, T, ylim=(5e-1, 1e5), xlim=(2.0, 7.0), hide_xticklabels=False)


        fig_Z.savefig(os.path.join(out, 'Z_shared.pdf'), dpi=200, transparent=True)

    if 'pb_shared' in figs:
        ax_pb.set_xlim(0.5, 7)
        ax_pb.tick_params(axis='both', labelsize=30)
        fig_pb.savefig(os.path.join(out, 'pb_shared.pdf'), dpi=200, transparent=True)

    if 'fig4d' in figs:
        for name, marker in [('primary', 'o'), ('recall', 's')]:
            pts = fig4d_pts[name]
            if not pts:
                continue
            P  = np.array([p[1] for p in pts])
            C  = [p[2] for p in pts]
            Da = np.array([p[3] for p in pts])
            # idx = [int(np.argmin(P))] if name == 'primary' else list(np.argsort(Da))
            idx = list(np.argsort(Da))
            # ax_4d.plot(Da[idx], P[idx], color='0.6', lw=1.5, ls='--', zorder=1)
            ax_4d.scatter(Da[idx], P[idx], facecolors=[C[i] for i in idx],
                          edgecolors='k', marker=marker,
                          s=130, linewidth=1.5, zorder=2, label=r'$\mathrm{%s}$' % name)
        rn0 = results.get(('null', cfg.get('null_sval'), 0))
        if rn0 is not None:
            _, P = mf.compute_potency_P(rn0, t_c_star=t_c_star)
            D_ana = mf.dkl_analytic(FIG4D_MU['null'], DG4D_STAR + mf.DG_mf(rn0),
                                    DG_star=DG4D_STAR + mf.DG_min(rn0), sigma2=DG4D_SIG2)
            if not (np.isnan(P) or np.isnan(D_ana)):
                ax_4d.scatter([D_ana], [P], facecolors=[mfp.my_grey], edgecolors='k',
                              marker='o', s=170, linewidth=1.5, zorder=3, label=r'$\mathrm{null}$')

        rn1 = results.get(('null', cfg.get('null_sval'), 1))
        if rn1 is not None:
            _, P = mf.compute_potency_P(rn1, t_c_star=t_c_star)
            D_ana = mf.dkl_analytic(FIG4D_MU['null']*2, DG4D_STAR + mf.DG_mf(rn1),
                                    DG_star=DG4D_STAR + mf.DG_min(rn1), sigma2=DG4D_SIG2)
            if not (np.isnan(P) or np.isnan(D_ana)):
                ax_4d.scatter([D_ana], [P], facecolors=[mfp.my_grey], edgecolors='k',
                                marker='s', s=170, linewidth=1.5, zorder=3)
                
        ax_4d.set_xlim(4.95, 18)
        ax_4d.set_yscale('log')
        ax_4d.tick_params(axis='both', labelsize=30)
        ax_4d.legend(fontsize=24, loc=3)
        fig_4d.savefig(os.path.join(out, 'potency_specificity_D.pdf'), dpi=200, transparent=True)

    if 'DG_shared' in figs:
        rn = results.get(('null', cfg.get('null_sval'), 0))
        if rn is not None:
            fr = mfp.compute_affinity_front(rn, on_slope=1.0)
            ax_DG.plot(fr['ton'], fr['DG_on'], color=mfp.my_grey, ls='--', lw=3)
            off = fr['DG_off'] > rn['params'].DG_min
            ax_DG.plot(fr['toff'][off], fr['DG_off'][off],
                       color=mfp.my_grey, ls='--', lw=3, label='null')
        rn = results.get(('null', cfg.get('null_sval'), 1))
        if rn is not None:
            fr = mfp.compute_affinity_front(rn, on_slope=1.0)
            ax_DG.plot(fr['ton'], fr['DG_on'], color=mfp.my_grey, ls=':', lw=3)
            off = fr['DG_off'] > rn['params'].DG_min
            ax_DG.plot(fr['toff'][off], fr['DG_off'][off],
                        color=mfp.my_grey, ls=':', lw=3)
        ax_DG.set_ylim(bottom=cfg['base']['DG_min'] - 0.1)
        ax_DG.set_xlim(0, T)
        ax_DG.legend(fontsize=14)
        fig_DG.savefig(os.path.join(out, 'DG_shared.pdf'), dpi=150)

    if 'n0_shared' in figs:
        rn = results.get(('null', cfg.get('null_sval'), 0))
        if rn is not None:
            try:
                expanded, produced = mf.memory_seed_from_primary(
                    rn, p=rn['params'], n_mem=int(rn['params'].n_mem))
                ax_n0.plot(rn['DG'], expanded, color=mfp.my_grey, lw=2, label='null')
                ax_n0.plot(rn['DG'], produced * rn['weights'], color=mfp.my_grey, lw=2)
            except ValueError:
                pass                      # null produced no activated cells to seed
        Emax = 12
        E = np.linspace(0, Emax, 10000)
        ax_n0.plot(E, 2e9*np.exp(-(E-16)**2/(2*6)), color='grey', alpha = 0.5, lw = 8)
        ax_n0.set_yscale('log')
        ax_n0.set_ylim(1, 2.5e9)
        ax_n0.legend(fontsize=16, title=r'$\pi_c$')
        fig_n0.savefig(os.path.join(out, 'n0_shared.pdf'), dpi=200, transparent=True)

        rn = results.get(('null', cfg.get('null_sval'), 1))
        if rn is not None:
            try:
                expanded, produced = mf.memory_seed_from_primary(
                    rn, p=rn['params'], n_mem=int(rn['params'].n_mem))
                ax_n02.plot(rn['DG'], expanded, color=mfp.my_grey, lw=2, label='null')
            except ValueError:
                pass                      # null produced no activated cells to seed
        ax_n02.plot(E, 2e9*np.exp(-(E-16)**2/(2*6)), color='grey', alpha = 0.5, lw = 8)
        ax_n02.set_yscale('log')
        ax_n02.set_ylim(1, 2.5e9)
        ax_n02.legend(fontsize=16, title=r'$\pi_c$')
        fig_n02.savefig(os.path.join(out, 'n0_shared2.pdf'), dpi=200, transparent=True)

    if 'P_pi' in figs:
        for name, marker in [('primary', 'o'), ('recall', 's')]:
            pts = fig4d_pts[name]
            if not pts:
                continue
            P  = np.array([p[1] for p in pts])
            C  = [p[2] for p in pts]
            Da = np.array([p[3] for p in pts])
            pi_cs = np.array([p[4] for p in pts])
            # idx = [int(np.argmin(P))] if name == 'primary' else list(np.argsort(Da))
            idx = list(np.argsort(Da))
            ax_P_pi.plot(pi_cs[idx], -P[idx], color='0.6', lw=1.5, ls='--', zorder=1)
            ax_P_pi.scatter(pi_cs[idx], -P[idx], facecolors=[C[i] for i in idx],
                            edgecolors='k', marker=marker,
                            s=130, linewidth=1.5, zorder=2, label=r'$\mathrm{%s}$' % name)
        ax_P_pi.tick_params(axis='both', labelsize=30)
        ax_P_pi.set_xscale('log')
        # ax_P_pi.legend(fontsize=16, title=r'$\pi_c$')
        fig_P_pi.savefig(os.path.join(out, 'P_pi_shared.pdf'), dpi=200, transparent=True)

    if not cfg['show']:
        for f in all_figs:
            plt.close(f)


def make_per_run_figures(cfg, results):
    """Standard per-run panels (antigen, pi, N_B, potency)."""
    T = cfg['T']
    sweep_name, sweep_values = cfg['sweep']
    # colour scale across the sweep (e.g. the green 'summer' scale for an h0 sweep)
    # sweep_cols_primary = mfp.sweep_colors(len(sweep_values), 'autumn')   # memory 0
    sweep_cols_primary =  mfp.colors_sweep   # memory 0
    # sweep_cols_memory  = mfp.sweep_colors(len(sweep_values), 'winter')   # memory 1
    sweep_cols_memory  = mfp.colors_sweep   # memory 1

    for (model, sval, memory), res in results.items():
        
        if len(sweep_values) > 1 and sval in sweep_values:
            pi_c = (_b0 / sval) ** (1 / res['params'].hill)  # for debugging: check that the label matches the run's pi_c
            if model == 'null':
                color, ls = mfp.colors_mem[memory], mfp.styles_sim[0]
            else:
                i_s = sweep_values.index(sval)
                color = sweep_cols_primary[i_s] if memory == 0 else sweep_cols_memory[i_s]
                ls = mfp.styles_mem[memory]
        else:
            pi_c = (_b0 / res['params'].h0) ** (1 / res['params'].hill)
            if model == 'null':
                color, ls = mfp.colors_mem[memory], mfp.styles_sim[0]
            else:
                color, ls = mfp.colors_mem[memory],  mfp.styles_sim[0]

        label = label_pi_star(res)
        out = mfp._outdir(cfg, model, _sweep_tag(cfg, sval), f'memory_{memory}')

        fig_NA, ax_NA = mfp.new_fig(wide=True)
        mfp.plot_NA(ax_NA, res, label=label)
        mfp.style_log_axis(ax_NA, T-1, ylim=(1e0, 1e11))
        fig_NA.savefig(os.path.join(out, 'N_A.pdf'), dpi=150)

        if 'pi' in res:
            fig_pi, ax_pi = mfp.new_fig(wide=True)
            mfp.plot_pi(ax_pi, res, color=color, label=label)
            ax_pi.axhline(pi_star_of(res['params']), lw=1, ls='--',
                          color='grey', alpha=0.8)
            mfp.style_log_axis(ax_pi, T-1, ylim=(5e-1, 1e4))
            fig_pi.savefig(os.path.join(out, 'pi.pdf'), dpi=150)

        fig_NB, ax_NB = mfp.new_fig(wide=True)
        mfp.plot_NB(ax_NB, res, color=color, label=label)
        ax_NB.axhline(1.0, color='k', ls='--', alpha=0.5)
        mfp.style_log_axis(ax_NB, T-1, ylim=(5e-1, 1e6))
        fig_NB.savefig(os.path.join(out, 'N_B.pdf'), dpi=150)

        fig_Z, ax_Z = mfp.new_fig(wide=True)
        mfp.plot_potency(ax_Z, res, color=color, label=label)
        ax_Z.axhline(res['params'].Z_c, lw=1, ls='--', color='k')
        mfp.style_log_axis(ax_Z, T-1, ylim=(5e-3, 1e7), hide_xticklabels=False)
        # ax_Z.set_xlabel('Time', fontsize=16)
        fig_Z.savefig(os.path.join(out, 'Z.pdf'), dpi=150)

        # DG: analytic on/off affinity front (skipped if the run never activates)
        try:
            fr = mfp.compute_affinity_front(res, on_slope=0.78)
            p = res['params']
            fig_DG, ax_DG = mfp.new_fig(wide=True)
            ax_DG.plot(fr['ton'], fr['DG_on'], lw=3, color=color, ls='dotted')
            off = fr['DG_off'] > p.DG_min
            ax_DG.plot(fr['toff'][off], fr['DG_off'][off], lw=3, color=color, ls='dashed')
            offn = fr['DG_off_null'] > p.DG_min
            ax_DG.plot(fr['toff'][offn], fr['DG_off_null'][offn], lw=3,
                       color=mfp.my_grey, ls='dashed')
            ax_DG.axhline(p.DG_min, lw=1, ls='--', color='grey', alpha=0.8)
            ax_DG.set_ylim(p.DG_min - 0.1, 6)
            ax_DG.set_xlim(0, T-1)
            ax_DG.tick_params(axis='both', labelsize=30)
            fig_DG.savefig(os.path.join(out, 'DG.pdf'), dpi=150)
        except (IndexError, ValueError):
            print(f"      [warn] no DG front for {model} "
                  f"{_sweep_tag(cfg, sval)} memory_{memory} (never activated)")

        if not cfg['show']:
            plt.close('all')


# ============================================================
# Entry point
# ============================================================
def main(cfg=CONFIG):
    cfg = normalize_config(cfg)
    mfp.use_style(usetex=cfg['usetex'])
    results = run_experiment(cfg)

    figs = cfg['figures']
    if any(f in figs for f in ('NA_shared', 'Z_shared', 'pb_shared',
                               'fig4d', 'DG_shared', 'n0_shared', 'P_pi')):
        make_shared_figures(cfg, results)
    if 'per_run' in figs:
        make_per_run_figures(cfg, results)

    if cfg['show']:
        plt.show()
    print("Done.")
    return results


if __name__ == '__main__':
    main()