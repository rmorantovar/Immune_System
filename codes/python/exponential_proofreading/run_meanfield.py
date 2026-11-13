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

Two-clone mode (`python run_meanfield.py two_clone`, see TWO_CLONE_CONFIG):
    Instead of the whole naive repertoire, the response starts from 2 clones,
        clone 1 : DG1,        initial size N1
        clone 2 : DG1 + gap,  initial size N1 * ratio
    so `gap` is the affinity gap (clone 1 is the stronger binder when gap > 0)
    and `ratio` = N2(0)/N1(0) the initial size ratio. Any of DG1, gap, N1,
    ratio can be swept like a Parameter, e.g. sweep=('gap', [0.5, 1, 2]).
    Figures: 'two_clone' (clone sizes, clone-1 fraction, pi, final fraction).

Boost mode (`python run_meanfield.py boost`, see BOOST_CONFIG), phase 2:
    The naive repertoire and the memory formed in the primary (phase 0, same
    model and sweep value) respond jointly -- sharing the T-cell pool -- to the
    primary's antigen trajectory, replayed with no feedback of the response.
    Memory uses the phase-1 activation threshold (h0_mem = b0 when
    memory_h0_equals_b0). Readout at t = boost_t_read (default T):
        R = memory-derived / naive-derived cells (expanded clones only),
        phi_naive = 1 / (1 + f R), f = boost_f (memory GC-entry bias),
    compared with R_th = kappa nu_mu Nbar_w and with Mesin et al. 2020.
    One boost run per memory capacity in boost_nmem_values.
    Figures: 'boost' (saved under outroot/.../boost/).
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../library'))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

import mf_lib as mf
import mf_plotting as mfp
from paths import DRIVE


# ============================================================
# CONFIG  -- edit this block, then run `python run_meanfield.py`
# ============================================================
_b0 = 2.0
# Fig 4D analytic-D parameters (naive N(0, sigma2); edge DG_star = -beta* sigma2)
DG4D_STAR = -16.0        # high-affinity edge (negative)
DG4D_SIG2 = 7.0          # naive variance sigma^2

BASE = dict(
    N_A0=1.0, lambda_innate=2.2, threshold_innate=1e4,
    lambda_A=6.2, delta_A=3.5, eta=1.0,
    k_on=1e2 * 2e5 * 1e6 * 24 * 3600 / mf.N_Avg, delta_pi=0.1,
    hill=2.0, beta_star=2.3, K_T=1e4,
    delta_T=0.00, Tcell_growth_factor=2.0,
    tau_eng=0.1, b0=_b0, delta_B=0.00, h0=_b0 / 1000.,
    DG_min=0.0, DG_max=9.0, M=32,
    omega_0=1.0, T_lim=True, N_T0=1e6,
    Z_c=1.e3, n_mem=4e4,
)

_alpha_o = BASE['b0'] / BASE['lambda_A'] + BASE['b0'] / BASE['delta_A']   # null EP exponent
FIG4D_MU  = {'primary': -8.0,    # mu = -zeta*|DG*|,  zeta ~ 0.5  (tilted once)
             'recall':  -16.0,   # mu ~ DG*,          zeta ~ 1    (tilted twice)
             'null': -5.16
            # 'null':    -(_alpha_o / BASE['beta_star']) * abs(DG4D_STAR),          # ~ -5.16
            }

CONFIG = dict(
    # --- what to run ---
    models=['semicomplete'],            # subset of mf.MODELS
    memory_phases=[0, 1],               # 0 = primary, 1 = memory (seeded from 0)
    sweep=('h0', np.flip(np.logspace(np.log10(_b0/1e8), np.log10(_b0/1e1), 16))),     # None -> no sweep (use BASE params as-is).
    # sweep=('h0', np.flip([_b0/1e7, _b0/1e5, _b0/1e3])),     # None -> no sweep (use BASE params as-is).                      
    # sweep = None, # None -> no sweep (use BASE params as-is).
                                        # To sweep: ('h0', [v1, v2, ...]) or any
                                        # Parameter name with a list/array of values.
    null_sval= 1.,      # single value: used to RUN null and to PLOT the null point                               
    sweep_cmap='autumn',                # colormap for multi-value sweeps
                                        #   ('summer' = your green scale; None = discrete)
    T=12.0,
    mode='grid',                        # 'grid' or 'stochastic'
    base=BASE,

    # --- behaviour toggles ---
    memory_h0_equals_b0=True,           # memory phase forces h0 = b0 (original behaviour)
    two_clone=None,                     # None -> whole naive repertoire; see TWO_CLONE_CONFIG

    # --- output / figures ---
    # figures=['NA_shared', 'Z_shared', 'pb_shared', 'DG_shared', 'n0_shared', 'P_shared', 'P_pi', 'per_run'],
    figures=['P_shared'],
    #   available: 'NA_shared', 'Z_shared', 'pb_shared', 'per_run', 'h0_extras'
    outroot=str(DRIVE) + '/documents/Science/research_projects/Immune_System/figures/exponential_proofreading/meanfield/figures',                  # base directory for saved PDFs
    usetex=True,                       # set True if you have a LaTeX toolchain
    show=False,                         # plt.show() at the end instead of/along saving
)

# Two-clone competition: run with `python run_meanfield.py two_clone`.
TWO_CLONE_CONFIG = dict(
    CONFIG,
    two_clone=dict(DG1=0.0,             # DG of clone 1
                   gap=1.0,             # DG2 - DG1 (affinity gap)
                   N1=1.0,              # initial size of clone 1
                   ratio=1.0),          # N2(0) / N1(0)
    sweep=('gap', [0.0, 0.5, 1.0, 2.0, 4.0]),   # or ('ratio', [...]), ('h0', [...]), ...
    # sweep=('ratio', np.logspace(-2, 2, 5)),
    memory_phases=[0, 1],
    figures=['two_clone'],
)

# Boost (phase 2): run with `python run_meanfield.py boost`.
BOOST_CONFIG = dict(
    CONFIG,
    memory_phases=[0, 2],
    boost_nmem_values=np.logspace(2, 6, 9),   # memory capacity n_mem, one boost run each
    boost_f=1.0,                        # memory GC-entry bias, applied post hoc in phi_naive
    boost_t_read=6.0,                   # readout time after the boost (None -> T); Mesin data: day 6
    figures=['boost'],
)
# Mesin et al. 2020: memory-derived fraction of the boosted GC
BOOST_DATA = dict(flu_HA_d6=0.044, lo=0.038, hi=0.11)


# ============================================================
# Parameter construction (fixes the in-place-mutation + p-before-assignment bugs)
# ============================================================
# Defaults filled in for any key omitted from a user config.
CONFIG_DEFAULTS = dict(
    models=['semicomplete'], memory_phases=[0], sweep=None, sweep_cmap='summer',
    T=10.0, mode='grid', base=BASE, memory_h0_equals_b0=True, two_clone=None,
    figures=['NA_shared', 'Z_shared', 'pb_shared', 'per_run'],
    outroot='figures', usetex=True, show=False,
    boost_nmem_values=None, boost_f=1.0, boost_t_read=None, boost_checks=True,
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
    tc = cfg.get('two_clone')
    if tc is not None:
        tc = dict(tc)
        if sweep_name in tc:                # sweeping a two-clone knob
            tc[sweep_name] = sweep_value
        d['clone_DG'] = (tc['DG1'], tc['DG1'] + tc['gap'])
        d['clone_N0'] = (tc['N1'], tc['N1'] * tc['ratio'])
    if sweep_name is not None and not (tc is not None and sweep_name in tc):
        d[sweep_name] = sweep_value
    d['memory'] = memory
    if memory == 1 and cfg.get('memory_h0_equals_b0', True):
        d['h0'] = d['b0']
    if memory == 2:                         # boost: naive h0 as swept, memory h0 as in phase 1
        d['memory'] = 0                     # antigen law without response feedback
        d['h0_mem'] = d['b0'] if cfg.get('memory_h0_equals_b0', True) else d['h0']
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
    sim_mode = 'clones' if cfg.get('two_clone') is not None else cfg['mode']
    results = {}
    memory_seeds = {}      # (model, sweep_value) -> formed_memory from the primary
    primary_params = {}    # (model, sweep_value) -> Parameters of the naive run

    print("Running experiment")
    for model in cfg['models']:
        print(f"    model = {model}")
        if model == 'null' and cfg.get('null_sval') is not None:
            model_vals = [cfg['null_sval']]                       # run null only here
        else:
            model_vals = cfg.get('model_sweep_values', {}).get(model, sweep_values)
        for sval in model_vals:
            if sweep_name is None:
                print("       (fixed parameters, no sweep)")
            else:
                print(f"        {sweep_name} = {sval:.4g}")
            for memory in cfg['memory_phases']:
                if memory == 2:
                    res0 = results.get((model, sval, 0))
                    if res0 is None or model != 'semicomplete' or sim_mode != 'grid':
                        print("      [skip] boost needs a phase-0 semicomplete run on the grid "
                              "(add 0 to memory_phases)")
                        continue
                    results[(model, sval, 2)] = run_boost(cfg, sval, res0)
                    continue
                p = build_params(cfg, sval, memory)
                print(f"            memory = {p.memory:.0g}")
                print(f"            pi_c = {(_b0 / p.h0) ** (1 / p.hill):.4g}")
                seed = memory_seeds.get((model, sval)) if memory == 1 else None
                if memory == 1 and seed is None:
                    print("      [skip] memory phase requested but no primary "
                          "(add 0 to memory_phases)")
                    continue

                print('             starting...', end='', flush=True)
                res = mf.run_simulation(model, p=p, t_span=(0.0, cfg['T']),
                                        mode=sim_mode, memory_seed=seed)
                print('done')
                # Attach the naive run's parameters so memory plots can recover
                # the primary pi_c (which differs from the memory run's own p).
                res['params_primary'] = primary_params.get((model, sval), p)
                results[(model, sval, memory)] = res

                if memory == 0:
                    primary_params[(model, sval)] = p
                    try:
                        # activated = at least doubled (naive grid: N0 = 1 -> 2)
                        _, formed = mf.memory_seed_from_primary(
                            res, p=p, n_mem=int(p.n_mem),
                            threshold=2.0 * mf._N_B(res)[:, 0])
                        memory_seeds[(model, sval)] = formed
                    except ValueError as e:
                        print(f"      [warn] no memory seed for {model}: {e}")
    return results


def run_boost(cfg, sval, res0):
    """
    Boost (phase 2) for one sweep value, one joint run per memory capacity n_mem.

    Naive repertoire + memory seed (expected draw from the primary res0) in one
    'boost' simulation, with N_A replayed from res0 (no feedback). Clones count
    as responding when at least doubled (naive: from 1; memory: from its seed).

    Returns a dict with the primary-derived constants (kappa, Nbar_w, Y, ...) and
    runs[n_mem] = N_naive_tot, N_mem_tot, R, R(t), N_naive(DG), N_mem(DG), r(DG),
    nu_mu, R_th, phi_naive.
    """
    T = cfg['T']
    p0, p = res0['params'], build_params(cfg, sval, 2)
    drive = dict(t=res0['t'], N_A=res0['N_A'])
    w, DG = res0['weights'], res0['DG']
    thr0 = 2.0 * mf._N_B(res0)[:, 0]

    if cfg['boost_checks']:
        rr = mf.run_simulation('semicomplete', p=p0, t_span=(0.0, T), mode='grid', drive=drive)
        ref = mf._N_B(res0)[:, -1] * w
        err = np.max(np.abs(mf._N_B(rr)[:, -1] * w - ref)) / ref.max()
        print(f"            check: phase 0 replayed on its own N_A, max rel. diff {err:.1e}")

    # primary-derived constants for R_th = kappa nu_mu Nbar_w
    N = mf._N_B(res0)[:, -1] * (mf._N_B(res0)[:, -1] > thr0)     # per-clone primary size, expanded clones
    Nbar_w = np.sum(w * N ** 2) / np.sum(w * N)                   # cell-weighted mean clone size
    pi_c, pi_c_mem = pi_star_of(p0), (p.b0 / p.h0_mem) ** (1 / p.hill)
    alpha = p.b0 * (1.0 / p.lambda_A + 1.0 / min(p.delta_A, p.b0))
    kappa = (pi_c / pi_c_mem) ** alpha

    t_read = T if cfg['boost_t_read'] is None else cfg['boost_t_read']
    out = dict(params=p, params_primary=p0, pi_c=pi_c, pi_c_mem=pi_c_mem, alpha=alpha,
               kappa=kappa, Nbar_w=Nbar_w, DG=DG, weights=w, t_read=t_read, runs={})
    n_mems = cfg['boost_nmem_values']
    for n_mem in ([p.n_mem] if n_mems is None else n_mems):
        expanded, formed = mf.memory_seed_expected(res0, n_mem=n_mem, threshold=thr0)
        rb = mf.run_simulation('boost', p=p, t_span=(0.0, T), mode='grid',
                               memory_seed=formed, drive=drive)
        i = int(np.argmin(np.abs(rb['t'] - t_read)))
        N_n = rb['N_Bo'] + rb['N_Ba']                             # per clone, (M, N_t)
        N_m = rb['N_Bo_m'] + rb['N_Ba_m']
        cells_n = w[:, None] * N_n * (N_n > 2.0 * rb['N_Bo'][:, :1])
        cells_m = w[:, None] * N_m * ((N_m > 2.0 * formed[:, None]) & (formed[:, None] > 0))
        R_t = cells_m.sum(0) / np.maximum(cells_n.sum(0), 1e-300)
        with np.errstate(divide='ignore', invalid='ignore'):
            r = np.where(formed > 0, (N_m[:, i] / formed) / N_n[:, i], np.nan)   # memory / naive fold expansion

        nu_mu = np.sum(formed * w) / np.sum(expanded)
        R = R_t[i]
        R_th = kappa * nu_mu * Nbar_w
        phi = 1.0 / (1.0 + cfg['boost_f'] * R)
        out['runs'][n_mem] = dict(
            N_naive_tot=cells_n[:, i].sum(), N_mem_tot=cells_m[:, i].sum(), R=R, t=rb['t'], R_t=R_t,
            N_naive=cells_n[:, i], N_mem=cells_m[:, i], r=r, Y=np.sum(expanded),
            nu_mu=nu_mu, R_th=R_th, phi_naive=phi)
        below = np.nanmin(r / kappa) if np.any(np.isfinite(r)) else np.nan
        print(f"            pi_c={pi_c:9.3g}  n_mem={n_mem:8.3g}  R_sim={R:9.3g}  R_th={R_th:9.3g}  "
              f"R_sim/R_th={R / R_th:7.3g}  phi_naive={phi:6.3f}  min r(DG)/kappa={below:6.3g}")
    out['Y'] = out['runs'][n_mem]['Y']
    return out


# ============================================================
# Plotting
# ============================================================

def _sweep_tag(cfg, sval):
    """Folder-name fragment for one sweep value ('fixed' when not sweeping)."""
    name = cfg['sweep'][0]
    return 'fixed' if name is None else f'{name}_{sval:.4g}'


def fig4d_theory_curve(cfg, results, pts, n_grid=200):
    """Analytic recall potency vs specificity across pi_c (SI, optimal threshold).

    Recall control time, up to a constant:
        lambda_A t_c = (1/r) [ -alpha E_max + log(1 + Y/M) ]
    (memory threshold fixed, so no timing term; memory = M/(Y+M) of the yield).
    P = -lambda_A t_c + const; the constant is fitted to the simulated recall points.
    Returns (exp(D), P) on a fine pi_c grid.
    """
    b = cfg['base']
    lamA, lamB, dA, beta = b['lambda_A'], b['b0'], b['delta_A'], b['beta_star']
    M = b['n_mem']                                  # memory capacity (gamma = 1)
    v_off = min(dA, lamB)
    alpha = lamB * (1.0 / lamA + 1.0 / v_off)       # EP exponent
    r = lamB / lamA

    # primary runs: pi_c, reach E_max, and total yield Y
    pic_sim, Emax_sim, Y_sim, DGmin = [], [], [], 0.0
    for (model, sval, memory), res in results.items():
        if model != 'semicomplete' or memory != 0:
            continue
        p = res['params']
        DGmin = mf.DG_min(res)
        pic_sim.append(pi_star_of(p))
        Emax_sim.append(mf.DG_mf(res) - DGmin)
        expanded, _ = mf.memory_seed_from_primary(res, p=p, n_mem=int(p.n_mem))
        Y_sim.append(np.sum(expanded))              # CHECK: total primary cells
    pic_sim, Emax_sim, Y_sim = map(np.asarray, (pic_sim, Emax_sim, Y_sim))

    # reach: slope -1 in log pi_c (theory), intercept from the runs
    E0 = np.mean(Emax_sim + np.log(pic_sim))
    Emax_of = lambda pic: E0 - np.log(pic)

    # yield: exact integral of the exponential tail; prefactor from the runs
    logY_shape = lambda Em: alpha * Em + np.log(np.expm1((beta - alpha) * Em) / (beta - alpha))
    logY0 = np.mean(np.log(Y_sim) - logY_shape(Emax_of(pic_sim)))

    def P_shape(pic):
        Em = Emax_of(pic)
        Y = np.exp(logY0 + logY_shape(Em))
        return (alpha * Em - np.log1p(Y / M)) / r   # = -lambda_A t_c + const

    pic_pts = np.array([pt[4] for pt in pts])
    P_pts = np.array([pt[1] for pt in pts])
    const = np.mean(P_pts - P_shape(pic_pts))

    pic = np.logspace(np.log10(pic_sim.min()), np.log10(pic_sim.max()), n_grid)
    Em = Emax_of(pic)
    D = np.array([mf.dkl_analytic(FIG4D_MU['recall'], DG4D_STAR + e + DGmin,
                                  DG_star=DG4D_STAR + DGmin, sigma2=DG4D_SIG2)
                  for e in Em])
    return np.exp(D), P_shape(pic) + const


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

    sweep_cols_primary = mfp.sweep_colors(len(sweep_values), 'turbo')   # memory 0
    # sweep_cols_primary =  mfp.colors_sweep   # memory 0
    sweep_cols_memory  = mfp.sweep_colors(len(sweep_values), 'turbo_r')   # memory 1
    # sweep_cols_memory  = mfp.colors_sweep   # memory 1

    P_shared_pts = {'primary': [], 'recall': []}
    P_shared_map = {('semicomplete', 0): 'primary', ('semicomplete', 1): 'recall'}

    opt_sval = None
    t_c_star = 0
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
        if memory == 2:                   # boost summaries: see make_boost_figures
            continue
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
                    marker='s', mark_Zc=True, mark_t = t_c_star)
                if opt_sval is None or sval == opt_sval:
                    ax_Z.axvline(t_c_star, lw=2, ls='--', color='k')
                    ax_Z.axhline(res['params'].Z_c, lw=2, ls='--', color='k')

        mfp.plot_pb_from_potency(ax_pb, res, color=color, ls=ls, label=label)

        if (model, memory) in P_shared_map:
            name = P_shared_map[(model, memory)]
            D, (_, P)  = mf.compute_specificity(res), mf.compute_potency_P(res, t_c_star=t_c_star)
            D_ana = mf.dkl_analytic(FIG4D_MU[name], DG4D_STAR + mf.DG_mf(res),
                                    DG_star=DG4D_STAR + mf.DG_min(res), sigma2=DG4D_SIG2)
            if not (np.isnan(D) or np.isnan(P)):
                if memory == 0:
                    if sval == opt_sval:
                        P_shared_pts[name].append((D, P, color, D_ana, pi_c))
                else:
                    P_shared_pts[name].append((D, P, color, D_ana, pi_c))
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
        fig_NA.savefig(os.path.join(out, 'N_A_shared.pdf'), dpi=800)

    if 'Z_shared' in figs:
        # ax_Z.axhline(cfg['base']['Z_c'], lw=1, ls='--', color='k')
        if sweep_name == 'h0':
            mfp.style_log_axis(ax_Z, T, ylim=(2e-1, 6e3), xlim=(2.45, 5.0))
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


        fig_Z.savefig(os.path.join(out, 'Z_shared.pdf'), dpi=800, transparent=True)

    if 'pb_shared' in figs:
        ax_pb.set_xlim(0.5, 7)
        ax_pb.tick_params(axis='both', labelsize=30)
        fig_pb.savefig(os.path.join(out, 'pb_shared.pdf'), dpi=800, transparent=True)

    if 'P_shared' in figs:
        # for name, marker in [('primary', 'o'), ('recall', 's')]:
        for name, marker in [('recall', 's')]:
            pts = P_shared_pts[name]
            if not pts:
                continue
            P  = np.array([pt[1] for pt in pts])
            C  = [pt[2] for pt in pts]
            Da = np.array([pt[3] for pt in pts])
            # idx = [int(np.argmin(P))] if name == 'primary' else list(np.argsort(Da))
            idx = list(np.argsort(Da))
            # ax_4d.plot(Da[idx], P[idx], color='0.6', lw=1.5, ls='--', zorder=1)
            ax_4d.scatter(np.exp(Da[idx]), P[idx], facecolors= [C[i] for i in idx],
                          edgecolors='k', marker=marker,
                          s=130, linewidth=1.5, zorder=2, label=r'$\mathrm{%s}$' % name)
        rn0 = results.get(('null', cfg.get('null_sval'), 0))
        if P_shared_pts['recall']:
            x_th, P_th = fig4d_theory_curve(cfg, results, P_shared_pts['recall'])
            ax_4d.plot(x_th, P_th, color='k', lw=2, ls='--', zorder=1,
                       label=r'$\mathrm{theory}$')
        if rn0 is not None:
            _, P = mf.compute_potency_P(rn0, t_c_star=t_c_star)
            D_ana = mf.dkl_analytic(FIG4D_MU['null'], DG4D_STAR + mf.DG_mf(rn0),
                                    DG_star=DG4D_STAR + mf.DG_min(rn0), sigma2=DG4D_SIG2)
            if not (np.isnan(P) or np.isnan(D_ana)):
                ax_4d.scatter([np.exp(D_ana)], [P], facecolors=[mfp.my_grey], edgecolors='k',
                              marker='o', s=170, linewidth=1.5, zorder=3, label=r'$\mathrm{null}$')

        rn1 = results.get(('null', cfg.get('null_sval'), 1))
        if rn1 is not None:
            _, P = mf.compute_potency_P(rn1, t_c_star=t_c_star)
            D_ana = mf.dkl_analytic(FIG4D_MU['null']*2, DG4D_STAR + mf.DG_mf(rn1),
                                    DG_star=DG4D_STAR + mf.DG_min(rn1), sigma2=DG4D_SIG2)
            if not (np.isnan(P) or np.isnan(D_ana)):
                ax_4d.scatter([np.exp(D_ana)], [P], facecolors=[mfp.my_grey], edgecolors='k',
                                marker='s', s=170, linewidth=1.5, zorder=3)
                
        ax_4d.set_xlim(9e5, 4e7)
        ax_4d.set_ylim(0, 11)
        ax_4d.set_xscale('log')
        # ax_4d.set_yscale('log')
        ax_4d.tick_params(axis='both', labelsize=44)
        ax_4d.legend(fontsize=24, loc=3)
        fig_4d.savefig(os.path.join(out, 'potency'
        '.pdf'), dpi=800, transparent=True)

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
        fig_DG.savefig(os.path.join(out, 'DG_shared.pdf'), dpi=800)

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
        fig_n0.savefig(os.path.join(out, 'n0_shared.pdf'), dpi=800, transparent=True)

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
        fig_n02.savefig(os.path.join(out, 'n0_shared2.pdf'), dpi=800, transparent=True)

    if 'P_pi' in figs:
        for name, marker in [('primary', 'o'), ('recall', 's')]:
            pts = P_shared_pts[name]
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
        fig_P_pi.savefig(os.path.join(out, 'P_pi_shared.pdf'), dpi=800, transparent=True)

    if not cfg['show']:
        for f in all_figs:
            plt.close(f)


def make_two_clone_figures(cfg, results):
    """Two-clone competition: clone sizes, clone-1 fraction, pi, final fraction.

    Colour = sweep value; solid = clone 1 (DG1), dashed = clone 2 (DG1 + gap).
    """
    T = cfg['T']
    sweep_name, sweep_values = cfg['sweep']
    cols = mfp.sweep_colors(len(sweep_values), cfg['sweep_cmap'])
    out = mfp._outdir(cfg, 'two_clone')
    name_tex = (sweep_name or '').replace('_', r'\_')

    def label_of(sval):
        return None if sweep_name is None else r'$\mathrm{%s} = %.3g$' % (name_tex, sval)

    final_frac = {}                       # (model, memory) -> [(sval, f1)]
    for model in cfg['models']:
        for memory in cfg['memory_phases']:
            fig_NB, ax_NB = mfp.new_fig()
            fig_f, ax_f = mfp.new_fig()
            fig_pi, ax_pi = mfp.new_fig()
            for i, sval in enumerate(sweep_values):
                res = results.get((model, sval, memory))
                if res is None:
                    continue
                t, N_B = res['t'], mf._N_B(res)
                c, lab = cols[i], label_of(sval)
                ax_NB.plot(t, N_B[0], lw=5, color=c, ls='-', label=lab)
                ax_NB.plot(t, N_B[1], lw=5, color=c, ls='--')
                f1 = N_B[0] / (N_B[0] + N_B[1])
                ax_f.plot(t, f1, lw=5, color=c, label=lab)
                if 'pi' in res:
                    ax_pi.plot(t, res['pi'][0], lw=5, color=c, ls='-', label=lab)
                    ax_pi.plot(t, res['pi'][1], lw=5, color=c, ls='--')
                    ax_pi.axhline(pi_star_of(res['params']), lw=1, ls=':', color=c)
                final_frac.setdefault((model, memory), []).append((sval, f1[-1]))

            tag = f'{model}_memory_{memory}'
            mfp.style_log_axis(ax_NB, T, ylim=(5e-1, None))
            ax_NB.legend(fontsize=20)
            fig_NB.savefig(os.path.join(out, f'N_B_{tag}.pdf'), dpi=800)

            ax_f.set_xlim(0, T)
            ax_f.set_ylim(-0.02, 1.02)
            ax_f.tick_params(axis='both', labelsize=44)
            ax_f.legend(fontsize=20)
            fig_f.savefig(os.path.join(out, f'frac_clone1_{tag}.pdf'), dpi=800)

            if ax_pi.lines:
                mfp.style_log_axis(ax_pi, T, ylim=(1e-2, None))
                ax_pi.legend(fontsize=20)
                fig_pi.savefig(os.path.join(out, f'pi_{tag}.pdf'), dpi=800)
            if not cfg['show']:
                for f in (fig_NB, fig_f, fig_pi):
                    plt.close(f)

    # final clone-1 fraction vs the swept knob (primary o, recall s)
    if sweep_name is not None and len(sweep_values) > 1:
        fig, ax = mfp.new_fig()
        for (model, memory), pts in final_frac.items():
            x, y = np.array(pts).T
            ax.plot(x, y, lw=3, ls=mfp.styles_mem[memory], color=mfp.colors_mem[memory],
                    marker='os'[memory], ms=14, mec='k',
                    label=r'$\mathrm{%s}$' % ('primary', 'recall')[memory])
            if np.all(x > 0) and x.max() / x.min() > 100:
                ax.set_xscale('log')
        ax.set_xlabel(r'$\mathrm{%s}$' % name_tex, fontsize=30)
        ax.set_ylim(-0.02, 1.02)
        ax.tick_params(axis='both', labelsize=30)
        ax.legend(fontsize=20)
        fig.savefig(os.path.join(out, 'final_frac_clone1.pdf'), dpi=800)
        if not cfg['show']:
            plt.close(fig)


def make_per_run_figures(cfg, results):
    """Standard per-run panels (antigen, pi, N_B, potency)."""
    T = cfg['T']
    sweep_name, sweep_values = cfg['sweep']
    # colour scale across the sweep (e.g. the green 'summer' scale for an h0 sweep)
    sweep_cols_primary = mfp.sweep_colors(len(sweep_values), 'autumn')   # memory 0
    # sweep_cols_primary =  mfp.colors_sweep   # memory 0
    sweep_cols_memory  = mfp.sweep_colors(len(sweep_values), 'winter')   # memory 1
    # sweep_cols_memory  = mfp.colors_sweep   # memory 1

    for (model, sval, memory), res in results.items():
        if memory == 2:                   # boost summaries: see make_boost_figures
            continue
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

        fig_NA, ax_NA = mfp.new_fig(wide=False)
        mfp.plot_NA(ax_NA, res, label=label)
        mfp.style_log_axis(ax_NA, T-1, ylim=(1e0, 1e11))
        ax_NA.tick_params(axis='y', which='both', labelsize=44, labelleft=True if memory == 0 else False)
        ax_NA.tick_params(axis='x', labelsize=44, labelbottom= False)
        fig_NA.savefig(os.path.join(out, 'N_A.pdf'), dpi=800)

        if 'pi' in res:
            fig_pi, ax_pi = mfp.new_fig(wide=False)
            mfp.plot_pi(ax_pi, res, color=color, label=label)
            ax_pi.axhline(pi_star_of(res['params']), lw=1, ls='--',
                          color='grey', alpha=0.8)
            mfp.style_log_axis(ax_pi, T-1, ylim=(5e-1, 1e4))
            ax_pi.tick_params(axis='y', which='both', labelsize=44, labelleft=True if memory == 0 else False)
            ax_pi.tick_params(axis='x', labelsize=44, labelbottom= False)
            fig_pi.savefig(os.path.join(out, 'pi.pdf'), dpi=800)

        fig_NB, ax_NB = mfp.new_fig(wide=False)
        mfp.plot_NB(ax_NB, res, color=color, label=label)
        ax_NB.axhline(1.0, color='k', ls='--', alpha=0.5)
        mfp.style_log_axis(ax_NB, T-1, ylim=(5e-1, 1e6))
        ax_NB.tick_params(axis='y', which='both', labelsize=44, labelleft=True if memory == 0 else False)
        ax_NB.tick_params(axis='x', labelsize=44, labelbottom= False)
        fig_NB.savefig(os.path.join(out, 'N_B.pdf'), dpi=800)

        fig_Z, ax_Z = mfp.new_fig(wide=False)
        mfp.plot_potency(ax_Z, res, color=color, label=label)
        ax_Z.axhline(res['params'].Z_c, lw=1, ls='--', color='k')
        mfp.style_log_axis(ax_Z, T-1, ylim=(5e-3, 1e7), hide_xticklabels=False)
        # ax_Z.set_xlabel('Time', fontsize=16)
        fig_Z.savefig(os.path.join(out, 'Z.pdf'), dpi=800)

        # DG: analytic on/off affinity front (skipped if the run never activates)
        try:
            fr = mfp.compute_affinity_front(res, on_slope=0.85)
            p = res['params']
            fig_DG, ax_DG = mfp.new_fig(wide=False)
            ax_DG.plot(fr['ton'], fr['DG_on'], lw=6, color=color, ls='-')
            off = fr['DG_off'] > p.DG_min
            ax_DG.plot(fr['toff'][off], fr['DG_off'][off], lw=6, color=color, ls='dashed')
            offn = fr['DG_off_null'] > p.DG_min
            # ax_DG.plot(fr['toff'][offn], fr['DG_off_null'][offn], lw=6,
            #            color=mfp.my_grey, ls='dashed')
            ax_DG.axhline(p.DG_min, lw=1, ls='--', color='grey', alpha=0.8)
            ax_DG.set_ylim(p.DG_min - 0.1, 7)
            ax_DG.set_xlim(0, T-1)
            ax_DG.tick_params(axis='y', which='both', labelsize=44, labelleft=True if memory == 0 else False)
            ax_DG.tick_params(axis='x', labelsize=44, labelbottom= True)
            fig_DG.savefig(os.path.join(out, 'DG.pdf'), dpi=800)
        except (IndexError, ValueError):
            print(f"      [warn] no DG front for {model} "
                  f"{_sweep_tag(cfg, sval)} memory_{memory} (never activated)")

        if not cfg['show']:
            plt.close('all')


def make_boost_figures(cfg, results):
    """Boost (phase 2): R vs pi_c, phi_naive map, r(DG), R vs n_mem."""
    boosts = sorted(((res['pi_c'], res) for (model, sval, memory), res in results.items()
                     if memory == 2), key=lambda x: x[0])
    if not boosts:
        print("[boost] no phase-2 results to plot")
        return
    out = mfp._outdir(cfg, 'boost', f"t_read_{boosts[0][1]['t_read']:g}")
    pics = np.array([pc for pc, _ in boosts])
    n_mems = sorted(boosts[0][1]['runs'])
    cols = mfp.sweep_colors(len(n_mems), 'winter')
    R_of = lambda m: m / (1.0 - m) / cfg['boost_f']           # memory fraction -> R (inverts phi_naive)
    nmem_label = lambda n: r'$10^{%.1f}$' % np.log10(n)
    # representative runs: lowest, middle and highest pi_c, at the middle memory capacity
    rep = [0, len(boosts) // 2, len(boosts) - 1] if len(boosts) > 2 else list(range(len(boosts)))
    rep_cols = mfp.sweep_colors(len(rep), 'autumn')
    n_rep = n_mems[len(n_mems) // 2]

    # 1. R vs pi_c: simulation (markers) and R_th (lines), one colour per n_mem; data band
    fig, ax = mfp.new_fig()
    for c, n in zip(cols, n_mems):
        ax.plot(pics, [b['runs'][n]['R'] for _, b in boosts], 'o', color=c, mec='k', ms=10,
                label=nmem_label(n))
        ax.plot(pics, [b['runs'][n]['R_th'] for _, b in boosts], '-', color=c, lw=2)
    ax.axhspan(R_of(BOOST_DATA['lo']), R_of(BOOST_DATA['hi']), color='grey', alpha=0.3)
    ax.axhline(R_of(BOOST_DATA['flu_HA_d6']), color='k', ls='--', lw=1.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\pi_c$', fontsize=30)
    ax.set_ylabel(r'$R$', fontsize=30)
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(fontsize=16, title=r'$n_{\mathrm{mem}}$', title_fontsize=18, ncol=2)
    fig.savefig(os.path.join(out, 'R_pi.pdf'), dpi=800, transparent=True)

    # 2. phi_naive over (pi_c, n_mem), contours at the data memory fractions
    phi = np.array([[b['runs'][n]['phi_naive'] for _, b in boosts] for n in n_mems])
    fig, ax = mfp.new_fig()
    mesh = ax.pcolormesh(pics, n_mems, phi, shading='nearest', cmap='viridis', vmin=0, vmax=1)
    fig.colorbar(mesh, ax=ax).set_label(r'$\phi_{\mathrm{naive}}$', fontsize=30)
    if phi.shape[0] > 1 and phi.shape[1] > 1:
        ax.contour(pics, n_mems, 1.0 - phi, levels=[BOOST_DATA['flu_HA_d6']], colors='w', linestyles='-')
        ax.contour(pics, n_mems, 1.0 - phi, levels=[BOOST_DATA['hi']], colors='w', linestyles='--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\pi_c$', fontsize=30)
    ax.set_ylabel(r'$n_{\mathrm{mem}}$', fontsize=30)
    ax.tick_params(axis='both', labelsize=30)
    fig.savefig(os.path.join(out, 'phi_naive_map.pdf'), dpi=800, transparent=True)

    # 3. r(DG) = memory / naive fold expansion, with kappa, for the representative runs
    fig, ax = mfp.new_fig()
    for c, k in zip(rep_cols, rep):
        pc, b = boosts[k]
        ax.plot(b['DG'], b['runs'][n_rep]['r'], 'o-', color=c, mec='k', ms=8, lw=2,
                label=mfp.pi_star_label(pc))
        ax.axhline(b['kappa'], color=c, ls='--', lw=1.5)
    ax.set_yscale('log')
    ax.set_xlabel(r'$\Delta G$', fontsize=30)
    ax.set_ylabel(r'$r(\Delta G)$', fontsize=30)
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(fontsize=18, title=r'$\pi_c$', title_fontsize=20)
    fig.savefig(os.path.join(out, 'r_DG.pdf'), dpi=800, transparent=True)

    # 4. R vs n_mem at fixed pi_c; saturation expected at n_mem ~ Y (dotted)
    fig, ax = mfp.new_fig()
    for c, k in zip(rep_cols, rep):
        pc, b = boosts[k]
        ax.plot(n_mems, [b['runs'][n]['R'] for n in n_mems], 'o', color=c, mec='k', ms=10,
                label=mfp.pi_star_label(pc))
        ax.plot(n_mems, [b['runs'][n]['R_th'] for n in n_mems], '-', color=c, lw=2)
        ax.axvline(b['Y'], color=c, ls=':', lw=2)
    ax.axhspan(R_of(BOOST_DATA['lo']), R_of(BOOST_DATA['hi']), color='grey', alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(min(n_mems) / 2, max(n_mems) * 2)
    ax.set_xlabel(r'$n_{\mathrm{mem}}$', fontsize=30)
    ax.set_ylabel(r'$R$', fontsize=30)
    ax.tick_params(axis='both', labelsize=30)
    ax.legend(fontsize=18, title=r'$\pi_c$', title_fontsize=20)
    fig.savefig(os.path.join(out, 'R_nmem.pdf'), dpi=800, transparent=True)

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
    if cfg['two_clone'] is not None:
        # the repertoire-level figures assume the full grid and an h0 sweep
        skipped = [f for f in figs if f != 'two_clone']
        if skipped:
            print(f"[two_clone] skipping repertoire figures: {skipped}")
        make_two_clone_figures(cfg, results)
        if cfg['show']:
            plt.show()
        print("Done.")
        return results

    if any(f in figs for f in ('NA_shared', 'Z_shared', 'pb_shared',
                               'P_shared', 'DG_shared', 'n0_shared', 'P_pi')):
        make_shared_figures(cfg, results)
    if 'per_run' in figs:
        make_per_run_figures(cfg, results)
    if 'boost' in figs:
        make_boost_figures(cfg, results)

    if cfg['show']:
        plt.show()
    print("Done.")
    return results


if __name__ == '__main__':
    args = sys.argv[1:]
    main(TWO_CLONE_CONFIG if 'two_clone' in args else BOOST_CONFIG if 'boost' in args else CONFIG)