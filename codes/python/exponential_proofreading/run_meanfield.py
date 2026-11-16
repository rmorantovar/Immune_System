"""
run_meanfield.py -- driver for the unified mean-field models (mf_lib.run_response).

Command line (see `python run_meanfield.py --help`; options override the preset, the preset
overrides CONFIG):
    python run_meanfield.py [preset] [--model M] [--repertoire R] [--responses naive [recall]]
                            [--antigen feedback|replay] [--with-naive] [--sweep NAME[=VALUES] ...]
                            [--set PARAM=VALUE ...] [--T T] [--t-obs T1 T2 ...] [--plot] [--dry-run]
    no --sweep (and no preset sweep) -> single run. Examples:
    python run_meanfield.py --model two_stage --with-naive --antigen replay --sweep pi_c --sweep N_T0 --plot
    python run_meanfield.py --model continuous --plot                      (single run: Fig. 2 panels)
    python run_meanfield.py --repertoire two_clone --sweep gap=0,0.5,1,2,4 --plot

Choose in CONFIG (or a PRESET, `python run_meanfield.py <preset>`):
    model       'continuous' : antigen loading, T-cell help and division at the same time
                'two_stage'  : loading for t < t_switch, then expansion limited by a finite
                               pool of N_T0 T cells (the only coupling between clones)
    repertoire  'polyclonal' (DG grid weighted by the density of states) | 'two_clone'
                (two explicit clones, see TWO_CLONE)
    responses   ['naive'] or ['naive', 'recall']; the recall is seeded with the memory
                formed in the naive response of the same parameters
    recall      with_naive : also include the naive repertoire in the recall
                antigen    : 'feedback' (the response controls the antigen through its
                             potency Z) | 'replay' (the naive response's antigen
                             trajectory, no feedback)
                h0_mem     : memory activation constant ('b0' -> h0_mem = b0; None -> h0)
                draw       : memory seed 'expected' (deterministic) | 'random' (multinomial)
                f          : memory GC-entry bias in phi_naive = 1 / (1 + f R)
    sweep       {} single run | {name: values} | {name1: values1, name2: values2} (grid).
                Names: any mf.Parameters field (h0, N_T0, n_mem, t_switch, h_max, ...)
                or a two-clone knob (DG1, gap, N1, ratio).
    t_obs       times at which the end-point observables are evaluated
    dynamics    time courses to keep (see response_dynamics); observables: see OBSERVABLES

Outputs, one folder per run (name from model, repertoire, recall and sweep):
    DATA_ROOT/<run>/observables.csv  one row per (sweep point, response, t_obs)
    DATA_ROOT/<run>/dynamics.pkl     {(point, response): {quantity: array}}
    DATA_ROOT/<run>/config.json
Figures: `python plot_meanfield.py <run>` (written under FIG_ROOT/<run>/).
The previous driver, with the manuscript figures, is run_meanfield_legacy.py.
"""

import os, sys, json, pickle, itertools, argparse
from dataclasses import fields
sys.path.append(os.path.join(os.path.dirname(__file__), '../../library'))

import numpy as np
import pandas as pd

import mf_lib as mf
from paths import DRIVE

DATA_ROOT = f'{DRIVE}/research_files/Immune_system/exponential_proofreading/meanfield/runs'
FIG_ROOT = f'{DRIVE}/documents/Science/research_projects/Immune_System/figures/exponential_proofreading/meanfield/runs'


# ============================================================
# CONFIG
# ============================================================
_b0 = 2.0
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
    t_switch=6.0, h_max=1e3,
)
TWO_CLONE = dict(DG1=0.0, gap=1.0, N1=1.0, ratio=1.0)     # clone 2 at DG1 + gap, size N1 * ratio
H0_SWEEP = np.flip(np.logspace(np.log10(_b0 / 1e8), np.log10(_b0 / 1e1), 16))   # pi_c from 3 to 1e4

CONFIG = dict(
    model='continuous',
    repertoire='polyclonal',
    responses=['naive', 'recall'],
    recall=dict(with_naive=False, antigen='feedback', h0_mem='b0', draw='expected', f=1.0),
    sweep={'h0': H0_SWEEP},
    T=12.0,
    t_obs=[12.0],
    dynamics=['N_A', 'Z', 'Y', 'L_act', 'T_free', 'DG_dist'],
    observables=['Y', 'L_act', 'N1', 'Z', 'D_KL', 'S2', 'alpha_fit', 'DG_front', 'DG_hi', 'DG_lo', 'f_best',
                 'N_A', 'T_free_min', 't_c', 'R', 'phi_naive', 'R_th', 'kappa'],
    base=BASE,
    two_clone=TWO_CLONE,
    mode='grid',            # polyclonal repertoire: 'grid' or 'stochastic' (with seed)
    seed=0,
    name=None,              # run folder name (None -> built from the config)
)

PRESETS = {
    # single primary + recall at the BASE parameters, with per-clone time courses (Fig. 2 panels)
    'single': dict(sweep={}, dynamics=['N_A', 'Z', 'Y', 'L_act', 'T_free', 'N_B_clones', 'pi_clones',
                                       'DG_dist']),
    # naive + memory responding together to the naive antigen trajectory (prime-boost)
    'boost': dict(recall=dict(CONFIG['recall'], with_naive=True, antigen='replay'),
                  sweep={'h0': H0_SWEEP, 'n_mem': np.logspace(2, 6, 9)}, t_obs=[6.0, 12.0],
                  dynamics=['N_A', 'Z', 'Y', 'L_act', 'T_free', 'N_B_clones']),
    # two-stage model: T-cell competition in the boost, pi_c x T-cell pool
    'two_stage': dict(model='two_stage', T=16.0, t_obs=[16.0],
                      recall=dict(CONFIG['recall'], with_naive=True, antigen='replay'),
                      sweep={'h0': H0_SWEEP, 'N_T0': np.logspace(3, 7, 9)}),
    # two-clone competition: affinity gap
    'two_clone': dict(repertoire='two_clone', sweep={'gap': [0.0, 0.5, 1.0, 2.0, 4.0]},
                      dynamics=['N_A', 'Z', 'Y', 'T_free', 'N_B_clones', 'pi_clones']),
}


# ============================================================
# Parameters
# ============================================================
PARAM_FIELDS = {f.name for f in fields(mf.Parameters)}


def full_config(cfg):
    """CONFIG defaults + cfg (base, recall and two_clone merged key by key); checks the choices."""
    out = {**CONFIG, **cfg}
    for key in ('base', 'recall', 'two_clone'):
        out[key] = {**CONFIG[key], **cfg.get(key, {})}
    out['sweep'] = {k: list(np.atleast_1d(v)) for k, v in out['sweep'].items()}
    checks = [(out['model'] in mf.RESPONSE_MODELS, f"model must be one of {list(mf.RESPONSE_MODELS)}"),
              (out['repertoire'] in ('polyclonal', 'two_clone'), "repertoire must be 'polyclonal' or 'two_clone'"),
              (set(out['responses']) <= {'naive', 'recall'} and 'naive' in out['responses'],
               "responses must be ['naive'] or ['naive', 'recall']"),
              (out['recall']['antigen'] in ('feedback', 'replay'), "recall antigen must be 'feedback' or 'replay'"),
              (out['recall']['draw'] in ('expected', 'random'), "recall draw must be 'expected' or 'random'"),
              (len(out['sweep']) <= 2, "sweep over at most 2 parameters")]
    for k in out['sweep']:
        checks.append((k in PARAM_FIELDS or k in out['two_clone'], f"unknown sweep parameter '{k}'"))
        checks.append((k not in out['two_clone'] or out['repertoire'] == 'two_clone',
                       f"'{k}' is a two-clone parameter: needs repertoire='two_clone'"))
    for ok, msg in checks:
        if not ok:
            raise ValueError(msg)
    return out


def build_params(cfg, point):
    """mf.Parameters for one sweep point (dict name -> value). Never mutates cfg."""
    d, tc = dict(cfg['base']), dict(cfg['two_clone'])
    for k, v in point.items():
        (tc if k in tc else d)[k] = v
    if cfg['repertoire'] == 'two_clone':
        d['clone_DG'] = (tc['DG1'], tc['DG1'] + tc['gap'])
        d['clone_N0'] = (tc['N1'], tc['N1'] * tc['ratio'])
    h0_mem = cfg['recall']['h0_mem']
    d['h0_mem'] = d['b0'] if h0_mem == 'b0' else h0_mem
    return mf.Parameters(**d)


def pi_c_of(p, h0=None):
    return (p.b0 / (p.h0 if h0 is None else h0)) ** (1 / p.hill)


def run_name(cfg):
    if cfg['name']:
        return cfg['name']
    name = f"{cfg['model']}_{cfg['repertoire']}"
    if 'recall' in cfg['responses']:
        rec = cfg['recall']
        name += f"_recall{'+naive' if rec['with_naive'] else ''}_{rec['antigen']}"
    name += '_single' if not cfg['sweep'] else '_sweep_' + '_'.join(cfg['sweep'])
    changed = [f'{k}={v:.3g}' for k, v in {**cfg['base'], **cfg['two_clone']}.items()
               if k not in cfg['sweep'] and v != {**BASE, **TWO_CLONE}.get(k)]
    return '_'.join([name] + changed)                 # changed base parameters keep runs apart


# ============================================================
# Observables  (v = mf.pop_view of one population; i = time index)
# ============================================================
def _expanded(v, i):
    """Clones that at least doubled since t = 0."""
    return v['N_B'][:, i] > 2.0 * v['N_B'][:, 0]


def _ever_active(v, i):
    return np.any(v['pi'][:, :i + 1] > pi_c_of(v['params']), axis=1)


def obs_D_KL(v, i):
    """Specificity D_KL(Omega || Omega_0) of the responding cells (as mf.compute_specificity)."""
    w = v['weights']
    Omega = w * v['N_B'][:, i] * _ever_active(v, i)
    if Omega.sum() <= 0:
        return np.nan
    q, r = w / w.sum(), Omega / Omega.sum()
    m = r > 0
    return float(np.sum(r[m] * np.log(r[m] / q[m])))


def obs_S2(v, i):
    """Renyi-2 clone entropy -ln sum_clones p^2 of the expanded clones (weights = clones per bin);
    as -ln(lambda2) in zipf.py, which small samples estimate without bias."""
    act = _expanded(v, i)
    w, n = v['weights'][act], v['N_B'][act, i]
    N = np.sum(w * n)
    return float(-np.log(np.sum(w * (n / N) ** 2))) if N > 0 else np.nan


def obs_alpha_fit(v, i):
    """-slope of ln N_B vs DG over the expanded clones (EP exponent)."""
    act = _expanded(v, i)
    if act.sum() < 3:
        return np.nan
    return float(-np.polyfit(v['DG'][act], np.log(v['N_B'][act, i]), 1)[0])


def _range_above_2(v, i, f):
    """f (min or max) of the DG of clones with more than 2 cells (legacy mf.DG_mf / mf.DG_min)."""
    act = v['N_B'][:, i] > 2.0
    return float(f(v['DG'][act])) if act.any() else np.nan


def _max_over_expanded(v, i, x):
    act = _expanded(v, i)
    return float(x[act].max()) if act.any() else np.nan


# per-population observables: name -> f(v, i)
POP_OBSERVABLES = {
    'Y': lambda v, i: float(np.sum(v['weights'] * v['N_B'][:, i] * _expanded(v, i))),
    'L_act': lambda v, i: float(np.sum(v['weights'] * _expanded(v, i))),
    'N1': lambda v, i: _max_over_expanded(v, i, v['N_B'][:, i]),
    'Z': lambda v, i: float(mf.compute_potency_t(v)[i]),
    'D_KL': obs_D_KL,
    'S2': obs_S2,
    'alpha_fit': obs_alpha_fit,
    'DG_front': lambda v, i: _max_over_expanded(v, i, v['DG']),
    'DG_hi': lambda v, i: _range_above_2(v, i, np.max),
    'DG_lo': lambda v, i: _range_above_2(v, i, np.min),
    'f_best': lambda v, i: float(v['weights'][0] * v['N_B'][0, i] / np.sum(v['weights'] * v['N_B'][:, i])),
}
POP_ONLY = ('Z', 'D_KL', 'S2', 'alpha_fit', 'DG_hi', 'DG_lo', 'f_best')   # not defined for the naive + memory sum
OBSERVABLES = list(POP_OBSERVABLES) + ['N_A', 'T_free_min', 't_c', 'R', 'phi_naive', 'R_th', 'kappa']


def response_observables(cfg, res, i, extra):
    """All requested observables of one response at time index i."""
    want, out = cfg['observables'], {}
    pops = res['pops'] + (['all'] if len(res['pops']) > 1 else [])
    views = {pop: mf.pop_view(res, pop) for pop in pops}
    for name, f in POP_OBSERVABLES.items():
        if name in want:
            for pop, v in views.items():
                if not (pop == 'all' and name in POP_ONLY):
                    out[f'{name}_{pop}'] = f(v, i)
    if 'N_A' in want:
        out['N_A'] = float(res['N_A'][i])
    if 'T_free_min' in want:
        out['T_free_min'] = float(np.min(res['T_free'][:i + 1]))
    if 't_c' in want:                          # first time the total potency reaches Z_c
        Z = sum(mf.compute_potency_t(views[pop]) for pop in res['pops'])
        hit = np.nonzero(Z[:i + 1] >= res['params'].Z_c)[0]
        out['t_c'] = float(res['t'][hit[0]]) if len(hit) else np.nan
    if len(res['pops']) == 2:                  # naive + memory
        Y = {pop: POP_OBSERVABLES['Y'](views[pop], i) for pop in res['pops']}
        R = Y['memory'] / Y['naive'] if Y['naive'] > 0 else np.nan
        if 'R' in want:
            out['R'] = R
        if 'phi_naive' in want:
            out['phi_naive'] = 1.0 / (1.0 + cfg['recall']['f'] * R)
        for name in ('R_th', 'kappa'):
            if name in want:
                out[name] = extra.get(name, np.nan)
    return out


def R_theory(cfg, res0, seed, p):
    """{R_th, kappa}: R_th = kappa nu_mu Nbar_w from the naive response (end-of-expansion boost
    prediction); kappa = (pi_c / pi_c_mem)^alpha is the per-cell expansion advantage of memory."""
    v = mf.pop_view(res0, 'naive')
    w, N = v['weights'], v['N_B'][:, -1] * _expanded(v, -1)
    alpha = 1.0 if cfg['model'] == 'two_stage' else p.b0 * (1 / p.lambda_A + 1 / min(p.delta_A, p.b0))
    kappa = (pi_c_of(p) / pi_c_of(p, p.h0_mem)) ** alpha
    if N.sum() <= 0:
        return dict(R_th=np.nan, kappa=kappa)
    Nbar_w = np.sum(w * N ** 2) / np.sum(w * N)              # cell-weighted mean clone size
    nu_mu = np.sum(seed * w) / np.sum(w * N)                  # memory cells per primary cell
    return dict(R_th=float(kappa * nu_mu * Nbar_w), kappa=kappa)


# ============================================================
# Dynamics  (time courses kept per run)
# ============================================================
def response_dynamics(cfg, res, extra=None):
    """Requested time courses: N_A, T_free, and per population Z, Y, L_act, N_B / pi per clone;
    'DG_dist' keeps the final cells per clone (and the memory seed of a recall) for the energy distributions."""
    want, out = cfg['dynamics'], {'t': res['t']}
    for q in ('N_A', 'T_free'):
        if q in want:
            out[q] = res[q]
    for pop in res['pops']:
        v = mf.pop_view(res, pop)
        act = v['N_B'] > 2.0 * v['N_B'][:, :1]
        if 'Z' in want:
            out[f'Z_{pop}'] = mf.compute_potency_t(v)
        if 'Y' in want:
            out[f'Y_{pop}'] = np.sum(v['weights'][:, None] * v['N_B'] * act, axis=0)
        if 'L_act' in want:
            out[f'L_act_{pop}'] = np.sum(v['weights'][:, None] * act, axis=0)
        if 'N_B_clones' in want:
            out[f'N_B_clones_{pop}'] = v['N_B']
        if 'pi_clones' in want:
            out[f'pi_clones_{pop}'] = v['pi']
        if 'DG_dist' in want:
            out[f'N_B_end_{pop}'] = v['N_B'][:, -1]
    if 'N_B_clones' in want or 'pi_clones' in want or 'DG_dist' in want:
        out['DG'] = res['DG']
    if 'DG_dist' in want:
        out['weights'] = res['weights']
        if extra and 'seed' in extra:
            out['seed'] = extra['seed']
    return out


# ============================================================
# Engine
# ============================================================
def sweep_points(cfg):
    names = list(cfg['sweep'])
    return [dict(zip(names, vals)) for vals in itertools.product(*cfg['sweep'].values())]


def run_point(cfg, point):
    """Naive response (+ recall) for one sweep point -> (params, {response: (res, extra)})."""
    p = build_params(cfg, point)
    mode = 'clones' if cfg['repertoire'] == 'two_clone' else cfg['mode']
    kw = dict(t_span=(0.0, cfg['T']), mode=mode, seed=cfg['seed'])
    res0 = mf.run_response(cfg['model'], p, **kw)
    out = {'naive': (res0, {})}
    if 'recall' in cfg['responses']:
        rec = cfg['recall']
        v0 = mf.pop_view(res0, 'naive')
        thr = 2.0 * v0['N_B'][:, 0]
        try:
            if rec['draw'] == 'expected':
                _, seed = mf.memory_seed_expected(v0, n_mem=p.n_mem, threshold=thr)
            else:
                _, seed = mf.memory_seed_from_primary(v0, p=p, n_mem=int(p.n_mem), threshold=thr)
        except ValueError:                     # no clone expanded: no memory, no recall
            print("  [skip recall] no expanded clone in the naive response at "
                  + ", ".join(f"{k}={v:.3g}" for k, v in point.items()))
            return p, out
        drive = dict(t=res0['t'], N_A=res0['N_A']) if rec['antigen'] == 'replay' else None
        res1 = mf.run_response(cfg['model'], p, naive=rec['with_naive'], memory_seed=seed,
                               antigen=rec['antigen'], drive=drive, **kw)
        out['recall'] = (res1, dict(R_theory(cfg, res0, seed, p), seed=seed))
    return p, out


def run(cfg=CONFIG, save=True):
    """Run the configured responses over the sweep; return (observables DataFrame, dynamics)."""
    cfg = full_config(cfg)
    name = run_name(cfg)
    rows, dynamics = [], {}
    points = sweep_points(cfg)
    print(f"[{name}] {len(points)} sweep point(s) x {cfg['responses']}")
    for k, point in enumerate(points):
        p, out = run_point(cfg, point)
        summary = []
        for resp, (res, extra) in out.items():
            dynamics[(tuple(point.items()), resp)] = response_dynamics(cfg, res, extra)
            for t_obs in cfg['t_obs']:
                i = int(np.argmin(np.abs(res['t'] - t_obs)))
                row = dict(point, pi_c=pi_c_of(p), response=resp, t_obs=float(res['t'][i]))
                row.update(response_observables(cfg, res, i, extra))
                rows.append(row)
            summary += [f"{c}[{resp}]={rows[-1][c]:.3g}" for c in ('Y_naive', 'Y_memory', 'R') if c in rows[-1]]
        print(f"  {k + 1:4d}/{len(points)}  " + "  ".join(f"{n}={v:.3g}" for n, v in point.items())
              + f"  pi_c={pi_c_of(p):.3g}  " + "  ".join(summary), flush=True)
    df = pd.DataFrame(rows)
    if save:
        out_dir = os.path.join(DATA_ROOT, name)
        os.makedirs(out_dir, exist_ok=True)
        df.to_csv(os.path.join(out_dir, 'observables.csv'), index=False)
        with open(os.path.join(out_dir, 'dynamics.pkl'), 'wb') as fh:
            pickle.dump(dynamics, fh)
        with open(os.path.join(out_dir, 'config.json'), 'w') as fh:
            json.dump(cfg, fh, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x))
        print(f"saved to {out_dir}")
    return df, dynamics


# ============================================================
# Command line
# ============================================================
# Default grid of each sweepable parameter (--sweep NAME without values)
SWEEP_DEFAULTS = {
    'h0': H0_SWEEP,
    'n_mem': np.logspace(2, 6, 9),
    'N_T0': np.logspace(3, 7, 9),
    't_switch': np.linspace(4.0, 8.0, 9),
    'h_max': np.logspace(2, 4, 5),
    'gap': [0.0, 0.5, 1.0, 2.0, 4.0],
    'ratio': np.logspace(-1, 1, 5),
}
MODEL_DEFAULTS = {'two_stage': dict(T=16.0)}       # expansion starts at t_switch and needs ~10 days


def parse_values(text):
    """'1,2,3' | 'log:a:b:n' (10^a .. 10^b, n points) | 'lin:a:b:n'."""
    if text.startswith(('log:', 'lin:')):
        kind, a, b, n = text.split(':')
        return list((np.logspace if kind == 'log' else np.linspace)(float(a), float(b), int(n)))
    return [float(x) for x in text.split(',')]


def parse_sweep(items, base):
    """--sweep entries -> {name: values}; pi_c is converted to h0 = b0 / pi_c^hill."""
    sweep = {}
    for item in items:
        name, _, text = item.partition('=')
        if name == 'pi_c':
            pics = parse_values(text) if text else list((base['b0'] / np.asarray(H0_SWEEP)) ** (1 / base['hill']))
            sweep['h0'] = [base['b0'] / pc ** base['hill'] for pc in pics]
        elif text:
            sweep[name] = parse_values(text)
        elif name in SWEEP_DEFAULTS:
            sweep[name] = list(SWEEP_DEFAULTS[name])
        else:
            raise SystemExit(f"--sweep {name}: give values, e.g. {name}=1,2,3 or {name}=log:0:2:5")
    return sweep


def cli_config(argv):
    """Config from the command line: CONFIG <- preset <- options. Returns (cfg, args)."""
    ap = argparse.ArgumentParser(description="Run the unified mean-field models; see the module docstring.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('preset', nargs='?', choices=list(PRESETS), help="start from a preset")
    ap.add_argument('--model', choices=list(mf.RESPONSE_MODELS))
    ap.add_argument('--repertoire', choices=['polyclonal', 'two_clone'])
    ap.add_argument('--responses', nargs='+', choices=['naive', 'recall'])
    ap.add_argument('--antigen', choices=['feedback', 'replay'], help="recall antigen law")
    ap.add_argument('--with-naive', action='store_true', default=None, help="include the naive repertoire in the recall")
    ap.add_argument('--no-naive', dest='with_naive', action='store_false', help="recall from memory only")
    ap.add_argument('--draw', choices=['expected', 'random'], help="memory seed")
    ap.add_argument('--sweep', action='append', default=None, metavar='NAME[=VALUES]',
                    help="swept parameter (up to 2): Parameters field, pi_c, or two-clone knob; "
                         "VALUES: 1,2,3 | log:a:b:n | lin:a:b:n (default grid if omitted)")
    ap.add_argument('--single', action='store_true', help="single run (drop any preset sweep)")
    ap.add_argument('--set', action='append', default=[], metavar='PARAM=VALUE',
                    help="override a base parameter or two-clone knob, e.g. --set N_T0=1e4")
    ap.add_argument('--T', type=float, help="end time")
    ap.add_argument('--t-obs', type=float, nargs='+', help="observation times (default T)")
    ap.add_argument('--name', help="run folder name")
    ap.add_argument('--plot', action='store_true', help="make the figures after the run")
    ap.add_argument('--dry-run', action='store_true', help="print the configuration and stop")
    args = ap.parse_args(argv)

    cfg = dict(CONFIG, **(PRESETS[args.preset] if args.preset else {}))
    for key in ('model', 'repertoire', 'responses', 'name'):
        if getattr(args, key) is not None:
            cfg[key] = getattr(args, key)
    recall = dict(cfg['recall'])
    for key, val in (('antigen', args.antigen), ('with_naive', args.with_naive), ('draw', args.draw)):
        if val is not None:
            recall[key] = val
    cfg['recall'] = recall
    if recall['with_naive'] and args.antigen is None and not args.preset:
        recall['antigen'] = 'replay'            # naive + memory: default to the boost setting

    base, tc = dict(cfg.get('base', {})), dict(cfg.get('two_clone', {}))
    for item in args.set:
        k, _, v = item.partition('=')
        (tc if k in TWO_CLONE else base)[k] = float(v)
    cfg['base'], cfg['two_clone'] = {**BASE, **base}, {**TWO_CLONE, **tc}

    if args.single:
        cfg['sweep'] = {}
    elif args.sweep is not None:
        cfg['sweep'] = parse_sweep(args.sweep, cfg['base'])
    elif not args.preset:
        cfg['sweep'] = {}                       # no preset and no --sweep: single run

    cfg.update({k: v for k, v in MODEL_DEFAULTS.get(cfg['model'], {}).items() if args.T is None})
    if args.T is not None:
        cfg['T'] = args.T
    cfg['t_obs'] = args.t_obs or ([cfg['T']] if args.T is not None or cfg['model'] in MODEL_DEFAULTS
                                  or not args.preset else cfg['t_obs'])

    # keep the per-clone time courses that the manuscript figures need
    dyn = list(cfg['dynamics'])
    needs = []
    if not cfg['sweep'] or cfg['repertoire'] == 'two_clone':
        needs = ['N_B_clones', 'pi_clones']
    elif cfg['recall']['with_naive'] and 'n_mem' in cfg['sweep']:
        needs = ['N_B_clones']
    cfg['dynamics'] = dyn + [q for q in needs if q not in dyn]
    try:
        return full_config(cfg), args
    except ValueError as e:
        ap.error(str(e))


def describe(cfg):
    rec = cfg['recall']
    sweep = ', '.join(f"{k} ({len(v)} values, {min(v):.3g} .. {max(v):.3g})" for k, v in cfg['sweep'].items())
    return (f"run       {run_name(cfg)}\n"
            f"model     {cfg['model']}   repertoire {cfg['repertoire']}   responses {cfg['responses']}\n"
            f"recall    antigen={rec['antigen']} with_naive={rec['with_naive']} draw={rec['draw']} h0_mem={rec['h0_mem']}\n"
            f"sweep     {sweep or 'single run'}\n"
            f"T         {cfg['T']}   t_obs {cfg['t_obs']}\n"
            f"dynamics  {cfg['dynamics']}\n"
            f"changed   " + (', '.join(f"{k}={v:g}" for k, v in cfg['base'].items() if v != BASE[k]) or 'none'))


def main(argv=None):
    cfg, args = cli_config(sys.argv[1:] if argv is None else argv)
    print(describe(cfg))
    if args.dry_run:
        return cfg
    run(cfg)
    if args.plot:
        import plot_meanfield
        plot_meanfield.main(run_name(cfg))
    return cfg


if __name__ == '__main__':
    main()
