"""
plot_meanfield.py -- generic figures for a run saved by run_meanfield.py.

    python plot_meanfield.py <run name>          (a folder under run_meanfield.DATA_ROOT)

Writes under run_meanfield.FIG_ROOT/<run>/:
    dynamics/<response>_<quantity>.pdf   time courses; one curve per sweep point
                                         (2-D sweeps: one figure per first / middle / last
                                         value of the second parameter)
    observables/<observable>.pdf         1-D sweep: observable vs the swept parameter,
                                         one line per (response, t_obs)
    observables/<observable>_<response>_t<t_obs>.pdf
                                         2-D sweep: heatmap over the two parameters
    zipf/zipf_<response>.pdf             1-D sweep: rank-size plot of the expanded clones at the end,
                                         one curve per sweep point (single runs: manuscript/zipf.pdf)
A swept h0 is shown as pi_c = (b0 / h0)^(1 / hill).
"""

import os, sys, json, pickle
from dataclasses import fields
sys.path.append(os.path.join(os.path.dirname(__file__), '../../library'))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mf_lib as mf
import mf_plotting as mfp
from run_meanfield import DATA_ROOT, FIG_ROOT

CLONE_ROWS = (0, 10, -10)       # clones drawn from per-clone arrays of a polyclonal run
CLONE_LS = ['-', '--', ':']
NOT_TIME_COURSES = ('t', 'DG', 'weights', 'seed')   # stored with 'DG_dist', not plotted against t


def tex(name):
    """Axis label for a column or quantity name."""
    special = {'pi_c': r'$\pi_c$', 'N_T0': r'$N_T^0$', 'n_mem': r'$n_{\mathrm{mem}}$', 'h0': r'$h_0$',
               't_switch': r'$t_s$', 'h_max': r'$h_{\max}$', 'phi_naive': r'$\phi_{\mathrm{naive}}$',
               'Pi': r'$\Pi$', 'H_T': r'$H_T$', 'K_Pi': r'$K_\Pi$'}
    return special.get(name, name.replace('_', r'\_'))


def axis_name(sweep_name):
    return 'pi_c' if sweep_name == 'h0' else sweep_name


def is_log(values):
    """Log scale when the positive values span > 50x and there are no negatives (zeros are masked)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    pos = v[v > 0]
    return bool(len(pos) > 1 and np.all(v >= 0) and pos.max() / pos.min() > 50)


def point_label(cfg, pt):
    """Legend label for a sweep point: its first parameter (a swept h0 shown as pi_c)."""
    if not pt:
        return None
    name, v = pt[0]
    if name == 'h0':
        b = cfg['base']
        return r'$\pi_c=%.3g$' % (b['b0'] / v) ** (1 / b['hill'])
    return f'{tex(name)}$={v:.3g}$'


def load(run):
    folder = os.path.join(DATA_ROOT, run)
    df = pd.read_csv(os.path.join(folder, 'observables.csv'))
    with open(os.path.join(folder, 'dynamics.pkl'), 'rb') as fh:
        dynamics = pickle.load(fh)
    with open(os.path.join(folder, 'config.json')) as fh:
        cfg = json.load(fh)
    return cfg, df, dynamics


def params_of(base, **over):
    """mf.Parameters from a saved base (fields no longer in Parameters, e.g. of removed models, are dropped)."""
    names = {f.name for f in fields(mf.Parameters)}
    return mf.Parameters(**{k: v for k, v in {**base, **over}.items() if k in names})


def memory_h0(cfg, base, h0):
    """h0 of the memory population: base h0_mem (None -> b0); older runs kept it in cfg['recall']
    ('b0' -> b0, None -> h0)."""
    if 'h0_mem' in cfg['recall']:
        hm = cfg['recall']['h0_mem']
        return base['b0'] if hm == 'b0' else (h0 if hm is None else hm)
    return base['b0'] if base.get('h0_mem') is None else base['h0_mem']


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, f'{name}.pdf'), dpi=800, transparent=True)
    plt.close(fig)


# ============================================================
# Dynamics
# ============================================================
def plot_dynamics(cfg, dynamics, out):
    names = list(cfg['sweep'])
    responses = list(dict.fromkeys(resp for _, resp in dynamics))
    keys = list(dict.fromkeys(point for point, _ in dynamics))
    groups = {'': keys}                                        # figure tag -> sweep points
    if len(names) == 2:                                        # fix the second parameter
        vals = cfg['sweep'][names[1]]
        groups = {f'_{names[1]}_{v:.3g}': [k for k in keys if dict(k)[names[1]] == v]
                  for v in [vals[0], vals[len(vals) // 2], vals[-1]]}
    for tag, pts_all in groups.items():
        for resp in responses:
            pts = [pt for pt in pts_all if (pt, resp) in dynamics]   # a recall can be missing
            if not pts:
                continue
            cols = mfp.sweep_colors(len(pts), 'turbo') if len(pts) > 1 else ['k']
            quantities = [q for q in dynamics[(pts[0], resp)] if q not in NOT_TIME_COURSES
                          and not q.startswith('N_B_end')]
            for q in quantities:
                fig, ax = mfp.new_fig()
                for c, pt in zip(cols, pts):
                    d = dynamics[(pt, resp)]
                    label = point_label(cfg, pt)
                    y = np.asarray(d[q])
                    if y.ndim == 2:                            # per-clone array
                        rows = range(len(y)) if len(y) <= 3 else CLONE_ROWS
                        for ls, r in zip(CLONE_LS, rows):
                            ax.plot(d['t'], y[r], color=c, ls=ls, lw=3, label=label if r == 0 else None)
                    else:
                        ax.plot(d['t'], y, color=c, lw=3, label=label)
                if q not in ('T_free',) and not q.startswith('f_'):
                    ax.set_yscale('log')
                    top = ax.get_ylim()[1]
                    ax.set_ylim(bottom=max(ax.get_ylim()[0], top * 1e-12))   # hide near-empty bins
                ax.set_xlabel(r'$t$', fontsize=30)
                ax.set_ylabel(tex(q), fontsize=30)
                ax.tick_params(axis='both', labelsize=26)
                if len(pts) > 1 and len(pts) <= 16:
                    ax.legend(fontsize=12, ncol=2)
                save(fig, os.path.join(out, 'dynamics'), f'{resp}_{q}{tag}')


# ============================================================
# Observables
# ============================================================
VALUE_COLS = ('response', 't_obs')


def observable_columns(cfg, df):
    meta = set(cfg['sweep']) | {'pi_c', 'response', 't_obs'}
    return [c for c in df.columns if c not in meta and df[c].notna().any()]


def plot_observables_1d(cfg, df, out):
    x = axis_name(list(cfg['sweep'])[0])
    combos = df[list(VALUE_COLS)].drop_duplicates().itertuples(index=False)
    combos = list(combos)
    cols = mfp.sweep_colors(len(combos), 'viridis') if len(combos) > 1 else ['k']
    for obs in observable_columns(cfg, df):
        fig, ax = mfp.new_fig()
        for c, (resp, t_obs) in zip(cols, combos):
            sub = df[(df.response == resp) & (df.t_obs == t_obs)].sort_values(x)
            if sub[obs].notna().any():
                ax.plot(sub[x], sub[obs], 'o-', color=c, mec='k', ms=9, lw=2,
                        label=fr'{resp}, $t={t_obs:.3g}$')
        if not ax.lines:
            plt.close(fig)
            continue
        if is_log(df[x]):
            ax.set_xscale('log')
        if is_log(df[obs].dropna()):
            ax.set_yscale('log')
        ax.set_xlabel(tex(x), fontsize=30)
        ax.set_ylabel(tex(obs), fontsize=30)
        ax.tick_params(axis='both', labelsize=26)
        ax.legend(fontsize=16)
        save(fig, os.path.join(out, 'observables'), obs)


def plot_observables_2d(cfg, df, out):
    n1, n2 = list(cfg['sweep'])
    x, y = axis_name(n1), n2
    for obs in observable_columns(cfg, df):
        for (resp, t_obs), sub in df.groupby(list(VALUE_COLS)):
            if not sub[obs].notna().any():
                continue
            grid = sub.pivot_table(index=y, columns=x, values=obs)
            Z = grid.to_numpy()
            log_z = is_log(Z[np.isfinite(Z)])
            fig, ax = mfp.new_fig()
            with np.errstate(divide='ignore'):
                Zp = np.where(Z > 0, np.log10(np.where(Z > 0, Z, 1.0)), np.nan) if log_z else Z
            mesh = ax.pcolormesh(grid.columns, grid.index, Zp,
                                 shading='nearest', cmap='viridis')
            fig.colorbar(mesh, ax=ax).set_label((r'$\log_{10}$ ' if log_z else '') + tex(obs), fontsize=24)
            if is_log(grid.columns):
                ax.set_xscale('log')
            if is_log(grid.index):
                ax.set_yscale('log')
            ax.set_xlabel(tex(x), fontsize=30)
            ax.set_ylabel(tex(y), fontsize=30)
            ax.tick_params(axis='both', labelsize=26)
            ax.set_title(fr'{resp}, $t={t_obs:.3g}$', fontsize=22)
            save(fig, os.path.join(out, 'observables'), f'{obs}_{resp}_t{t_obs:.3g}')


# ============================================================
# Energy distributions ('DG_dist' dynamics)
# ============================================================
def _draw_DG_dist(ax, cfg, d0, d1, fs=1.0):
    """Cell density (cells per unit DG) of the expanded clones: primary end, memory seed, recall end.
    Theory shapes over each range, scaled to the same total cells: EP exponential Omega_0 e^{-a DG}
    (a = alpha, or 2 alpha for the recall) and the tilted Gaussian of Fig. 4D for the recall."""
    b = cfg['base']
    alpha = b['b0'] * (1.0 / b['lambda_A'] + 1.0 / min(b['delta_A'], b['b0']))
    DG, w = d0['DG'], d0['weights']
    dDG = DG[1] - DG[0]
    curves = [(r'$\mathrm{primary}$', d0['N_B_end_naive'], d0['N_B_end_naive'] > 2.0, mfp.my_orange, alpha)]
    if d1 is not None and 'seed' in d1:
        curves.append((r'$\mathrm{memory}$', d1['seed'], d1['seed'] > 0, mfp.my_purple, alpha))
        if 'N_B_end_memory' in d1:
            n = d1['N_B_end_memory']
            curves.append((r'$\mathrm{recall}$', n, n > 2.0 * d1['seed'], mfp.my_blue, 2 * alpha))
    for label, n, m, c, a in curves:
        if m.sum() < 2:
            continue
        rho = w[m] * n[m] / dDG
        total = np.sum(w[m] * n[m])
        seed = label == r'$\mathrm{memory}$'            # the primary's shape, scaled by n_draw / Y
        ax.plot(DG[m], rho, 'o--' if seed else 'o-', color=c, ms=(11 if seed else 8) * fs, lw=2 * fs,
                mfc='none' if seed else c, mew=2 * fs, label=label, zorder=1 if seed else 2)
        x = np.linspace(DG[m].min(), DG[m].max(), 200)
        ep = np.exp((b['beta_star'] - a) * (x - x[0]))
        ax.plot(x, total * ep / np.trapezoid(ep, x), ':', color=c, lw=3 * fs)
        if label == r'$\mathrm{recall}$':
            g = np.exp(-(x - x[0] + DG4D_STAR - FIG4D_MU_RECALL) ** 2 / (2 * DG4D_SIG2))
            ax.plot(x, total * g / np.trapezoid(g, x), '--', color='grey', lw=2 * fs, label=r'$\mathrm{Gaussian}$')
    ax.set_yscale('log')


def plot_DG_dist(cfg, dynamics, out):
    """One figure per sweep point (and a grid for 1-D sweeps) of the cell densities in DG."""
    pts = [pt for pt in dict.fromkeys(pt for pt, _ in dynamics)
           if 'N_B_end_naive' in dynamics.get((pt, 'naive'), {})]
    if not pts or len(cfg['sweep']) > 1:
        return
    o = os.path.join(out, 'DG_dist')
    for pt in pts:
        fig, ax = mfp.new_fig()
        _draw_DG_dist(ax, cfg, dynamics[(pt, 'naive')], dynamics.get((pt, 'recall')))
        ax.set_xlabel(r'$\Delta G$', fontsize=30)
        ax.set_ylabel(r'$\mathrm{cells\ per\ unit}\ \Delta G$', fontsize=30)
        ax.tick_params(axis='both', labelsize=26)
        if pt:
            ax.set_title(point_label(cfg, pt), fontsize=24)
        ax.legend(fontsize=18)
        b = cfg['base']
        tag = '_'.join(f'pi_c_{(b["b0"] / v) ** (1 / b["hill"]):.3g}' if k == 'h0' else f'{k}_{v:.3g}' for k, v in pt)
        save(fig, o, 'DG_dist_' + (tag or 'single'))
    if len(pts) > 1:
        ncol = int(np.ceil(np.sqrt(len(pts))))
        nrow = int(np.ceil(len(pts) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.2 * nrow), sharex=True, sharey=True)
        for ax, pt in zip(axs.flat, pts):
            _draw_DG_dist(ax, cfg, dynamics[(pt, 'naive')], dynamics.get((pt, 'recall')), fs=0.5)
            ax.set_title(point_label(cfg, pt), fontsize=14)
        for ax in axs.flat[len(pts):]:
            ax.axis('off')
        axs.flat[0].legend(fontsize=9)
        fig.supxlabel(r'$\Delta G$', fontsize=18)
        fig.supylabel(r'$\mathrm{cells\ per\ unit}\ \Delta G$', fontsize=18)
        fig.tight_layout()
        save(fig, o, 'DG_dist_grid')


# ============================================================
# Zipf (rank-size) plots ('DG_dist' dynamics)
# ============================================================
def zipf_curve(d, pop, n_ranks=200):
    """Rank-size curve of the expanded clones at the end of a response (doubled since t = 0, as L).
    Each DG bin holds `weights` clones of equal size, so rank k is the size of the bin where the
    cumulative clone count reaches k. Returns (ranks, sizes / size at rank 1), or None if < 1 clone."""
    n_end = d[f'N_B_end_{pop}']
    n0 = d['seed'] if pop == 'memory' else d[f'N_B_clones_{pop}'][:, 0] if f'N_B_clones_{pop}' in d \
        else np.ones_like(n_end)
    m = (n0 > 0) & (n_end > 2.0 * n0)
    order = np.argsort(n_end[m])[::-1]
    sizes, cum = n_end[m][order], np.cumsum(d['weights'][m][order])
    if not len(cum) or cum[-1] < 1:
        return None
    ranks = np.unique(np.logspace(0, np.log10(cum[-1]), n_ranks))
    s = sizes[np.minimum(np.searchsorted(cum, ranks), len(sizes) - 1)]
    return ranks, s / s[0]


def plot_zipf(cfg, dynamics, out):
    """Rank-size plots of the expanded clones at the end of the primary and the recall, with the
    power laws k^-zeta (primary) and k^-2 zeta (recall), zeta = eta (b0/lambda_A + b0/min(b0, delta_A)) / beta*.
    Single runs: one figure with both responses, in <run>/manuscript/zipf.pdf.
    1-D sweeps: one figure per response, one curve per sweep point, in <run>/zipf/."""
    if 'DG_dist' not in cfg['dynamics'] or len(cfg['sweep']) > 1:
        return
    b = cfg['base']
    zeta = b['eta'] * (b['b0'] / b['lambda_A'] + b['b0'] / min(b['b0'], b['delta_A'])) / b['beta_star']
    resps = [(r, pop, c, k * zeta) for r, pop, c, k in
             [('naive', 'naive', mfp.my_orange, 1), ('recall', 'memory', mfp.my_blue, 2)]]

    def style(ax, curves, laws):
        top = max(r[-1] for r, _ in curves)
        k = np.logspace(0, np.log10(top), 50)
        for c, z in laws:
            ax.plot(k, k ** -z, '--', color=c, lw=2)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(bottom=max(ax.get_ylim()[0], 1e-4), top=1.5)
        ax.set_xlabel(r'$\mathrm{rank}\ k$', fontsize=30)
        ax.set_ylabel(r'$N_k / N_1$', fontsize=30)
        ax.tick_params(axis='both', labelsize=26)

    pts = list(dict.fromkeys(pt for pt, _ in dynamics))
    if not cfg['sweep']:
        fig, ax = mfp.new_fig()
        curves = []
        for r, pop, c, z in resps:
            zc = zipf_curve(dynamics[((), r)], pop) if ((), r) in dynamics else None
            if zc is not None:
                ax.plot(*zc, '-', color=c, lw=3,
                        label=r'$\mathrm{%s},\ %s=%.2g$' % ('primary' if r == 'naive' else 'recall',
                                                          r'\zeta' if r == 'naive' else r'2\zeta', z))
                curves.append(zc)
        if curves:
            style(ax, curves, [(c, z) for _, _, c, z in resps])
            ax.legend(fontsize=20)
            save(fig, os.path.join(out, 'manuscript'), 'zipf')
        else:
            plt.close(fig)
        return
    cols = mfp.sweep_colors(len(pts), 'turbo')
    for r, pop, _, z in resps:
        fig, ax = mfp.new_fig()
        curves = []
        for c, pt in zip(cols, pts):
            zc = zipf_curve(dynamics[(pt, r)], pop) if (pt, r) in dynamics else None
            if zc is not None:
                ax.plot(*zc, '-', color=c, lw=2, label=point_label(cfg, pt))
                curves.append(zc)
        if curves:
            style(ax, curves, [('grey', z)])
            if len(curves) <= 16:
                ax.legend(fontsize=12, ncol=2)
            save(fig, os.path.join(out, 'zipf'), f'zipf_{r}')
        else:
            plt.close(fig)


# ============================================================
# Manuscript figures (ported from run_meanfield_legacy.py), under <run>/manuscript/
# ============================================================
DG4D_STAR, DG4D_SIG2 = -16.0, 7.0     # naive N(0, sigma2) with high-affinity edge -beta* sigma2
FIG4D_MU_RECALL = -16.0               # tilt of the recall responding Gaussian (zeta ~ 1)
P_T_REF = 4.4                         # potency P = lambda_A (P_T_REF - t_c), as mf.compute_potency_P
BOOST_DATA = dict(flu_HA_d6=0.044, lo=0.038, hi=0.11)   # Mesin et al. 2020: memory fraction
MESIN_S2_RECALL = dict(CGG=1.76, HA=1.93)   # Mesin et al. 2020, recall GC+PB: -ln(mean lambda2) (zipf.py)


def find_dynamics(dynamics, point, resp):
    """Time courses of one (sweep point, response); values matched with a tolerance
    (they come back from observables.csv slightly rounded)."""
    for (key, r), d in dynamics.items():
        if r == resp and len(key) == len(point) and all(
                n == k and np.isclose(v, point[k], rtol=1e-9, atol=0) for (n, v), k in zip(key, point)):
            return d
    raise KeyError((point, resp))


def last_t(df, resp):
    sub = df[df.response == resp]
    return sub[sub.t_obs == sub.t_obs.max()]


def fig_potency(cfg, df, out):
    """Recall potency vs specificity across pi_c, with the analytic curve (legacy 'P_shared' -> potency.pdf).

    Specificity: analytic D_KL of the tilted, truncated Gaussian between the recall's activated
    DG range. Theory (SI, optimal threshold): lambda_A t_c = (1/r)[-alpha E_max + log(1 + Y/M)],
    reach E_max and yield Y of the naive response fitted across pi_c, constant fitted to the points.
    """
    b = cfg['base']
    lamA, lamB, dA, beta, M = b['lambda_A'], b['b0'], b['delta_A'], b['beta_star'], b['n_mem']
    alpha = lamB * (1.0 / lamA + 1.0 / min(dA, lamB))
    r = lamB / lamA
    rec, nav = last_t(df, 'recall'), last_t(df, 'naive')

    # recall points, coloured by sweep position (legacy 'turbo_r')
    cols = mfp.sweep_colors(len(rec), 'turbo_r')
    pts = []
    for c, row in zip(cols, rec.itertuples()):
        P = lamA * (P_T_REF - row.t_c)
        D_ana = mf.dkl_analytic(FIG4D_MU_RECALL, DG4D_STAR + row.DG_hi_memory,
                                DG_star=DG4D_STAR + row.DG_lo_memory, sigma2=DG4D_SIG2)
        if not (np.isnan(row.D_KL_memory) or np.isnan(P) or np.isnan(D_ana)):
            pts.append((D_ana, P, c, row.pi_c))
    D_pts, P_pts = np.array([q[0] for q in pts]), np.array([q[1] for q in pts])
    pic_pts = np.array([q[3] for q in pts])

    # theory curve from the naive responses
    nav = nav.dropna(subset=['DG_hi_naive', 'DG_lo_naive', 'Y_naive'])
    pic_sim = nav.pi_c.to_numpy()
    Emax_sim = (nav.DG_hi_naive - nav.DG_lo_naive).to_numpy()
    DGmin = nav.DG_lo_naive.iloc[-1]
    E0 = np.mean(Emax_sim + np.log(pic_sim))
    Emax_of = lambda pic: E0 - np.log(pic)
    logY_shape = lambda Em: alpha * Em + np.log(np.expm1((beta - alpha) * Em) / (beta - alpha))
    logY0 = np.mean(np.log(nav.Y_naive.to_numpy()) - logY_shape(Emax_of(pic_sim)))
    P_shape = lambda pic: (alpha * Emax_of(pic) - np.log1p(np.exp(logY0 + logY_shape(Emax_of(pic))) / M)) / r
    const = np.mean(P_pts - P_shape(pic_pts))
    pic = np.logspace(np.log10(pic_sim.min()), np.log10(pic_sim.max()), 200)
    D_th = np.array([mf.dkl_analytic(FIG4D_MU_RECALL, DG4D_STAR + e + DGmin, DG_star=DG4D_STAR + DGmin,
                                     sigma2=DG4D_SIG2) for e in Emax_of(pic)])

    P_th = P_shape(pic) + const
    ylim = (np.min(P_pts) - 1.0, max(np.max(P_pts), np.nanmax(P_th)) + 1.0)   # P is defined up to P_T_REF
    fig, ax = mfp.new_fig()
    idx = np.argsort(D_pts)
    ax.scatter(np.exp(D_pts[idx]), P_pts[idx], facecolors=[pts[i][2] for i in idx], edgecolors='k',
               marker='s', s=130, linewidth=1.5, zorder=2, label=r'$\mathrm{recall}$')
    ax.plot(np.exp(D_th), P_shape(pic) + const, color='k', lw=2, ls='--', zorder=1, label=r'$\mathrm{theory}$')
    ax.set_xlim(np.exp(D_pts.min()) / 1.5, np.exp(D_pts.max()) * 1.5)
    ax.set_ylim(*ylim)
    ax.set_xscale('log')
    ax.tick_params(axis='both', labelsize=44)
    ax.legend(fontsize=24, loc=0)
    save(fig, out, 'potency')

    # consistency check: potency vs pi_c, simulation and theory. Second curve: the same recall step
    # (best memory clone N_1 / (1 + Y/M), finite reach) fed with the measured primary N_1, Y and reach.
    k = alpha + 1.0 + r * b['eta'] - beta
    nav_pic = nav.set_index('pi_c').reindex(pic_pts)
    E_pts = (nav_pic.DG_hi_naive - nav_pic.DG_lo_naive).to_numpy()
    P_meas = (np.log(nav_pic.N1_naive.to_numpy()) - np.log1p(nav_pic.Y_naive.to_numpy() / M)
              + np.log(-np.expm1(-k * E_pts))) / r
    ok = np.isfinite(P_meas)
    P_meas = P_meas + np.mean(P_pts[ok] - P_meas[ok])
    fig, ax = mfp.new_fig()
    ax.scatter(pic_pts, P_pts, facecolors=[q[2] for q in pts], edgecolors='k', marker='s', s=130,
               linewidth=1.5, zorder=2, label=r'$\mathrm{recall}$')
    ax.plot(pic, P_th, color='k', lw=2, ls='--', zorder=1, label=r'$\mathrm{theory}$')
    ax.plot(pic_pts[ok], P_meas[ok], color='grey', lw=2, ls=':', zorder=1,
            label=r'$\mathrm{theory,\ measured}\ N_1, Y$')
    ax.set_xscale('log')
    ax.set_ylim(*ylim)
    ax.set_xlabel(r'$\pi_c$', fontsize=36)
    ax.set_ylabel(r'$P$', fontsize=36)
    ax.tick_params(axis='both', labelsize=44)
    ax.legend(fontsize=24, loc=0)
    save(fig, out, 'potency_pi_c')

    # same potencies against the Renyi-2 clone entropy, with the data (zipf.py) as a band.
    # Like D, S2 is evaluated analytically on a reach (monotonic in pi_c): recall clone sizes
    # ~ exp(-2 alpha DG) (primary EP x recall EP) over the recall's DG range (points) or E_max(pi_c) (theory).
    S2_of = lambda E: mf.renyi2_analytic(2 * alpha, E, beta, b['omega_0'])
    rec_pic = rec.set_index('pi_c').reindex(pic_pts)
    S2_pts = np.array([S2_of(E) for E in (rec_pic.DG_hi_memory - rec_pic.DG_lo_memory)])
    S2_th = np.array([S2_of(e) for e in Emax_of(pic)])
    S2_sim = rec_pic.S2_memory.to_numpy() if 'S2_memory' in rec_pic else None
    variants = [('potency_entropy', S2_pts, True)]
    if S2_sim is not None:               # check: S2 measured on the simulated clone sizes
        variants.append(('potency_entropy_sim', S2_sim, False))
    for name, x, with_theory in variants:
        fig, ax = mfp.new_fig()
        ax.axvspan(min(MESIN_S2_RECALL.values()), max(MESIN_S2_RECALL.values()), color='grey', alpha=0.3,
                   label=r'$\mathrm{data}$')
        ax.scatter(x[idx], P_pts[idx], facecolors=[pts[i][2] for i in idx], edgecolors='k',
                   marker='s', s=130, linewidth=1.5, zorder=2, label=r'$\mathrm{recall}$')
        if with_theory:
            ax.plot(S2_th, P_shape(pic) + const, color='k', lw=2, ls='--', zorder=1, label=r'$\mathrm{theory}$')
        ax.set_ylim(*ylim)
        ax.set_xlabel(r'$S_2\,[\mathrm{nats}]$', fontsize=36)
        ax.tick_params(axis='both', labelsize=44)
        ax.legend(fontsize=24, loc=0)
        save(fig, out, name)
    return dict(D=D_pts, P=P_pts, pi_c=pic_pts, x_th=np.exp(D_th), P_th=P_shape(pic) + const, P_meas=P_meas,
                S2=S2_pts, S2_th=S2_th, S2_sim=S2_sim)


def fig_boost(cfg, df, dynamics, out):
    """R vs pi_c, phi_naive map, r(DG) and R vs n_mem, one folder per t_obs (legacy boost figures)."""
    f = cfg['recall']['f']
    R_of = lambda m: m / (1.0 - m) / f                  # memory fraction -> R (inverts phi_naive)
    nmem_label = lambda n: r'$10^{%.1f}$' % np.log10(n)
    Y = last_t(df, 'naive').set_index(['h0', 'n_mem']).Y_naive          # primary yield (end of run)
    for t_obs, rec in df[df.response == 'recall'].groupby('t_obs'):
        o = os.path.join(out, f't_read_{round(t_obs, 1):g}')
        pics = np.sort(rec.pi_c.unique())
        n_mems = np.sort(rec.n_mem.unique())
        grid = {k: rec.pivot_table(index='n_mem', columns='pi_c', values=k).reindex(index=n_mems, columns=pics)
                for k in ('R', 'R_th', 'phi_naive')}
        cols = mfp.sweep_colors(len(n_mems), 'winter')
        rep = [0, len(pics) // 2, len(pics) - 1] if len(pics) > 2 else list(range(len(pics)))
        rep_cols = mfp.sweep_colors(len(rep), 'autumn')
        n_rep = n_mems[len(n_mems) // 2]

        # 1. R vs pi_c
        fig, ax = mfp.new_fig()
        for c, n in zip(cols, n_mems):
            ax.plot(pics, grid['R'].loc[n], 'o', color=c, mec='k', ms=10, label=nmem_label(n))
            ax.plot(pics, grid['R_th'].loc[n], '-', color=c, lw=2)
        ax.axhspan(R_of(BOOST_DATA['lo']), R_of(BOOST_DATA['hi']), color='grey', alpha=0.3)
        ax.axhline(R_of(BOOST_DATA['flu_HA_d6']), color='k', ls='--', lw=1.5)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$\pi_c$', fontsize=30)
        ax.set_ylabel(r'$R$', fontsize=30)
        ax.tick_params(axis='both', labelsize=30)
        ax.legend(fontsize=16, title=r'$n_{\mathrm{mem}}$', title_fontsize=18, ncol=2)
        save(fig, o, 'R_pi')

        # 2. phi_naive over (pi_c, n_mem), contours at the data memory fractions
        phi = grid['phi_naive'].to_numpy()
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
        save(fig, o, 'phi_naive_map')

        # 3. r(DG) = memory / naive fold expansion, with kappa
        fig, ax = mfp.new_fig()
        for c, k in zip(rep_cols, rep):
            row = rec[(rec.pi_c == pics[k]) & (rec.n_mem == n_rep)].iloc[0]
            d = find_dynamics(dynamics, {'h0': row.h0, 'n_mem': row.n_mem}, 'recall')
            i = int(np.argmin(np.abs(d['t'] - t_obs)))
            Nm, Nn = d['N_B_clones_memory'], d['N_B_clones_naive']
            with np.errstate(divide='ignore', invalid='ignore'):
                r = np.where(Nm[:, 0] > 0, Nm[:, i] / Nm[:, 0] / Nn[:, i], np.nan)
            ax.plot(d['DG'], r, 'o-', color=c, mec='k', ms=8, lw=2, label=mfp.pi_star_label(pics[k]))
            ax.axhline(row.kappa, color=c, ls='--', lw=1.5)
        ax.set_yscale('log')
        ax.set_xlabel(r'$\Delta G$', fontsize=30)
        ax.set_ylabel(r'$r(\Delta G)$', fontsize=30)
        ax.tick_params(axis='both', labelsize=30)
        ax.legend(fontsize=18, title=r'$\pi_c$', title_fontsize=20)
        save(fig, o, 'r_DG')

        # 4. R vs n_mem at fixed pi_c; saturation expected at n_mem ~ Y (dotted)
        fig, ax = mfp.new_fig()
        for c, k in zip(rep_cols, rep):
            ax.plot(n_mems, grid['R'][pics[k]], 'o', color=c, mec='k', ms=10, label=mfp.pi_star_label(pics[k]))
            ax.plot(n_mems, grid['R_th'][pics[k]], '-', color=c, lw=2)
            h0 = rec[rec.pi_c == pics[k]].h0.iloc[0]
            ax.axvline(Y.loc[(h0, n_mems[0])], color=c, ls=':', lw=2)
        ax.axhspan(R_of(BOOST_DATA['lo']), R_of(BOOST_DATA['hi']), color='grey', alpha=0.3)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(n_mems.min() / 2, n_mems.max() * 2)
        ax.set_xlabel(r'$n_{\mathrm{mem}}$', fontsize=30)
        ax.set_ylabel(r'$R$', fontsize=30)
        ax.tick_params(axis='both', labelsize=30)
        ax.legend(fontsize=18, title=r'$\pi_c$', title_fontsize=20)
        save(fig, o, 'R_nmem')


def fig_two_clone(cfg, df, dynamics, out):
    """Clone sizes, clone-1 fraction, pi and final clone-1 fraction vs the swept knob (legacy two-clone figures).
    Solid = clone 1 (DG1), dashed = clone 2 (DG1 + gap); colour = sweep value."""
    name = list(cfg['sweep'])[0]
    vals = cfg['sweep'][name]
    T, b = cfg['T'], cfg['base']
    cols = mfp.sweep_colors(len(vals), 'autumn')
    name_tex = name.replace('_', r'\_')
    responses = list(dict.fromkeys(resp for _, resp in dynamics))
    final = {}
    for resp in responses:
        pop = 'naive' if resp == 'naive' else 'memory'
        fig_NB, ax_NB = mfp.new_fig()
        fig_f, ax_f = mfp.new_fig()
        fig_pi, ax_pi = mfp.new_fig()
        for c, v in zip(cols, vals):
            try:
                d = find_dynamics(dynamics, {name: v}, resp)
            except KeyError:
                continue
            t, N_B, pi = d['t'], d[f'N_B_clones_{pop}'], d[f'pi_clones_{pop}']
            lab = r'$\mathrm{%s} = %.3g$' % (name_tex, v)
            ax_NB.plot(t, N_B[0], lw=5, color=c, ls='-', label=lab)
            ax_NB.plot(t, N_B[1], lw=5, color=c, ls='--')
            f1 = N_B[0] / (N_B[0] + N_B[1])
            ax_f.plot(t, f1, lw=5, color=c, label=lab)
            ax_pi.plot(t, pi[0], lw=5, color=c, ls='-', label=lab)
            ax_pi.plot(t, pi[1], lw=5, color=c, ls='--')
            h0 = v if name == 'h0' else b['h0']
            if pop == 'memory':
                h0 = v if name == 'h0_mem' else memory_h0(cfg, b, h0)
            ax_pi.axhline((b['b0'] / h0) ** (1 / b['hill']), lw=1, ls=':', color=c)
            final.setdefault(resp, []).append((v, f1[-1]))
        tag = f"{cfg['model']}_{resp}"
        mfp.style_log_axis(ax_NB, T, ylim=(5e-1, None))
        ax_NB.legend(fontsize=20)
        save(fig_NB, out, f'N_B_{tag}')
        ax_f.set_xlim(0, T)
        ax_f.set_ylim(-0.02, 1.02)
        ax_f.tick_params(axis='both', labelsize=44)
        ax_f.legend(fontsize=20)
        save(fig_f, out, f'frac_clone1_{tag}')
        mfp.style_log_axis(ax_pi, T, ylim=(1e-2, None))
        ax_pi.legend(fontsize=20)
        save(fig_pi, out, f'pi_{tag}')

    if len(vals) > 1:
        fig, ax = mfp.new_fig()
        for k, (resp, pts) in enumerate(final.items()):
            x, y = np.array(pts).T
            ax.plot(x, y, lw=3, ls=mfp.styles_mem[k], color=mfp.colors_mem[k], marker='os'[k], ms=14,
                    mec='k', label=r'$\mathrm{%s}$' % ('primary' if resp == 'naive' else 'recall'))
            if np.all(x > 0) and x.max() / x.min() > 100:
                ax.set_xscale('log')
        ax.set_xlabel(r'$\mathrm{%s}$' % name_tex, fontsize=30)
        ax.set_ylim(-0.02, 1.02)
        ax.tick_params(axis='both', labelsize=30)
        ax.legend(fontsize=20)
        save(fig, out, 'final_frac_clone1')


def fig_ep(cfg, dynamics, out):
    """Fig. 2 panels for a single run, one folder per response (legacy per-run figures, no null model):
    N_A (A), pi of three clones (B), N_B of three clones (C), potency Z, and the affinity fronts (D;
    analytic continuous-model fronts, so only for model='continuous')."""
    T, base = cfg['T'], cfg['base']
    p0 = params_of(base)
    weights = mf.build_context(p0, mode='grid')['weights']
    pi_c_primary = (p0.b0 / p0.h0) ** (1 / p0.hill)
    label = mfp.pi_star_label(pi_c_primary)                 # both responses labelled by the primary pi_c
    for k, resp in enumerate(['naive', 'recall']):
        if ((), resp) not in dynamics:
            continue
        d = dynamics[((), resp)]
        pop = 'naive' if resp == 'naive' else 'memory'
        h0 = memory_h0(cfg, base, p0.h0) if pop == 'memory' else p0.h0
        res = dict(t=d['t'], N_A=d['N_A'], pi=d[f'pi_clones_{pop}'], N_B=d[f'N_B_clones_{pop}'], Z=d[f'Z_{pop}'],
                   DG=d['DG'], weights=weights, model=cfg['model'],
                   params=params_of(base, h0=h0, memory=(pop == 'memory')))
        color = mfp.colors_mem[k]
        o = os.path.join(out, resp)
        left = dict(labelleft=(k == 0))                     # y tick labels on the primary only

        fig, ax = mfp.new_fig(wide=False)
        mfp.plot_NA(ax, res, label=label)
        mfp.style_log_axis(ax, T - 1, ylim=(1e0, 1e11))
        ax.tick_params(axis='y', which='both', labelsize=44, **left)
        ax.tick_params(axis='x', labelsize=44, labelbottom=False)
        save(fig, o, 'N_A')

        fig, ax = mfp.new_fig(wide=False)
        mfp.plot_pi(ax, res, color=color, label=label)
        ax.axhline((p0.b0 / h0) ** (1 / p0.hill), lw=1, ls='--', color='grey', alpha=0.8)
        mfp.style_log_axis(ax, T - 1, ylim=(5e-1, 1e4))
        ax.tick_params(axis='y', which='both', labelsize=44, **left)
        ax.tick_params(axis='x', labelsize=44, labelbottom=False)
        save(fig, o, 'pi')

        fig, ax = mfp.new_fig(wide=False)
        mfp.plot_NB(ax, res, color=color, label=label)
        ax.axhline(1.0, color='k', ls='--', alpha=0.5)
        mfp.style_log_axis(ax, T - 1, ylim=(5e-1, 1e6))
        ax.tick_params(axis='y', which='both', labelsize=44, **left)
        ax.tick_params(axis='x', labelsize=44, labelbottom=False)
        save(fig, o, 'N_B')

        fig, ax = mfp.new_fig(wide=False)
        mfp.plot_potency(ax, res, color=color, label=label)
        ax.axhline(p0.Z_c, lw=1, ls='--', color='k')
        mfp.style_log_axis(ax, T - 1, ylim=(5e-3, 1e7), hide_xticklabels=False)
        save(fig, o, 'Z')

        if cfg['model'] == 'continuous':
            fr = mfp.compute_affinity_front(res, on_slope=0.85)
            fig, ax = mfp.new_fig(wide=False)
            ax.plot(fr['ton'], fr['DG_on'], lw=6, color=color, ls='-')
            off = fr['DG_off'] > p0.DG_min
            ax.plot(fr['toff'][off], fr['DG_off'][off], lw=6, color=color, ls='dashed')
            ax.axhline(p0.DG_min, lw=1, ls='--', color='grey', alpha=0.8)
            ax.set_ylim(p0.DG_min - 0.1, 7)
            ax.set_xlim(0, T - 1)
            ax.tick_params(axis='y', which='both', labelsize=44, **left)
            ax.tick_params(axis='x', labelsize=44, labelbottom=True)
            save(fig, o, 'DG')


SNAPSHOT_TIMES = {'naive': (3.5, 4.5, 5.5, 6.5),        # panel A, primary (days)
                  'recall': (3, 4, 5, 6)}      # panel B, recall (days)
COLOR_L, COLOR_Z, COLOR_Y = mfp.my_red, mfp.my_green, mfp.my_purple
GOLDEN = (1 + 5 ** 0.5) / 2                     # width / height of the fig_energy_dynamics panels


def _energy_panels(cfg, d, resp, pop, axA, axC, zaxC):
    """Draw one response's panels of fig_energy_dynamics: cells per unit DG of activated clones at the
    snapshot times (axA), and L (axC) with Z (zaxC, right axis)."""
    T = cfg['T']
    t, DG, w = d['t'], d['DG'], d['weights']
    dDG = DG[1] - DG[0]
    N_B = d[f'N_B_clones_{pop}']
    times = [ts for ts in SNAPSHOT_TIMES[resp] if ts <= T]
    tcols = plt.get_cmap('Blues')(np.linspace(0.35, 1.0, len(times)))
    for c, ts in zip(tcols, times):
        i = int(np.argmin(np.abs(t - ts)))
        n = N_B[:, i]
        axA.plot(DG, np.where((N_B[:, 0] > 0) & (n > 2.0 * N_B[:, 0]), w * n / dDG, np.nan), 'o-', color=c, ms=4, lw=2, label=r'$t=%g$' % ts)
    axA.set_yscale('log')
    axA.set_xlim(right=5)
    axA.set_xlabel(r'$\Delta G$', fontsize=22)

    L, Y, Z = d[f'L_act_{pop}'], d[f'Y_{pop}'], d[f'Z_{pop}']
    axC.plot(t, np.where(L > 0, L, np.nan), color=COLOR_L, lw=3, label=r'$L$')
    axC.plot(t, np.where(Y > 0, Y, np.nan), color=COLOR_Y, lw=3, label=r'$Y$')
    zaxC.plot(t, np.where(Z > 1e-2, Z, np.nan), color=COLOR_Z, lw=3, label=r'$Z$')   # hide the switch's leak
    zaxC.axhline(cfg['base']['Z_c'], color=COLOR_Z, ls='--', lw=1)
    zaxC.set_yscale('log')
    axC.set_yscale('log')

    # for c, ts in zip(tcols, times):
    #     axC.axvline(ts, color=c, ls=':', lw=1.5)
    axC.set_xlim(0, T)
    axC.set_xlabel(r'$t$', fontsize=22)


def _mean_affinity_panel(cfg, d, pop, ax, lax):
    """Mean affinity z = e^-DG of the activated clones (doubled since t = 0, as L(t)) over time, weighted per
    clone (L) and per cell (Y; this one is Z/Y), with
    the activation front (z of the weakest activated clone). The number of activated clones L(t) on the right axis (lax)."""
    t, DG, w = d['t'], d['DG'], d['weights']
    N_B = d[f'N_B_clones_{pop}']
    act = (N_B[:, :1] > 0) & (N_B > 2.0 * N_B[:, :1])
    weights = {'L': w[:, None] * act, 'Y': w[:, None] * N_B * act}
    z = np.exp(-DG)[:, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        for q, ls in (('Y', '-'), ('L', '--')):
            ax.plot(t, np.sum(z * weights[q], axis=0) / np.sum(weights[q], axis=0), color=COLOR_Z, lw=3, ls=ls,
                    label=r'$\langle z \rangle_%s$' % q)
    front = np.where(act.any(axis=0), np.max(np.where(act, DG[:, None], -np.inf), axis=0), np.nan)
    ax.plot(t, np.exp(-front), color='grey', lw=1.5, ls=':', label=r'$\mathrm{front}$')
    ax.set_yscale('log')
    L = d[f'L_act_{pop}']
    lax.plot(t, np.where(L > 0, L, np.nan), color=COLOR_L, lw=3, label=r'$L$')
    lax.set_yscale('log')
    ax.set_xlim(0, cfg['T'])
    ax.set_ylim(top=1.5)
    ax.set_xlabel(r'$t$', fontsize=22)


def fig_energy_dynamics(cfg, dynamics, out):
    """Density in DG of the cells in activated clones (doubled since t = 0, as L(t)) at several times (A primary,
    B recall) and the number of activated clones L(t) and their cells Y(t) (left axis) with the potency Z(t) (right axis)
    (C primary, D recall), and the mean affinity <z> = <e^-DG> of the activated clones per cell (= Z/Y, from the
    clone sizes with Y's doubling rule) and per clone, with L(t) (E primary, F recall; _mean_affinity_panel).
    Each response has its own snapshot times (SNAPSHOT_TIMES), shown in A-B.
    Saved as energy_dynamics.pdf and, one panel per file with the same axis limits, energy_dynamics/<A-F>.pdf."""
    RESP = [('naive', 'naive', 'primary'), ('recall', 'memory', 'recall')]
    fig, axs = plt.subplots(3, 2, figsize=(14, 12.5), sharey='row',
                            gridspec_kw=dict(left=.1, right=.9, bottom=.07, top=.95, hspace=.35, wspace=.08))
    zaxs, laxs = [], []
    for k, (resp, pop, title) in enumerate(RESP):
        zax, lax = axs[1, k].twinx(), axs[2, k].twinx()
        zaxs.append(zax)
        laxs.append(lax)
        _energy_panels(cfg, dynamics[((), resp)], resp, pop, axs[0, k], axs[1, k], zax)
        _mean_affinity_panel(cfg, dynamics[((), resp)], pop, axs[2, k], lax)
        zax.tick_params(axis='y', which='both', labelsize=18, colors=COLOR_Z, labelright=(k == 1))
        lax.tick_params(axis='y', which='both', labelsize=18, colors=COLOR_L, labelright=(k == 1))
        axs[0, k].set_title(r'$\mathrm{%s}$' % title, fontsize=24)
    axs[0, 0].set_ylabel(r'$\mathrm{cells\ per\ unit}\ \Delta G$', fontsize=22)
    zlim = [f(z.get_ylim()[j] for z in zaxs) for j, f in enumerate((min, max))]
    zaxs[1].sharey(zaxs[0])
    zaxs[0].set_ylim(*zlim)
    llim = [f(a.get_ylim()[j] for a in laxs) for j, f in enumerate((min, max))]
    laxs[1].sharey(laxs[0])
    laxs[0].set_ylim(*llim)
    laxs[1].set_ylabel(r'$L(t)$', fontsize=22, color=COLOR_L)
    axs[1, 0].set_ylabel(r'$L(t)$', fontsize=22, color=COLOR_L)
    zaxs[1].set_ylabel(r'$Z(t)$', fontsize=22, color=COLOR_Z)
    axs[2, 0].set_ylabel(r'$\langle z \rangle = \langle e^{-\Delta G} \rangle$', fontsize=22)
    axs[2, 0].legend(axs[2, 0].lines + laxs[0].lines, [l.get_label() for l in axs[2, 0].lines + laxs[0].lines],
                     fontsize=16, loc='upper left')
    for ax in axs[0]:
        ax.legend(fontsize=14, loc='lower right')
    axs[1, 0].legend(axs[1, 0].lines[:1] + zaxs[0].lines[:1], [r'$L$', r'$Z$'], fontsize=16, loc='upper left')
    for ax, letter in zip(axs.flat, 'ABCDEF'):
        ax.tick_params(axis='both', labelsize=18)
        ax.text(-0.02, 1.04, r'$\mathbf{%s}$' % letter, transform=ax.transAxes, fontsize=26, ha='right')
    for ax in list(axs.flat) + zaxs + laxs:
        ax.set_box_aspect(1 / GOLDEN)
    lims = {(r, k): (axs[r, k].get_xlim(), axs[r, k].get_ylim()) for r in range(3) for k in range(2)}
    save(fig, out, 'energy_dynamics')

    ylabels = [r'$\mathrm{cells\ per\ unit}\ \Delta G$', r'$L(t)$', r'$\langle z \rangle = \langle e^{-\Delta G} \rangle$']
    for k, (resp, pop, title) in enumerate(RESP):
        figs = [plt.subplots(figsize=(7.5, 5.2), gridspec_kw=dict(left=.15, right=.85, bottom=.15, top=.9)) for _ in range(3)]
        (fA, aA), (fC, aC), (fE, aE) = figs
        zC, lE = aC.twinx(), aE.twinx()
        _energy_panels(cfg, dynamics[((), resp)], resp, pop, aA, aC, zC)
        _mean_affinity_panel(cfg, dynamics[((), resp)], pop, aE, lE)
        aE.legend(aE.lines + lE.lines, [l.get_label() for l in aE.lines + lE.lines], fontsize=14, loc='upper left')
        lE.set_ylim(*llim)
        lE.set_ylabel(r'$L(t)$', fontsize=22, color=COLOR_L)
        lE.tick_params(axis='y', which='both', labelsize=18, colors=COLOR_L)
        zC.set_ylim(*zlim)
        zC.set_ylabel(r'$Z(t)$', fontsize=22, color=COLOR_Z)
        zC.tick_params(axis='y', which='both', labelsize=18, colors=COLOR_Z)
        aA.legend(fontsize=14, loc='lower right')
        aC.legend(aC.lines[:1] + zC.lines[:1], [r'$L$', r'$Z$'], fontsize=16, loc='upper left')
        for a in (aA, aC, aE, zC, lE):
            a.set_box_aspect(1 / GOLDEN)
        for r, ((f, ax), ylab) in enumerate(zip(figs, ylabels)):
            ax.set_xlim(*lims[(r, k)][0])
            ax.set_ylim(*lims[(r, k)][1])
            ax.set_ylabel(ylab, fontsize=22, color=COLOR_L if r == 1 else 'k')
            ax.set_title(r'$\mathrm{%s}$' % title, fontsize=24)
            ax.tick_params(axis='both', labelsize=18)
            save(f, os.path.join(out, 'energy_dynamics'), 'ABCDEF'[2 * r + k])


def manuscript_figures(cfg, df, dynamics, out):
    """The legacy figure set that matches this run's configuration (nothing if none applies)."""
    rec, sweep = cfg['recall'], list(cfg['sweep'])
    recall = 'recall' in cfg['responses']
    if (recall and rec['antigen'] == 'feedback' and not rec['with_naive']
            and cfg['repertoire'] == 'polyclonal' and sweep == ['h0']):
        fig_potency(cfg, df, out)
    if recall and rec['with_naive'] and set(sweep) == {'h0', 'n_mem'} and 'N_B_clones' in cfg['dynamics']:
        fig_boost(cfg, df, dynamics, out)
    if cfg['repertoire'] == 'two_clone' and len(sweep) == 1 and 'N_B_clones' in cfg['dynamics']:
        fig_two_clone(cfg, df, dynamics, out)
    if (not sweep and cfg['repertoire'] == 'polyclonal'
            and {'N_B_clones', 'pi_clones'} <= set(cfg['dynamics'])):
        fig_ep(cfg, dynamics, out)
        if ((), 'recall') in dynamics and 'L_act' in cfg['dynamics']:
            fig_energy_dynamics(cfg, dynamics, out)


def main(run):
    cfg, df, dynamics = load(run)
    mfp.use_style(usetex=True)
    out = os.path.join(FIG_ROOT, run)
    plot_dynamics(cfg, dynamics, out)
    if len(cfg['sweep']) == 1:
        plot_observables_1d(cfg, df, out)
    elif len(cfg['sweep']) == 2:
        plot_observables_2d(cfg, df, out)
    plot_DG_dist(cfg, dynamics, out)
    plot_zipf(cfg, dynamics, out)
    manuscript_figures(cfg, df, dynamics, os.path.join(out, 'manuscript'))
    print(f"figures in {out}")


if __name__ == '__main__':
    main(sys.argv[1])
