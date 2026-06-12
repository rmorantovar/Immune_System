import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import os
import sys
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
    f'Immune_System/_Repository/Figures/{project}/{subproject}/mesin2020/zipf/all'
)
os.makedirs(output_plot, exist_ok=True)

EXCEL_FILE   = root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx"
CLONE_COLS   = ['V', 'J', 'D']
# CLONE_COLS   = ['CDR3:']
MAX_RANK     = 100
MAX_RANK_FIT = 15
N_ENSEMBLE   = 400
RUN_EXPERIMENTS = {'1', '2', '3a', '3b', '2-3a', '2-3b', '4a', '4b'}
RUN_EXPERIMENTS = {'1', '2-3a','4a'}
# RUN_EXPERIMENTS = {'1', '2'}

color_vals      = np.linspace(0, 2, 200)
my_colors_alpha = [plt.get_cmap('autumn_r')(v) for v in color_vals]
my_colors2      = [my_purple, my_purple, my_purple, my_cyan, my_purple, my_blue2]


def model(x, m):
    return m * x

def barN(L, Z):
    return L ** Z / (1 - Z) * (L ** (1 - Z) - 1)

def S_th(L_e, Z):
    return (-Z * ((-L_e ** (1 - Z) * np.log(L_e)) / (L_e ** (1 - Z) - 1) + 1 / (1 - Z))
            + np.log((L_e ** (1 - Z) - 1) / (1 - Z)))


# ── Data ──────────────────────────────────────────────────────────────────────

def load_and_group(sheet_name, header, group_cols, filter_specs=None):
    """Load Excel sheet, apply filters, and return clone-count grouped DataFrame."""
    data = pd.read_excel(EXCEL_FILE, sheet_name=sheet_name, header=header)
    if filter_specs:
        for col, val in filter_specs.items():
            mask = data[col].isin(val) if isinstance(val, list) else data[col] == val
            data = data[mask]
    return data.groupby(group_cols).size().reset_index(name='count')


# ── Bootstrap ─────────────────────────────────────────────────────────────────

def run_bootstrap(data_grouped, mice, mouse_col,
                  ax_r=None, line_color=None, scaling_info=None):
    """Bootstrap Zipf exponent zeta over mice.

    On the last repetition (the non-bootstrap one), individual mouse curves are
    drawn on ax_r and per-mouse stats are appended to scaling_info['scaling_dict'].

    scaling_info keys: experiment, response, phenotype, scaling_dict

    Returns: x_avg_eff, max_rank_eff, zetas (array), zetas_mice (array)
    """
    zetas, zetas_mice = [], []

    for rep in tqdm(range(N_ENSEMBLE)):
        is_last            = rep == N_ENSEMBLE - 1
        mice_rep           = mice if is_last else np.random.choice(mice, len(mice), replace=True)
        x_avg              = np.zeros(MAX_RANK)
        counts_per_ranking = np.zeros(MAX_RANK)

        for mouse in mice_rep:
            counts = data_grouped[data_grouped[mouse_col] == mouse]['count'].to_numpy()
            if len(counts) == 0 or np.sum(counts) == 0:
                continue
            
            N       = np.sum(counts)
            S_i     = -np.sum((counts / N) * np.log(counts / N))
            largest = np.max(counts)
            x       = np.flip(counts[counts.argsort()])
            n_clones = len(x)
            pseudo_E = np.log(np.arange(1, len(counts) + 1)) / 2.5
            Z       = np.sum(x * np.exp(-pseudo_E))

            if is_last:
                if ax_r is not None:
                    ax_r.step(range(1, n_clones + 1), x / largest,
                              color=line_color, alpha=.5, lw=0.5)
                n_fit = min(n_clones, MAX_RANK_FIT)
                params_m, _ = curve_fit(model, np.log(range(1, n_fit + 1)), np.log((x / largest)[:n_fit]))
                zetas_mice.append(-params_m[0])
                if scaling_info is not None:
                    d = scaling_info['scaling_dict']
                    d['experiment'].append(scaling_info['experiment'])
                    d['response'].append(scaling_info['response'])
                    d['phenotype'].append(scaling_info['phenotype'])
                    d['N1'].append(largest)
                    d['L_act'].append(len(counts))
                    d['barN'].append(N)
                    d['S'].append(S_i)
                    d['zeta'].append(-params_m[0])
                    d['Z'].append(Z)
                    if -params_m[0] < 1:
                        d['SlogLact'].append(abs(- S_i + np.log(len(counts))))
                    else:
                        d['SlogLact'].append(abs(- S_i + np.log(len(counts))))

            x = (x[:MAX_RANK] if n_clones > MAX_RANK
                 else np.pad(x, (0, MAX_RANK - n_clones), mode='constant'))
            for k in range(MAX_RANK):
                if x[k] > 0:
                    counts_per_ranking[k] += 1
                    x_avg[k] += x[k] / largest

        max_rank_eff = len(counts_per_ranking[counts_per_ranking > 1])
        x_avg_eff    = x_avg[:max_rank_eff] / counts_per_ranking[:max_rank_eff]
        params, _    = curve_fit(
            model,
            np.log(range(1, max_rank_eff + 1))[:MAX_RANK_FIT],
            np.log(x_avg_eff[:])[:MAX_RANK_FIT],
        )
        zetas.append(-params[0])

    return x_avg_eff, max_rank_eff, np.array(zetas), np.array(zetas_mice)


# ── Plot helpers ──────────────────────────────────────────────────────────────

def recolor_last_lines(ax, n, color):
    for j in range(n):
        ax.lines[-(j + 1)].set_color(color)


def plot_ranking_result(ax_r, x_avg, max_rank_eff, zeta_mean, color, marker, label, ms=12):
    ax_r.plot(range(1, max_rank_eff + 1), x_avg,
              color=color, markerfacecolor='None', ms=ms, alpha=1, ls='', marker=marker, label=label)
    ax_r.plot(np.arange(1, max_rank_eff + 1), np.arange(1, max_rank_eff + 1) ** (-zeta_mean),
              color=color, alpha=.8, lw=3)


def plot_zeta_violin(ax_zeta, zetas, zetas_mice, position, color):
    parts = ax_zeta.violinplot([zetas], positions=[position], showmeans=True, showextrema=False)
    for body in parts['bodies']:
        body.set_facecolor(color)
        body.set_edgecolor('black')
        body.set_alpha(0.5)
    parts['cmeans'].set_color(color)
    ax_zeta.scatter(np.random.normal(position, 0.04, len(zetas_mice)), zetas_mice,
                    color=color, edgecolor='k', s=80, alpha=.5)
    ax_zeta.scatter(position, np.mean(zetas_mice), color=color, edgecolor='k', s=150)
    # ax_zeta.scatter(position, np.median(zetas_mice), color=color, edgecolor='k', s=150, marker='X')


def plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, zeta_mean, zeta_std, color):
    L   = np.linspace(1, 400, 100)
    L_e = np.linspace(1, 200, 100)
    z1 = np.linspace(0.01, 0.99, 100)
    z2 = np.linspace(1.01, 2.5, 100)
    Z   = zeta_mean
    Z_min = Z - 2*zeta_std
    Z_max = Z + 2*zeta_std
    
    S1_z = z1/(1-z1) + np.log(1-z1)
    S2_z = z2/(z2-1) - np.log(z2-1)
    
    ax_scaling1.plot(barN(L, Z), L ** Z, color=color, linestyle='--')
    ax_scaling1.fill_betweenx(L ** Z, barN(L, Z_min), barN(L, Z_max), color=color, alpha=0.1)
    
    ax_scaling2.plot(barN(L, Z), L,      color=color, linestyle='--')
    ax_scaling2.fill_betweenx(L, barN(L, Z_min), barN(L, Z_max), color=color, alpha=0.1)
    
    # ax_scaling3.plot(L ** Z, L,      color=color, linestyle='--')
    # ax_scaling3.fill_betweenx(L, L ** Z_min, L ** Z_max, color=color, alpha=0.1)
    ax_scaling3.plot(L ** Z, L,      color=color, linestyle='--')
    ax_scaling3.fill_betweenx(L, L, 200, color=my_blue, alpha=0.05)
    ax_scaling3.fill_betweenx(L, 0, L, color=my_red, alpha=0.05)
    
    
    ax_entropy.plot(L_e, S_th(L_e, Z),    color=color, linestyle='--')
    ax_entropy.fill_between(L_e, S_th(L_e, Z_min), S_th(L_e, Z_max), color=color, alpha=0.1)
    
    # ax_entropy2.plot(L_e, -S_th + np.log(L_e),    color=color, linestyle='--')
    # ax_entropy_logL_zeta.plot(z, abs(-S1_z),    color=color, linestyle='--')
    # ax_entropy_logL_zeta.plot(z, abs(-S2_z),    color=color, linestyle='--')


def apply_ranking_layout(fig_r, ax_r, suffix):
    ax_r.set_ylim(bottom=2e-2, top=1.1)
    ax_r.set_xlim(left = 0.9, right=5e1)
    # ax_r.set_xlabel(r'$\mathrm{Rank, } k $', fontsize=30)
    # ax_r.set_ylabel(r'$\mathrm{Relative clone size, } x_k $', fontsize=30)
    ax_r.tick_params(axis='y', labelsize=30)
    ax_r.tick_params(axis='x', labelsize=30)
    ax_r.set_yscale('log')
    ax_r.set_xscale('log')
    fig_r.savefig(output_plot + f'/ranking_B_cells_{suffix}.pdf', bbox_inches='tight', transparent=.5)


def apply_zeta_layout(fig_zeta, ax_zeta, suffix, tick_positions, tick_labels, x_labelsize=14):
    my_plot_layout(ax=ax_zeta, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax_zeta.set_xticks(tick_positions, tick_labels)
    ax_zeta.set_ylabel(r'$\zeta$', fontsize=30)
    ax_zeta.tick_params(axis='y', labelsize=30)
    ax_zeta.tick_params(axis='x', labelsize=x_labelsize, rotation=45)
    fig_zeta.savefig(output_plot + f'/zetas_{suffix}.pdf', bbox_inches='tight', transparent=.5)


# ── Figure setup ──────────────────────────────────────────────────────────────

fig_kw = dict(figsize=(8 * 1.62, 8), gridspec_kw={'left': .12, 'right': .95, 'bottom': .15, 'top': .94})
fig_r,        ax_r        = plt.subplots(**fig_kw)
fig_zeta,     ax_zeta     = plt.subplots(**fig_kw)
fig_scaling1, ax_scaling1 = plt.subplots(**fig_kw)
fig_scaling2, ax_scaling2 = plt.subplots(**fig_kw)
fig_scaling3, ax_scaling3 = plt.subplots(**fig_kw)
fig_scaling4, ax_scaling4 = plt.subplots(**fig_kw)
fig_entropy,  ax_entropy  = plt.subplots(**fig_kw)
fig_entropy_zeta,  ax_entropy_zeta  = plt.subplots(**fig_kw)
fig_entropy_logL_zeta,  ax_entropy_logL_zeta  = plt.subplots(**fig_kw)

scaling_dict = defaultdict(list)
violin_stats = []
apply_ranking_layout(fig_r, ax_r, '0')


# ── Experiment 1: naive GC (Figure 1D) ───────────────────────────────────────

if '1' in RUN_EXPERIMENTS:
    print("Experiment 1")
    data = load_and_group('Photoactivation CGG', header=1,
                          group_cols=['Mouse'] + CLONE_COLS,
                          filter_specs={'Figure': 1})
    mice = data['Mouse'].unique()

    x_avg, mre, zetas, zetas_mice = run_bootstrap(
        data, mice, 'Mouse', ax_r=ax_r, line_color=my_red,
        scaling_info=dict(experiment='1', response='naive', phenotype='GC', scaling_dict=scaling_dict),
    )
    print(np.mean(zetas), np.std(zetas))
    violin_stats.append(dict(experiment='1', response='naive', phenotype='GC',
                             violin_zeta_mean=np.mean(zetas), violin_zeta_std=np.std(zetas),
                             mice_zeta_mean=np.mean(zetas_mice), mice_zeta_std=np.std(zetas_mice)))
    recolor_last_lines(ax_r, len(mice), my_red)
    plot_ranking_result(ax_r, x_avg, mre, np.mean(zetas), my_red, '*',
                        r'$%.2f$' % np.mean(zetas) + ' ; GC', ms=18)
    plot_zeta_violin(ax_zeta, zetas, zetas_mice, position=0, color=my_red)
    plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, np.mean(zetas), np.std(zetas), my_red)
    apply_ranking_layout(fig_r, ax_r, '1')
    apply_zeta_layout(fig_zeta, ax_zeta, '1', [0, 1], ['', ''])


# ── Experiment 2: recall GC+fm (Figure 4A) ───────────────────────────────────

if '2' in RUN_EXPERIMENTS:
    print("Experiment 2 (Figure 4A)")
    data = load_and_group('Fate-mapping CGG', header=1,
                          group_cols=['Mouse'] + CLONE_COLS,
                          filter_specs={'Figure': '4A'})
    mice = data['Mouse'].unique()

    x_avg, mre, zetas, zetas_mice = run_bootstrap(
        data, mice, 'Mouse', ax_r=ax_r, line_color=my_blue,
        scaling_info=dict(experiment='2', response='recall', phenotype='GC + fm', scaling_dict=scaling_dict),
    )
    print(np.mean(zetas), np.std(zetas))
    violin_stats.append(dict(experiment='2', response='recall', phenotype='GC + fm',
                             violin_zeta_mean=np.mean(zetas), violin_zeta_std=np.std(zetas),
                             mice_zeta_mean=np.mean(zetas_mice), mice_zeta_std=np.std(zetas_mice)))
    recolor_last_lines(ax_r, len(mice), my_blue)
    plot_ranking_result(ax_r, x_avg, mre, np.mean(zetas), my_blue, 'o',
                        r'$%.2f$' % np.mean(zetas) + ' ; GC + fm')
    plot_zeta_violin(ax_zeta, zetas, zetas_mice, position=1, color=my_blue)
    # plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, np.mean(zetas), np.std(zetas), my_blue)
    apply_ranking_layout(fig_r, ax_r, '2')
    apply_zeta_layout(fig_zeta, ax_zeta, '2', [0, 1], ['', ''])


# ── Experiment 3a: recall per phenotype (Figure 4C) ──────────────────────────

if '3a' in RUN_EXPERIMENTS:
    print("Experiment 3 (Figure 4C)")
    data = load_and_group('Fate-mapping CGG', header=1,
                          group_cols=['Mouse', 'Phenotype'] + CLONE_COLS,
                          filter_specs={'Figure': '4C-H'})
    mice       = data['Mouse'].unique()
    phenotypes = data['Phenotype'].unique()
    colors_ph  = [my_blue, my_cyan, my_purple]
    print(phenotypes)

    for i_ph, ph in enumerate(phenotypes):
        x_avg, mre, zetas, zetas_mice = run_bootstrap(
            data[data['Phenotype'] == ph], mice, 'Mouse', ax_r=ax_r, line_color=colors_ph[i_ph],
            scaling_info=dict(experiment='3', response='recall', phenotype=ph, scaling_dict=scaling_dict),
        )
        violin_stats.append(dict(experiment='3', response='recall', phenotype=ph,
                                 violin_zeta_mean=np.mean(zetas), violin_zeta_std=np.std(zetas),
                                 mice_zeta_mean=np.mean(zetas_mice), mice_zeta_std=np.std(zetas_mice)))
        recolor_last_lines(ax_r, len(mice), colors_ph[i_ph])
        plot_ranking_result(ax_r, x_avg, mre, np.mean(zetas), colors_ph[i_ph], '^',
                            r'$%.2f$' % np.mean(zetas) + ' ; ' + ph)
        plot_zeta_violin(ax_zeta, zetas, zetas_mice, position=i_ph + 2, color=colors_ph[i_ph])
        # plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, np.mean(zetas), np.std(zetas), colors_ph[i_ph])

    apply_ranking_layout(fig_r, ax_r, '3')
    apply_zeta_layout(fig_zeta, ax_zeta, '3', [0, 1, 2, 3], ['', '', 'GC + fm', 'PB + fm'])


# ── Experiment 3b: recall combined (Figure 4C) ───────────────────────────────

if '3b' in RUN_EXPERIMENTS:
    print("Experiment 3 (Figure 4C) - 2")
    data = load_and_group('Fate-mapping CGG', header=1,
                          group_cols=['Mouse'] + CLONE_COLS,
                          filter_specs={'Figure': '4C-H'})
    mice = data['Mouse'].unique()

    x_avg, mre, zetas, zetas_mice = run_bootstrap(
        data, mice, 'Mouse', ax_r=ax_r, line_color=my_colors2[0],
        scaling_info=dict(experiment='3', response='recall', phenotype='combined', scaling_dict=scaling_dict),
    )
    print(np.mean(zetas), np.std(zetas))
    violin_stats.append(dict(experiment='3', response='recall', phenotype='combined',
                             violin_zeta_mean=np.mean(zetas), violin_zeta_std=np.std(zetas),
                             mice_zeta_mean=np.mean(zetas_mice), mice_zeta_std=np.std(zetas_mice)))
    recolor_last_lines(ax_r, len(mice), my_colors_alpha[int(np.mean(zetas) * 100)])
    # ranking point omitted (combined lines already visible on ax_r)
    plot_zeta_violin(ax_zeta, zetas, zetas_mice, position=4, color=my_blue2)
    # plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, np.mean(zetas), np.std(zetas), my_blue2)
    apply_ranking_layout(fig_r, ax_r, '3b')
    apply_zeta_layout(fig_zeta, ax_zeta, '3b',
                      [0, 1, 2, 3, 4], ['', '', 'GC + fm', 'GC + fm + recall', 'combined'])


# ── Experiments 2-3: shared data load ────────────────────────────────────────

if '2-3a' in RUN_EXPERIMENTS or '2-3b' in RUN_EXPERIMENTS:
    data_2_gcfm = load_and_group('Fate-mapping CGG', header=1,
                                  group_cols=['Mouse'] + CLONE_COLS,
                                  filter_specs={'Figure': '4A'})
    data_2_gcfm['Mouse'] = 'exp2_' + data_2_gcfm['Mouse'].astype(str)


# ── Experiment 2-3a: pooled recall, GC+fm only (Exp2 + Exp3) ─────────────────

if '2-3a' in RUN_EXPERIMENTS:
    print("Experiment 2-3a (pooled GC+fm from Exp2 and Exp3)")
    data_3_gcfm = load_and_group('Fate-mapping CGG', header=1,
                                  group_cols=['Mouse'] + CLONE_COLS,
                                  filter_specs={'Figure': '4C-H', 'Phenotype': 'GC + fm'})
    data_3_gcfm['Mouse'] = 'exp3_' + data_3_gcfm['Mouse'].astype(str)
    data_23a  = pd.concat([data_2_gcfm, data_3_gcfm], ignore_index=True)
    mice_23a  = data_23a['Mouse'].unique()

    x_avg, mre, zetas, zetas_mice = run_bootstrap(
        data_23a, mice_23a, 'Mouse', ax_r=ax_r, line_color=my_blue,
        scaling_info=dict(experiment='2-3a', response='recall', phenotype='GC + fm', scaling_dict=scaling_dict),
    )
    print(np.mean(zetas), np.std(zetas))
    violin_stats.append(dict(experiment='2-3a', response='recall', phenotype='GC + fm',
                             violin_zeta_mean=np.mean(zetas), violin_zeta_std=np.std(zetas),
                             mice_zeta_mean=np.mean(zetas_mice), mice_zeta_std=np.std(zetas_mice)))
    recolor_last_lines(ax_r, len(mice_23a), my_blue)
    plot_ranking_result(ax_r, x_avg, mre, np.mean(zetas), my_blue, 's',
                        r'$%.2f$' % np.mean(zetas) + ' ; GC + fm (pooled 2\&3)')
    plot_zeta_violin(ax_zeta, zetas, zetas_mice, position=5, color=my_blue)
    plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, np.mean(zetas), np.std(zetas), my_blue)
    apply_ranking_layout(fig_r, ax_r, '23a')
    apply_zeta_layout(fig_zeta, ax_zeta, '23a',
                      [0, 1, 2, 3, 4, 5],
                      ['', '', 'GC + fm', 'PB + fm', 'combined', 'GC+fm pool'])


# ── Experiment 2-3b: pooled recall, all phenotypes (Exp2 + Exp3) ─────────────

if '2-3b' in RUN_EXPERIMENTS:
    print("Experiment 2-3b (pooled all phenotypes from Exp2 and Exp3)")
    data_3_all = load_and_group('Fate-mapping CGG', header=1,
                                 group_cols=['Mouse'] + CLONE_COLS,
                                 filter_specs={'Figure': '4C-H'})
    data_3_all['Mouse'] = 'exp3_' + data_3_all['Mouse'].astype(str)
    data_23b   = pd.concat([data_2_gcfm, data_3_all], ignore_index=True)
    mice_23b   = data_23b['Mouse'].unique()

    x_avg, mre, zetas, zetas_mice = run_bootstrap(
        data_23b, mice_23b, 'Mouse', ax_r=ax_r, line_color=my_blue2,
        scaling_info=dict(experiment='2-3b', response='recall', phenotype='all', scaling_dict=scaling_dict),
    )
    print(np.mean(zetas), np.std(zetas))
    violin_stats.append(dict(experiment='2-3b', response='recall', phenotype='all',
                             violin_zeta_mean=np.mean(zetas), violin_zeta_std=np.std(zetas),
                             mice_zeta_mean=np.mean(zetas_mice), mice_zeta_std=np.std(zetas_mice)))
    recolor_last_lines(ax_r, len(mice_23b), my_blue2)
    plot_ranking_result(ax_r, x_avg, mre, np.mean(zetas), my_blue2, 'P',
                        r'$%.2f$' % np.mean(zetas) + ' ; all (pooled 2\&3)')
    plot_zeta_violin(ax_zeta, zetas, zetas_mice, position=6, color=my_blue2)
    # plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, np.mean(zetas), np.std(zetas), my_blue2)
    apply_ranking_layout(fig_r, ax_r, '23b')
    apply_zeta_layout(fig_zeta, ax_zeta, '23b',
                      [0, 1, 2, 3, 4, 5, 6],
                      ['', '', 'GC + fm', 'PB + fm', 'combined', 'GC+fm pool', 'all pool'])


# ── Experiment 4a: infection per sort (Figure 5) ─────────────────────────────

if '4a' in RUN_EXPERIMENTS:
    print("Experiment 4 (Figure 5)")
    colors_ph = [my_red, my_blue, my_cyan, my_purple]
    data = load_and_group('Influenza_IGH', header=2,
                          group_cols=['Experiment / Mouse', 'Sort'] + CLONE_COLS,
                          filter_specs={'Sort': ['GC', 'GC + fm', 'PB + fm', 'M']})
    mice       = data['Experiment / Mouse'].unique()
    phenotypes = data['Sort'].unique()
    print(mice)
    print(phenotypes)

    for i_ph, ph in enumerate(phenotypes[[0, 1]]):#, 3, 2]]):
        x_avg, mre, zetas, zetas_mice = run_bootstrap(
            data[data['Sort'] == ph], mice, 'Experiment / Mouse', ax_r=ax_r, line_color=colors_ph[i_ph],
            scaling_info=dict(experiment='4', response='recall', phenotype=ph, scaling_dict=scaling_dict),
        )
        print(np.mean(zetas), np.std(zetas))
        violin_stats.append(dict(experiment='4', response='recall', phenotype=ph,
                                 violin_zeta_mean=np.mean(zetas), violin_zeta_std=np.std(zetas),
                                 mice_zeta_mean=np.mean(zetas_mice), mice_zeta_std=np.std(zetas_mice)))
        recolor_last_lines(ax_r, len(mice), colors_ph[i_ph])
        plot_ranking_result(ax_r, x_avg, mre, np.mean(zetas), colors_ph[i_ph], 'D',
                            r'$%.2f$' % np.mean(zetas) + ' ; ' + ph)
        plot_zeta_violin(ax_zeta, zetas, zetas_mice, position=i_ph + 7, color=colors_ph[i_ph])
        plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, np.mean(zetas), np.std(zetas), colors_ph[i_ph])

    apply_ranking_layout(fig_r, ax_r, '4')
    apply_zeta_layout(fig_zeta, ax_zeta, '4',
                      list(range(11)),
                      ['exp1/naive', 'exp2/recall', 'exp3/recall', 'exp3/recall', 'exp3/recall',
                       'GC+fm pool', 'all pool',
                       'exp4/naive', 'exp4/recall', 'exp4/recall', 'exp4/M'],
                      x_labelsize=24)


# ── Experiment 4b: infection combined fm (Figure 5) ──────────────────────────

if '4b' in RUN_EXPERIMENTS:
    print("Experiment 4 (Figure 5) - 2")
    colors_ph = [my_blue2, my_blue, my_blue, my_blue, my_red]
    data = load_and_group('Influenza_IGH', header=2,
                          group_cols=['Experiment / Mouse', 'Sort2'] + CLONE_COLS,
                          filter_specs={'Sort2': 'fm'})
    mice       = data['Experiment / Mouse'].unique()
    phenotypes = data['Sort2'].unique()
    print(mice)
    print(phenotypes)

    for i_ph, ph in enumerate(phenotypes):
        x_avg, mre, zetas, zetas_mice = run_bootstrap(
            data[data['Sort2'] == ph], mice, 'Experiment / Mouse', ax_r=ax_r, line_color=colors_ph[i_ph],
            scaling_info=dict(experiment='4', response='recall', phenotype='combined', scaling_dict=scaling_dict),
        )
        print(np.mean(zetas), np.std(zetas))
        violin_stats.append(dict(experiment='4', response='recall', phenotype='combined',
                                 violin_zeta_mean=np.mean(zetas), violin_zeta_std=np.std(zetas),
                                 mice_zeta_mean=np.mean(zetas_mice), mice_zeta_std=np.std(zetas_mice)))
        recolor_last_lines(ax_r, len(mice), colors_ph[i_ph])
        # ranking point omitted (combined lines already visible on ax_r)
        plot_zeta_violin(ax_zeta, zetas, zetas_mice, position=i_ph + 11, color=colors_ph[i_ph])
        # plot_theory_curves(ax_scaling1, ax_scaling2, ax_scaling3, ax_entropy, ax_entropy_zeta, ax_entropy_logL_zeta, np.mean(zetas), np.std(zetas), my_blue2)

    apply_ranking_layout(fig_r, ax_r, '4b')
    apply_zeta_layout(fig_zeta, ax_zeta, '4b',
                      list(range(12)),
                      ['exp1/naive', 'exp2/recall', 'exp3/recall', 'exp3/recall', 'exp3/recall',
                       'GC+fm pool', 'all pool',
                       'exp4/naive', 'exp4/recall', 'exp4/recall', 'exp4/M', 'exp4/recall combined'],
                      x_labelsize=24)


# ── Scaling & entropy scatter plots ───────────────────────────────────────────

scaling_results = pd.DataFrame(scaling_dict)
palette = [my_red, my_blue, my_cyan, my_purple]

sns.scatterplot(data=scaling_results, x='barN', y='N1',
                hue='phenotype', style='experiment', ax=ax_scaling1,
                s=150, palette=palette, edgecolors='black', alpha=0.8)
my_plot_layout(ax=ax_scaling1, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax_scaling1.set_xlabel(r'$N_B^{\mathrm{tot}}$', fontsize=30)
ax_scaling1.set_ylabel(r'$N_1$', fontsize=30)
ax_scaling1.set_ylim(bottom=0, top=80)
ax_scaling1.set_xlim(left=1, right=400)
ax_scaling1.tick_params(axis='both', labelsize=30)
ax_scaling1.legend(title='Response', title_fontsize=20, fontsize=15, loc=4)
fig_scaling1.savefig(output_plot + '/size_scaling_1_linear.pdf', bbox_inches='tight', transparent=.5)

sns.scatterplot(data=scaling_results, x='barN', y='L_act',
                hue='phenotype', style='experiment', ax=ax_scaling2,
                s=150, palette=palette, edgecolors='black', alpha=0.8)
my_plot_layout(ax=ax_scaling2, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax_scaling2.set_xlabel(r'$N_B^{\mathrm{tot}}$', fontsize=30)
ax_scaling2.set_ylabel(r'$L_{act}$', fontsize=30)
ax_scaling2.set_ylim(bottom=0, top=125)
ax_scaling2.set_xlim(left=1, right=200)
ax_scaling2.tick_params(axis='both', labelsize=30)
ax_scaling2.legend(title='Response', title_fontsize=20, fontsize=15, loc=2)
fig_scaling2.savefig(output_plot + '/size_scaling_2_linear.pdf', bbox_inches='tight', transparent=.5)

sns.scatterplot(data=scaling_results, x='N1', y='L_act',
                hue='phenotype', style='experiment', ax=ax_scaling3,
                s=150, palette=palette, edgecolors='black', alpha=0.8)
x = np.linspace(1, 240, 100)
ax_scaling3.plot(x, x, 'k--')  # Example line plot, replace with actual function if needed
my_plot_layout(ax=ax_scaling3, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
# ax_scaling3.set_xlabel(r'$N_1$', fontsize=30)
# ax_scaling3.set_ylabel(r'$L_{act}$', fontsize=30)
# ax_scaling3.set_xscale('log')
# ax_scaling3.set_yscale('log')
ax_scaling3.set_ylim(bottom=0, top=125)
ax_scaling3.set_xlim(left=0, right=60)
ax_scaling3.tick_params(axis='both', labelsize=30)
ax_scaling3.legend(title='Response', title_fontsize=20, fontsize=15, loc=1)
fig_scaling3.savefig(output_plot + '/size_scaling_3_linear.pdf', bbox_inches='tight', transparent=.5)

scaling_results['barN_prediction'] = scaling_results['N1']/(1-scaling_results['zeta'])*(scaling_results['L_act']**(1-scaling_results['zeta']) - 1)
sns.scatterplot(data=scaling_results, x='barN_prediction', y='barN',
                style='experiment', hue='phenotype', ax=ax_scaling4,
                s=150, palette=palette, edgecolors='black', alpha=0.8)
ax_scaling4.plot(x, x, 'k--')  # Example line plot, replace with actual function if needed
my_plot_layout(ax=ax_scaling4, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
# ax_scaling4.set_xlabel(r'$N_1$', fontsize=30)
# ax_scaling4.set_ylabel(r'$L_{act}$', fontsize=30)
# ax_scaling4.set_ylim(bottom=0, top=125)
# ax_scaling4.set_xlim(left=0, right=80)
ax_scaling4.tick_params(axis='both', labelsize=30)
ax_scaling4.legend(title='Response', title_fontsize=20, fontsize=15, loc=4)
fig_scaling4.savefig(output_plot + '/size_scaling_4.pdf', bbox_inches='tight', transparent=.5)

sns.scatterplot(data=scaling_results, x='L_act', y='S',
                hue='phenotype', style='experiment', ax=ax_entropy,
                s=150, palette=palette, edgecolors='black', alpha=0.8)
my_plot_layout(ax=ax_entropy, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax_entropy.set_xlabel(r'$L_{act}$', fontsize=30)
ax_entropy.set_ylabel(r'$S$', fontsize=30)
ax_entropy.set_xlim(left=2, right=160)
ax_entropy.set_ylim(bottom=0.5, top=4.8)
# ax_entropy.set_xscale('log')
ax_entropy.tick_params(axis='both', labelsize=30)
ax_entropy.legend(title='Response', title_fontsize=20, fontsize=15, loc=4)
fig_entropy.savefig(output_plot + '/size_scaling_entropy.pdf', bbox_inches='tight', transparent=.5)


sns.scatterplot(data=scaling_results, x='N1', y='Z',
                hue='phenotype', style='experiment', ax=ax_entropy_zeta,
                s=150, palette=palette, edgecolors='black', alpha=0.8)
ax_entropy_zeta.plot(x, x/(1.05+1/2.3-1), 'k--')  # Example line plot, replace with actual function if needed
my_plot_layout(ax=ax_entropy_zeta, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
# ax_entropy_zeta.set_xlabel(r'$\zeta$', fontsize=30)
# ax_entropy_zeta.set_ylabel(r'$S$', fontsize=30)
ax_entropy_zeta.set_xlim(left=1, right=60)
ax_entropy_zeta.set_ylim(bottom=10, top=130)
ax_entropy_zeta.tick_params(axis='both', labelsize=30)
ax_entropy_zeta.set_xscale('log')
ax_entropy_zeta.set_yscale('log')
ax_entropy_zeta.legend(title='Response', title_fontsize=20, fontsize=15, loc=4)
fig_entropy_zeta.savefig(output_plot + '/size_scaling_entropy_zeta.pdf', bbox_inches='tight', transparent=.5)

sns.scatterplot(data=scaling_results, x='L_act', y='SlogLact',
                hue='phenotype', style='experiment', ax=ax_entropy_logL_zeta,
                s=150, palette=palette, edgecolors='black', alpha=0.8)
my_plot_layout(ax=ax_entropy_logL_zeta, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
ax_entropy_logL_zeta.set_xlabel(r'$L_{act}$', fontsize=30)
ax_entropy_logL_zeta.set_ylabel(r'$|\log L_{act}-S|$', fontsize=30)
# ax_entropy_logL_zeta.set_xlim(left=2, right=160)
# ax_entropy_logL_zeta.set_ylim(bottom=-0.05, top=4)
ax_entropy_logL_zeta.tick_params(axis='both', labelsize=30)
ax_entropy_logL_zeta.legend(title='Response', title_fontsize=20, fontsize=15, loc=4)
fig_entropy_logL_zeta.savefig(output_plot + '/size_scaling_entropy-logL_zeta.pdf', bbox_inches='tight', transparent=.5)


# ── Save zeta statistics ──────────────────────────────────────────────────────

pd.DataFrame(violin_stats).to_csv(output_plot + '/zeta_stats.csv', index=False)