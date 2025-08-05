import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Import custom colors and plot layout if available
try:
    from funcs import my_plot_layout, my_blue2, my_purple, my_cyan, my_blue
except ImportError:
    def my_plot_layout(**kwargs): pass
    my_blue2 = 'blue'
    my_purple = 'purple'
    my_cyan = 'cyan'
    my_blue = 'deepskyblue'

def load_and_group_data(filepath, sheet_name, filter_col=None, filter_val=None, groupby_cols=None, header=1):
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=header)
    if filter_col and filter_val is not None:
        df = df[df[filter_col] == filter_val]
    grouped = df.groupby(groupby_cols).size().reset_index(name='count')
    return grouped

def bootstrap_clonal_ranks(data_grouped, mice, max_rank, n_ensemble, ax_step=None, color=None):
    zetas = []
    for rep in tqdm(range(n_ensemble)):
        mice_rep = mice if rep == n_ensemble - 1 else np.random.choice(mice, len(mice), replace=True)
        x_avg = np.zeros(max_rank)
        counts_per_ranking = np.zeros(max_rank)
        for mouse in mice_rep:
            counts = data_grouped[data_grouped['Mouse'] == mouse]['count'].to_numpy()
            largest = np.max(counts)
            x = np.flip(np.sort(counts))
            x = x[:max_rank] if len(x) > max_rank else np.pad(x, (0, max_rank - len(x)), mode='constant')
            for k in range(max_rank):
                if x[k] > 0:
                    counts_per_ranking[k] += 1
                    x_avg[k] += x[k] / largest
            if rep == n_ensemble - 1 and ax_step is not None:
                ax_step.step(range(1, len(x) + 1), x / largest, color=color, alpha=.5, lw=0.5)
        max_rank_eff = len(counts_per_ranking[counts_per_ranking > 2])
        x_avg = x_avg[:max_rank_eff] / counts_per_ranking[:max_rank_eff]
        params, _ = curve_fit(lambda x, m: m * x, np.log(range(1, max_rank_eff + 1))[:10], np.log(x_avg)[:10])
        zetas.append(-params[0])
    return zetas, x_avg, max_rank_eff

def fit_power_law(x, y, fit_range=10):
    params, _ = curve_fit(lambda x, m: m * x, np.log(x[:fit_range]), np.log(y[:fit_range]))
    return params[0]

def plot_rank_distribution(ax, x, y, color, marker, label):
    ax.plot(range(1, len(y) + 1), y, color=color, markerfacecolor="None", ms=12, alpha=1, ls='', marker=marker, label=label)

def plot_histogram(ax, zetas, bins, color, label):
    ax.hist(zetas, bins=bins, alpha=.7, label=label, color=color, density=True, histtype='stepfilled', edgecolor='k')

def save_figure(fig, path, transparent=True):
    fig.savefig(path, transparent=transparent)

def process_experiment(
    data_path, sheet_name, filter_col, filter_val, groupby_cols, 
    max_rank, n_ensemble, ax_r, ax_r2, ax_zeta, ax_zeta2, color, marker, label, output_prefix
):
    data_grouped = load_and_group_data(data_path, sheet_name, filter_col, filter_val, groupby_cols)
    mice = data_grouped[groupby_cols[0]].unique()
    zetas, x_avg, max_rank_eff = bootstrap_clonal_ranks(data_grouped, mice, max_rank, n_ensemble, ax_r, color)
    plot_rank_distribution(ax_r, range(1, max_rank_eff + 1), x_avg, color, marker, label % np.mean(zetas))
    plot_rank_distribution(ax_r2, range(1, max_rank_eff + 1), x_avg, color, marker, label % np.mean(zetas))
    ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0) * np.arange(1, max_rank_eff + 1) ** (-np.mean(zetas)), color=color, alpha=.8, lw=3)
    ax_r2.plot(np.arange(1, max_rank_eff + 1), np.exp(0) * np.arange(1, max_rank_eff + 1) ** (-np.mean(zetas)), color=color, alpha=.8, lw=3)
    plot_histogram(ax_zeta, zetas, bins=np.linspace(0.2, 1.6, 20), color=color, label=label % np.mean(zetas))
    plot_histogram(ax_zeta2, zetas, bins=np.linspace(0.2, 1.6, 20), color=color, label=label % np.mean(zetas))
    return zetas

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True

    project = 'memory_response'
    subproject = 'data'
    experiment = 0
    root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/mesin2020"
    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/' + project + '/' + subproject + '/' + str(experiment) + '/mesin2020'
    os.makedirs(output_plot, exist_ok=True)

    n_ensemble = 1000
    max_rank = 100

    color_vals = np.linspace(0, 2, 200)
    cmap = plt.get_cmap('managua_r')
    my_colors_alpha = [cmap(val) for val in color_vals]

    # Set up figures and axes
    fig_r, ax_r = plt.subplots(figsize=(8*1.62, 8), gridspec_kw={'left': 0.12, 'right': .95, 'bottom': .15, 'top': 0.94})
    fig_r2, ax_r2 = plt.subplots(figsize=(8, 8), gridspec_kw={'left': 0.15, 'right': .98, 'bottom': .15, 'top': 0.98})
    fig_zeta, ax_zeta = plt.subplots(figsize=(10*1.62, 8), gridspec_kw={'left': 0.12, 'right': .8, 'bottom': .15, 'top': 0.94})
    fig_zeta2, ax_zeta2 = plt.subplots(figsize=(10*1.62, 8), gridspec_kw={'left': 0.12, 'right': .8, 'bottom': .15, 'top': 0.94})

    # Example: Process Experiment 1 (Primary)
    zetas_primary = process_experiment(
        data_path=root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx",
        sheet_name='Photoactivation CGG',
        filter_col='Figure',
        filter_val=1,
        groupby_cols=['Mouse', 'Sequence'],
        max_rank=max_rank,
        n_ensemble=n_ensemble,
        ax_r=ax_r,
        ax_r2=ax_r2,
        ax_zeta=ax_zeta,
        ax_zeta2=ax_zeta2,
        color=my_blue2,
        marker='*',
        label=r'$%.2f$',
        output_prefix='primary'
    )

    # Example: Process Experiment 2 (Recall)
    zetas_recall = process_experiment(
        data_path=root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx",
        sheet_name='Fate-mapping CGG',
        filter_col='Figure',
        filter_val='4A',
        groupby_cols=['Mouse', 'Sequence'],
        max_rank=max_rank,
        n_ensemble=n_ensemble,
        ax_r=ax_r,
        ax_r2=ax_r2,
        ax_zeta=ax_zeta,
        ax_zeta2=ax_zeta2,
        color=my_purple,
        marker='o',
        label=r'$%.2f$',
        output_prefix='recall'
    )

    # Add more experiments as needed...

    # Layout and save figures
    my_plot_layout(ax=ax_r, yscale='log', xscale='log', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax_r.set_ylim(bottom=2e-2, top=1.1)
    ax_r.set_xlim(right=5e1)
    ax_r.legend(title=r'$\zeta$', fontsize=30, title_fontsize=30, loc=3)
    save_figure(fig_r, os.path.join(output_plot, 'ranking_B_cells.pdf'), transparent=.5)

    my_plot_layout(ax=ax_zeta, yscale='linear', xscale='linear', ticks_labelsize=40, x_fontsize=30, y_fontsize=30)
    ax_zeta.set_xlim(left=0.2, right=1.6)
    ax_zeta.legend(title=r'$\mathrm{sub-pop}$', fontsize=30, title_fontsize=30, loc=(1, 0))
    save_figure(fig_zeta, os.path.join(output_plot, 'zetas.pdf'), transparent=.5)