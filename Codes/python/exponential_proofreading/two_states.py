from email import parser
import sys
sys.path.append('../../my_lib/')
from funcs import*

energy_model = 'TCRen'
project = 'exponential_proofreading'
subproject = 'two_state'

output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
os.makedirs(output_plot, exist_ok=True)

lambdaA = 2.
lambdaB = 1.


def step_sample_from_breaks(A, B, n_points=10000):
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=int)  # 0/1 expected
    assert A.ndim == 1 and B.ndim == 1 and len(A) == len(B)
    assert np.all(np.diff(A) > 0), "A must be strictly increasing"

    # Sample from min(A) to max(A), inclusive
    t = np.linspace(A[0], 15, n_points)

    # Find the surrounding breakpoints for each t
    left  = np.searchsorted(A, t, side='left') - 1
    right = np.searchsorted(A, t, side='right')

    # Open intervals: (A[i], A[i+1]) -> take B[i]
    mask = (left >= 0) & (right < len(A)) & (left == right - 1)

    values = np.zeros_like(t, dtype=int)
    values[mask] = B[left[mask]]
    return t, values

for delta in [1e-1, 1, 1e1]:

    K = 1e-9
    fig, ax = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
    fig2, ax2 = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})

    times = [0.]
    state = [0.]
    # for i in range(0, 10):
    while times[-1] < 15:

        S = np.exp(lambdaA*15)/N_A
        alpha = S / (K + S)
        alpha = 1.
        rates = [alpha * (1 - state[-1]),  delta * (state[-1])]
        tot_rates = sum(rates)
        r = np.random.rand()
        t = np.log(1. / r) / tot_rates
        times.append(times[-1] + t)
        state.append(1.-state[-1])
        
    times = np.asarray(times)
    state = np.asarray(state)

    times_array = np.linspace(0, 10, 100000)
    times_array, growth_rate = step_sample_from_breaks(times, state)

    N = np.exp(lambdaB*np.cumsum(growth_rate[:-1]*np.diff(times_array)))

    ax.step(times[:], state[:], where='post')
    ax2.plot(times_array[:-1], N, label=r'$N(t)$')

    my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30)
    # ax.set_xticks([])
    ax.set_yticks([0, 1])
    # ax.set_ylim(bottom = 1e4, top = 2e13)
    ax.set_xlim(right = 15, left = 0)
    # ax.legend(fontsize = 30, loc = 2)
    fig.savefig(output_plot + '/state_delta_%.1e.pdf'%(delta))

    my_plot_layout(ax = ax2, xscale='linear', yscale= 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30)
    # ax.set_xticks([])
    ax.set_yticks([0, 1])
    # ax.set_ylim(bottom = 1e4, top = 2e13)
    ax2.set_xlim(right = 15, left = 0)
    # ax.legend(fontsize = 30, loc = 2)
    fig2.savefig(output_plot + '/N_delta_%.1e.pdf'%(delta))


