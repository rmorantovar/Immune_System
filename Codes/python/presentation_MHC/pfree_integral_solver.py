# pfree_integral_solver.py
# ------------------------
# Solve the self-consistent integral equation for the free-MHC fraction P_free
# in a continuum (mean-field) model of peptide-MHC competition.
#
# Equation solved:
#   1 - P_free = ∫ [A0(k) * k_er_plus * k_plus * P_free * M_tot] / [ (k_step + k) *
#                    (k_er_minus + M_tot * P_free * k_plus * (k_step/(k_step + k))**2 ) ] * Omega(k) * 1e7 dk
#
# where k == k_m_minus is the peptide-specific unbinding rate, Omega(k) is the pdf over k,
# and A0(k) is the peptide-specific source.
#
import sys
from turtle import color
sys.path.append('../../my_lib/')
from funcs import*
import argparse
import sys
from scipy.optimize import root_scalar
from scipy.integrate import simpson


HAS_MPL = True

def make_omega(E, dist='lognormal', **kwargs):
    """Build a normalized Omega(k) over the grid k.
    Supported dist:
      - 'lognormal' with kwargs: mu (default 0.0), sigma (default 1.0)
      - 'uniform'  with kwargs: kmin (default min(k)), kmax (default max(k))
      - 'gamma'    with kwargs: shape (alpha), scale (theta)
    Returns: Omega(k) array, normalized so integral Omega(k) dk = 1
    """
    if dist == 'lognormal':
        # mu = float(kwargs.get('mu', 0.0))
        # sigma = float(kwargs.get('sigma', 1.0))
        # k_safe = np.maximum(k, 1e-12)
        # pdf = (1.0 / (k_safe * sigma * np.sqrt(2*np.pi))) * np.exp(-(np.log(k_safe) - mu)**2 / (2*sigma**2))
        mu = float(kwargs.get('mu', 0.0))
        sigma = float(kwargs.get('sigma', 1.0))
        # k_safe = np.maximum(k, 1e-12)
        pdf = (1.0 / (sigma * np.sqrt(2*np.pi))) * np.exp(-(E - mu)**2 / (2*sigma**2))
    elif dist == 'uniform':
        kmin = float(kwargs.get('kmin', float(np.min(k))))
        kmax = float(kwargs.get('kmax', float(np.max(k))))
        width = max(kmax - kmin, 1e-12)
        pdf = np.where((k >= kmin) & (k <= kmax), 1.0/width, 0.0)
    elif dist == 'gamma':
        alpha = float(kwargs.get('shape', 2.0))
        theta = float(kwargs.get('scale', 1.0))
        from math import lgamma
        k_safe = np.maximum(k, 1e-12)
        logpdf = (alpha - 1.0) * np.log(k_safe) - (k_safe / theta) - lgamma(alpha) - alpha * np.log(theta)
        pdf = np.exp(logpdf)
    else:
        raise ValueError(f'Unsupported dist {dist!r}')
    norm = simpson(pdf, x=E)
    if norm <= 0 or not np.isfinite(norm):
        raise ValueError('Invalid distribution normalization (check parameters).')
    return pdf / norm


def make_A0(k, form='constant', **kwargs):
    """Build A0(k) on the grid k.
    Supported form:
      - 'constant': A0(k) = a0 (default 1.0)
      - 'powerlaw': A0(k) = c * k^p  (kwargs: c, p)
      - 'exp':      A0(k) = c * exp(-k / tau) (kwargs: c, tau)
    """
    if form == 'constant':
        a0 = float(kwargs.get('a0', 1.0))
        return np.full_like(k, a0, dtype=float)
    elif form == 'powerlaw':
        c = float(kwargs.get('c', 1.0))
        p = float(kwargs.get('p', 0.0))
        return c * np.power(k, p)
    elif form == 'exp':
        c = float(kwargs.get('c', 1.0))
        tau = float(kwargs.get('tau', 1.0))
        return c * np.exp(-k / tau)
    else:
        raise ValueError(f'Unsupported A0 form {form!r}')


def build_integrand(k, P_free, M_tot, k_plus, k_step, k_er_plus, k_er_minus, Omega, A0):
    """Compute integrand(k) for the RHS at a given P_free.
    integrand(k) = [A0(k) * k_er_plus * k_plus * P_free * M_tot] / [ (k_step + k) *
                    (k_er_minus + M_tot * P_free * k_plus * (k_step/(k_step + k))**2 ) ] * Omega(k) * 1e7
    """
    M_free = M_tot * P_free
    numerator = A0 * k_er_plus * k_plus * P_free# * M_tot
    denom_term = k_er_minus + M_tot * P_free * k_plus * (k_step / (k_step + k))**2
    integrand_vals = (numerator / ((k_step + k) * denom_term)) * Omega * 1e3
    return integrand_vals

def build_integrand2(k, P_free, M_tot, k_plus, k_step, k_er_plus, k_er_minus, Omega, A0):
    """Compute integrand(k) for the RHS at a given P_free.
    integrand(k) = [A0(k) * k_er_plus * k_plus * P_free * M_tot] / [ (k_step + k) *
                    (k_er_minus + M_tot * P_free * k_plus * (k_step/(k_step + k))**2 ) ] * Omega(k) * 1e3
    """
    M_free = M_tot * P_free
    numerator = A0 * k_er_plus * k_plus * P_free# * M_tot
    denom_term = k_er_minus + M_tot * k_plus * (k_step / (k_step + k))**2
    integrand_vals = (numerator / ((k_step + k) * denom_term)) * Omega * 1e3
    return integrand_vals


def solve_pfree(k_grid, Omega, A0, M_tot, k_plus, k_step, k_er_plus, k_er_minus, bracket=(1e-20, 1-1e-20)):
    """Solve 1 - P_free = integral integrand(k, P_free) dk using Brent's method on [bracket].
    Returns: (P_free, root_result)
    """
    def self_consistency(P_free):
        vals = build_integrand(k_grid, P_free, M_tot, k_plus, k_step, k_er_plus, k_er_minus, Omega, A0)
        rhs = 1 - simpson(vals, x=k_grid)
        return P_free - rhs

    res = root_scalar(self_consistency, bracket=bracket, method='brentq')
    if not res.converged:
        raise RuntimeError('Root finding did not converge.')
    return float(res.root), res

def solve_pfree2(k_grid, Omega, A0, M_tot, k_plus, k_step, k_er_plus, k_er_minus, bracket=(1e-20, 1-1e-20)):
    """Solve 1 - P_free = integral integrand(k, P_free) dk using Brent's method on [bracket].
    Returns: (P_free, root_result)
    """
    def self_consistency2(P_free):
        vals = build_integrand2(k_grid, P_free, M_tot, k_plus, k_step, k_er_plus, k_er_minus, Omega, A0)
        rhs = 1 - simpson(vals, x=k_grid)
        return P_free - rhs

    res = root_scalar(self_consistency2, bracket=bracket, method='brentq')
    if not res.converged:
        raise RuntimeError('Root finding did not converge.')
    return float(res.root), res


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description='Solve the integral equation for P_free.')
    # Model parameters
    ap.add_argument('--M_tot', type=float, default=10000.0, help='Total MHCs (default: 100.0)')
    ap.add_argument('--k_plus', type=float, default=1, help='MHC binding rate (default: 1.0)')
    ap.add_argument('--k_step', type=float, default=1e4, help='Stepping rate (default: 1.0)')
    ap.add_argument('--k_er_plus', type=float, default=1e6, help='ER entry rate (default: 1.0)')
    ap.add_argument('--k_er_minus', type=float, default=1e-1, help='ER exit rate (default: 0.1)')

    # k-grid
    ap.add_argument('--k_min', type=float, default=1e-2, help='min k_m_minus (default: 0.01)')
    ap.add_argument('--k_max', type=float, default=1e6, help='max k_m_minus (default: 10.0)')
    ap.add_argument('--num_points', type=int, default=10000, help='grid points (default: 1000)')

    # Distribution Omega(k)
    ap.add_argument('--dist', type=str, default='lognormal', choices=['lognormal', 'uniform', 'gamma'], help='distribution for Omega(k)')
    ap.add_argument('--mu', type=float, default=6.0, help='lognormal mu')
    ap.add_argument('--sigma', type=float, default=2.5, help='lognormal sigma')
    ap.add_argument('--kmin_omega', type=float, help='uniform lower bound (if set)')
    ap.add_argument('--kmax_omega', type=float, help='uniform upper bound (if set)')
    ap.add_argument('--gamma_shape', type=float, default=2.0, help='gamma shape (alpha)')
    ap.add_argument('--gamma_scale', type=float, default=1.0, help='gamma scale (theta)')

    # A0(k)
    ap.add_argument('--A0_form', type=str, default='constant', choices=['constant', 'powerlaw', 'exp'], help='functional form of A0(k)')
    ap.add_argument('--A0_a0', type=float, default=1e4, help='A0 constant value (for constant)')
    ap.add_argument('--A0_c', type=float, default=1.0, help='A0 coefficient (powerlaw/exp)')
    ap.add_argument('--A0_p', type=float, default=0.0, help='A0 power (for powerlaw)')
    ap.add_argument('--A0_tau', type=float, default=1.0, help='A0 decay constant (for exp)')

    # Plotting
    ap.add_argument('--plot', action='store_true', help='Plot integrand at the solved P_free')
    ap.add_argument('--save_plot', type=str, default='', help='If set, save plot to this path')

    return ap.parse_args(argv)

def main(argv=None):

    #----------------------------------------------------------------
    energy_model = 'TCRen'
    project = 'presentation_MHC'
    subproject = 'integral_solver'

    output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
    os.makedirs(output_plot, exist_ok=True)

    args = parse_args(argv)

    # Build k-grid
    k = np.logspace(np.log10(args.k_min), np.log10(args.k_max), args.num_points)

    # Build Omega(k)
    if args.dist == 'lognormal':
        Omega = make_omega(np.log(k), dist='lognormal', mu=args.mu, sigma=args.sigma)
    elif args.dist == 'uniform':
        kmin = args.kmin_omega if args.kmin_omega is not None else args.k_min
        kmax = args.kmax_omega if args.kmax_omega is not None else args.k_max
        Omega = make_omega(k, dist='uniform', kmin=kmin, kmax=kmax)
    elif args.dist == 'gamma':
        Omega = make_omega(k, dist='gamma', shape=args.gamma_shape, scale=args.gamma_scale)
    else:
        raise ValueError('Unsupported dist')
    
    fig, ax = plt.subplots(figsize=(8*1.62, 8), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    ax.plot(k, Omega, marker='.', label='Omega(k)', color='blue')

    my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    ax.set_xlabel('k')
    ax.set_ylabel('Omega(k)')
    ax.set_yscale('log')
    ax.set_xscale('log')

    fig.savefig(output_plot + '/Omega.pdf', dpi=150)
    print(np.sum(Omega[:-1]*np.diff(np.log(k))))
    
    # Build A0(k)
    if args.A0_form == 'constant':
        A0 = make_A0(k, form='constant', a0=args.A0_a0)
    elif args.A0_form == 'powerlaw':
        A0 = make_A0(k, form='powerlaw', c=args.A0_c, p=args.A0_p)
    elif args.A0_form == 'exp':
        A0 = make_A0(k, form='exp', c=args.A0_c, tau=args.A0_tau)
    else:
        raise ValueError('Unsupported A0_form')

    fig, ax = plt.subplots(figsize=(8*1.62, 8), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    M_tot_value = np.logspace(7, 11, num=150)

    cmap = plt.colormaps.get_cmap('turbo').resampled(5)  # Use 5 colors for 5 k_er_minus values

    for idx, k_er_minus in enumerate(tqdm(np.logspace(3, 4, num=5))):
        P_free_values = []
        P_free_values2 = []
        
        for M_tot in M_tot_value:

            # Solve for P_free
            P_free, res = solve_pfree(
                k_grid=np.log(k),
                Omega=Omega,
                A0=A0,
                M_tot=M_tot,
                k_plus=args.k_plus,
                k_step=args.k_step,
                k_er_plus=args.k_er_plus,
                k_er_minus=k_er_minus,
            )
            # Solve for P_free
            P_free2, res = solve_pfree2(
                k_grid=np.log(k),
                Omega=Omega,
                A0=A0,
                M_tot=M_tot,
                k_plus=args.k_plus,
                k_step=args.k_step,
                k_er_plus=args.k_er_plus,
                k_er_minus=k_er_minus,
            )
            P_free_values.append(P_free)
            P_free_values2.append(P_free2)

        color = cmap(idx)
        ax.plot(M_tot_value*1e-5, np.array(P_free_values), lw = 3, ls = '-', ms = 10, marker='', color = color, label=r'$%.1e$'%(k_er_minus))
        ax.plot(M_tot_value*1e-5, np.array(P_free_values2), lw = 3, ls = '--', ms = 10, marker='', color = ax.lines[-1].get_color())

    my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30)
    # ax.set_xlabel(r'$M_{\mathrm{tot}}$')
    # ax.set_ylabel('P_free')
    ax.legend(fontsize=24, loc=0, ncol=1, title=r'$k_{\mathrm{er}^-}/k_{\mathrm{er}^+}$', title_fontsize=24, frameon=False)
    fig.savefig(output_plot + '/P_free.pdf', dpi=150)



if __name__ == '__main__':
    main()