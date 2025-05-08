import sys
sys.path.append('../../lib/')
from funcs import*
from scipy.stats import entropy

# Function to generate x values from Gaussian distribution
def generate_x(g, sigma):
    return np.exp(np.random.normal(0, sigma, g))

# Function to compute y from x
def compute_y(x):
    return x / np.sum(x)

# Function to calculate Shannon entropy of y
def compute_entropy(y):
    y_nonzero = y[y > 0]  # Avoid log(0)
    return entropy(y_nonzero, base=None)

# Parameters
g_values = [5, 7, 9]  # Different values of g
sigma_values = np.logspace(-1.2, 2.3, 30)  # Range of sigma^2
sigma_values = np.logspace(-1.5, 2.0, 30)  # Range of sigma^2
n_samples = 1000  # Number of samples for averaging

# Storage for results
average_entropies = {g: [] for g in g_values}

# Main computation
for g in g_values:
    for sigma in sigma_values:
        entropies = []
        for _ in range(n_samples):
            x = generate_x(g, sigma)
            y = compute_y(x)
            entropies.append(compute_entropy(y))
        average_entropies[g].append(np.mean(entropies))

markers = ['o', '*', '^', 's', '.']
colors = [my_blue, my_green, my_red]
beta_values = np.sqrt(np.pi**2/(6*sigma_values))
cost_values = beta_values**2/5
cost_values = np.exp(beta_values)

# Plot results
fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
for i_g, g in enumerate(g_values):
    # plot = ax.plot(sigma_values, (average_entropies[g])/np.log(g), color = colors[i_g], ls = '', marker = markers[i_g], ms = 8, alpha = .6, label=fr'${g}$')
    plot = ax.plot(np.exp(average_entropies[g]), -cost_values, color = colors[i_g], ls = '', marker = markers[i_g], ms = 8, alpha = .8, label=fr'${g}$')
    # ax.plot(sigma_values, 1/(1+ sigma_values), color = plot[-1].get_color(), ls = '--')

# ax.plot(sigma_values[sigma_values>0.8], 1.4*sigma_values[sigma_values>0.8]**-1, color = 'grey', ls = '--')
# ax.plot(sigma_values[sigma_values<3], np.ones_like(sigma_values)[sigma_values<3], color = 'grey', ls = '--')

# ax.plot(sigma_values, 1/(1+ sigma_values + 1/(np.abs(np.log(np.sqrt(sigma_values)))+1)**2), color = 'k', ls = '-')
# ax.plot(sigma_values, 1/(1+ sigma_values), color = 'k', ls = '--', label = r'$\frac{1}{1+\textrm{Var} [\log{Z}]}$')
#
project = 'epitope_complexity'
subproject = 'entropy'
experiment = 1
output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
os.makedirs(output_plot, exist_ok=True)

my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=28, y_fontsize=28, bottom = -3., top = -.99, right = 6)
# ax.set_xlabel(r'$\textrm{Var} [\log{Z}]$')
# ax.set_ylabel(r'$\langle S(\textbf{z}) \rangle / \log{g}$')
# ax.set_xlabel(r'$g_r$')
# ax.set_ylabel(r'$\textrm{F}_c$')
ax.legend(fontsize = 18, title = r'$g$', title_fontsize = 20, loc = 0)
fig.savefig(output_plot+'/fitness_c.pdf')
# fig.savefig(output_plot+'/ID.pdf')
# 
