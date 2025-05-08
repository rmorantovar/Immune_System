import sys
sys.path.append('../../lib/')
from funcs import*
from scipy.stats import gamma, expon
plt.rcParams['text.usetex'] = True
project = 'memory_response'
subproject = 'plots'
my_colors = [my_cyan, my_green, my_brown, my_red]

fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
output_plot = '../../../Figures/'+project+'/'+subproject
os.makedirs(output_plot, exist_ok=True)

t = np.linspace(0, 5, 101)
g = 5

ax.plot(t, expon.pdf(t, scale = 1), color = 'grey', lw = 2, label = r'$1$', ls = '--')
mean_expon = expon.stats(scale = 1, moments='m')
ax.vlines(mean_expon, 0, 1, color = 'k', ls = ':')

for i_g, g in enumerate([2, 3, 4, 5]):
	if g==2:
		ax.plot(t, gamma.pdf(t, a = g, scale = 1/g), lw = 2, color = my_colors[i_g], label = r'$%d$'%g)
	else:
		ax.plot(t, gamma.pdf(t, a = g, scale = 1/g), lw = 2, color = my_colors[i_g], label = r'$%d$'%g)
		mean_gamma = gamma.stats(a = g, scale = 1/g, moments='m')
		ax.vlines(mean_gamma, 0, 1, color = 'k', ls = ':')

my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
ax.set_xticks([])
# ax.set_yticks([])
ax.legend(fontsize = 18, title = r'$g$', title_fontsize = 20)
fig.savefig(output_plot + '/distributions_t.pdf')