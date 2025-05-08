import sys
sys.path.append('../../lib/')
from funcs import*
plt.rcParams['text.usetex'] = True
project = 'memory_response'
subproject = 'schematic'
output_plot = '../../../Figures/'+project+'/'+subproject
os.makedirs(output_plot, exist_ok=True)


#--------------------------------------------------------------------------------------------------------------

my_colors = [my_green_a, my_green_b, my_green_c]

for j in range(1, 3):
	fig, ax = plt.subplots(figsize=(9, 3), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})

	position = 4
	for i in range(3):
		circle = plt.Circle((position, 3.5), 1.5*(3/(i+2))**j, color=my_colors[i])
		ax.add_patch(circle)
		position = position + 1.5*(3/(i+2))**j + 1.5*(3/(i+3))**j + 2.5


	my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	ax.set_yticks([])
	ax.set_xticks([])
	ax.set_xlim(left = 0, right = 21)
	ax.set_ylim(bottom = 0, top = 7)
	# ax.legend(fontsize = 16)
	ax.axis('off')
	fig.savefig(output_plot + '/reweights%d.pdf'%j)

#--------------------------------------------------------------------------------------------------------------


fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})

tau = np.linspace(0, 10, 101)
titers = -tau + 10

ax.fill_between(tau[:31], -tau[:31] + 10, color = my_blue, alpha = .4)
ax.fill_between(tau[30:61], -tau[30:61] + 10, color = my_purple, alpha = .4)
ax.fill_between(tau[60:], -tau[60:] + 10, color = my_red, alpha = .4)

ax.plot(tau, titers, color = 'k', lw = 2, label = r'$\mathrm{Antiserum titers}$')

ax.hlines(7, 0, 10, color = 'darkgrey', lw = 1, ls = '--')
ax.hlines(4, 0, 10, color = 'k', lw = 1, ls = '--')

ax.vlines(3, 0, 7, color = 'darkgrey', lw = 1, ls = '--')
ax.vlines(6, 0, 4, color = 'k', lw = 1, ls = '--')


my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
ax.set_xticks([])
ax.set_yticks([7, 4], ['', ''])
ax.legend(fontsize = 16)
fig.savefig(output_plot + '/Titers_tau.pdf')

#--------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})

ax.plot(titers, (1+np.exp(titers - 2.5)**3)**(-1), color = my_red, lw = 2)
ax.hlines(1, 0, 10, color = 'k', lw = 1, ls = '--')

my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
ax.set_yticks([1], [r'$1$'])
ax.set_xticks([2.5], [''])
fig.savefig(output_plot + '/lambda_titers.pdf')

#--------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
lambda_ = 6.0   # growth rate
delta = 1     # decay rate
N0 = 1.0 
# Time span
t_span = (0, 10)
t_eval = np.linspace(*t_span, 500)

# Sigmoidal function: logistic
def pb(t, t0):
    # return 1 / (1 + np.exp(-k * (t - t0)))
    return (1+np.exp(-2*(t - t0)))**(-1)

# ODE builder with t0
def make_dNdt(t0):
    return lambda t, N: (lambda_ * (1 - pb(t, t0=t0)) - delta) * N

# Solve the ODE

for t0 in range(0, 7, 2):
	dNdt = make_dNdt(t0)
	sol = solve_ivp(dNdt, t_span, [N0], t_eval=t_eval)
	ax.plot(sol.t, sol.y[0], color = my_yellow2, lw = 5)

ax.hlines(1e8, 0, 10, color = 'k', lw = 1, ls = '--')

my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30)
ax.set_yticks([1e8], [r'$N_A^c$'])
ax.set_ylim(bottom = 0.8)
# ax.legend(fontsize = 16)
fig.savefig(output_plot + '/N_A_time.pdf')

#--------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})

ax.plot(tau, np.exp(np.log(1e8)/3*(tau - 3)), color = my_red, lw = 2)

ax.hlines(1e8, 0, 10, color = 'k', lw = 1, ls = '--')

my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
ax.set_yticks([1e8], [r'$N_A^c$'])
ax.set_xticks([])
ax.set_ylim(bottom = 0.8, top = 1e11)
# ax.legend(fontsize = 16)
fig.savefig(output_plot + '/N_A_peak_tau.pdf')



