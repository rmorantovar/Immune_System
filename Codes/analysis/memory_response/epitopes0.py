import sys
sys.path.append('../../lib/')
from funcs import*

project = 'memory_response'
subproject = 'multi-epitope'
experiment = 1
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
energy_model = 'TCRen'

my_colors = [my_blue, my_green]
L0 = 1e8
output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
os.makedirs(output_plot, exist_ok=True)

# Set random seed for reproducibility
np.random.seed(42)
mu_xs = [-1, -1]
mu_ys = [-1, -4]

for j in range(2):
	mu_x = mu_xs[j]
	mu_y = mu_ys[j]

	K_stars = []

	figOmega, axOmega = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.25, 'right':.95, 'bottom':.25, 'top': 0.95})
	fig, ax = plt.subplots(figsize=(5, 5), gridspec_kw={'left':0.25, 'right':.95, 'bottom':.25, 'top': 0.95})
	# Parameters for the Gaussian distributions
	sigma = 1  # Same standard deviation

	# Generate 100 random samples from each distribution
	x_samples = -np.random.gumbel(loc = -mu_x, scale = sigma, size = 10000)
	y_samples = -np.random.gumbel(loc = -mu_y, scale = sigma, size = 10000)

	x = np.linspace(-7, 2, 10)

	# Create the scatter plot
	sns.kdeplot(ax = ax, x=x_samples, y=y_samples, cmap="Greys", fill=False, thresh=0, levels = 10)
	ax.scatter(x_samples, y_samples, alpha=0.1, edgecolors='k', facecolors = 'grey', marker = '.')
	ax.plot(x, x, ls = '--', color = 'k')

	my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 20, x_fontsize=30, y_fontsize=30)
	ax.set_xticks(range(-6, 2), ['', '', '', '', '', r"$\log K^{*(1)}$", '', ''])
	ax.set_yticks(range(-6, 2), ['']*(6-abs(mu_y)) + [r"$\log K^{*(2)}$"] + ['']*(abs(mu_y)+1))
	ax.set_xlim(left = -7, right = 2)
	ax.set_ylim(bottom = -7, top = 2)
	fig.savefig(output_plot + '/Ks_2d_' + str(abs(mu_x)) + '_' + str(abs(mu_y)) + '.pdf')

	antigen = 'ARTWMNLKPRTSW'
	l = len(antigen)
	antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
	#--------------------------Energy Motif--------------------------
	motif = get_motif(antigen_seq, energy_model, '../../')*1.2

	#Change values by the minimum
	E_m = -3 + 2*mu_y
	for i in np.arange(l):
		E_m+=np.min(motif[:,i], axis=0)
		motif[:,i]-=np.min(motif[:,i], axis=0)

	Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
	Kds = np.exp(Es[:-1])
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
	K_stars.append(E_r)
	
	axOmega.plot(Es[:-1], L0*Q0, color = my_green, lw = 3)

	#--------------------------Energy Motif--------------------------
	motif = get_motif(antigen_seq, energy_model, '../../')*1.2

	#Change values by the minimum
	E_m = -3 + 2*mu_x
	for i in np.arange(l):
		E_m+=np.min(motif[:,i], axis=0)
		motif[:,i]-=np.min(motif[:,i], axis=0)

	Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
	Kds = np.exp(Es[:-1])
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
	K_stars.append(E_r)
	
	axOmega.plot(Es[:-1], L0*Q0, color = my_blue, lw = 3)

	axOmega.hlines(1, -26, -4, ls = '--', color = 'k')

	my_plot_layout(ax = axOmega, xscale='linear', yscale= 'log', ticks_labelsize= 20, x_fontsize=30, y_fontsize=30)
	if j == 0:
		axOmega.set_xticks([K_stars[0]], [r"$\log K^{*(1)}=\log K^{*(2)}$"])
	else:
		axOmega.set_xticks(K_stars, [r"$\log K^{*(2)}$", r"$\log K^{*(1)}$"])
	axOmega.set_xlim(left = -26, right = -4)
	# axOmega.set_ylim(bottom = -7, top = 2)
	figOmega.savefig(output_plot + '/Omega_' + str(abs(mu_x)) + '_' + str(abs(mu_y)) + '.pdf')



