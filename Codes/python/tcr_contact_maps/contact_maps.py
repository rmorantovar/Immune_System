import sys
sys.path.append('../../lib/')
from functions_cm import*
from classes import*
plt.rcParams['text.usetex'] = True
import time
import itertools

def generate_sequences(l, d):
    return list(itertools.product(range(0, d), repeat=l))

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Immune_system/tcr_contact_maps/'

start_time = time.time()

# Your code here

l = 6  # Length of the sequences
d = 20  # Maximum integer value

#TCRs = generate_sequences(l, d)
TCRs = [list(seq) for seq in generate_sequences(l, d)]
#print(len(TCRs), (d**l))
E_models = ['TCRen']#, 'gaussian', MJ2]
ps = [0, .5, .2, .1]#
peptide = np.random.randint(0, d, l)


for E_model in E_models:
	fig, ax = plt.subplots(figsize = (10, 8))
	print('Energy model = ' + E_model)
	if E_model == 'TCRen':
		E_matrix = np.loadtxt("../../in/TCRen.txt")
	elif E_model == 'gaussian':
		E_matrix = np.random.normal(loc = 0, scale = 1, size = (d, d))
	motif = E_matrix[:, peptide]
	for p in ps:
		print('p=%.2f'%p)
		energies = []
		if p == 0:
			W_matrix = np.identity(l)
			for t, tcr in enumerate(tqdm(TCRs[::1])):
				E = energy(E_matrix, W_matrix, peptide, tcr, l)
				energies.append(E)
			E_m = np.sum(np.min(motif, axis = 0))
			E_m_temp = 0
			for i_l in range(l):
				E_m_temp += np.min(motif[:, i_l])
				motif[:, i_l]-=np.min(motif[:, i_l])
			Es, dE, Q0, betas = calculate_Q0(0.01, 50, 1000000, motif, E_m, l)
			ax.plot(Es, Q0*(1), ls = '--', color = 'black')
		else:
			for t, tcr in enumerate(tqdm(TCRs[::10])):
				W_matrix = np.random.choice([1, 0], size=(l, l), p=[p, 1-p])
				E = energy(E_matrix, W_matrix, peptide, tcr, l)
				energies.append(E)

		data_hist = np.histogram(energies, bins = np.linspace(-10, 20, 100), density = False)
		E_hist = data_hist[1][:-1]
		Omega_0 = np.diff(np.cumsum(data_hist[0])/len(energies))/np.diff(E_hist)
		ax.plot(E_hist[:-1], Omega_0, label = '%.2f'%p, ls = '', marker = '^')

	my_plot_layout(ax = ax, yscale = 'log')
	ax.legend(title = r'$p$', fontsize = 20, title_fontsize = 22)
	ax.set_xlim(left = -30, right = 30)
	fig.savefig('../../../Figures/tcr_contact_map/p_E_'+ E_model+'_l-%d.pdf'%l)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time/60:.3f} minutes")


