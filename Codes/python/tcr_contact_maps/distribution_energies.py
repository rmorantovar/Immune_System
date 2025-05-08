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

#---------- import matrix ----------
E_model = 'TCRen'
#------- TCRen ------
potential_csv = pd.read_csv('../../in/TCRen_potential.csv')
Alphabet_T = np.unique(potential_csv['residue.aa.from'])
Alphabet_p = np.unique(potential_csv['residue.aa.to'])
#print(Alphabet_T, Alphabet_p)

potential_csv = pd.read_csv('../../in/TCRen_potential.csv')
Alphabet_T = np.unique(potential_csv['residue.aa.from'])
Alphabet_p = np.unique(potential_csv['residue.aa.to'])
#print(Alphabet_T, Alphabet_p)

potential_dict = {}
for aa_p in Alphabet_p:
	#print(aa_p)
	potential_dict[aa_p] = {}
	A = potential_csv[potential_csv['residue.aa.to']==aa_p]
	for aa_T in Alphabet_T:
		#print(aa_T)
		B = A[A['residue.aa.from']==aa_T]
		potential_dict[aa_p][aa_T] = float(B['TCRen'])

potential_df = pd.DataFrame(potential_dict)

#---------- ---------- ----------
l = 4
ps = [.1, .2, .5, .8]

antigen_random = np.random.choice(Alphabet_p, l)
antigens = [['TTMM'], ['AAAR'], ['WWWR'], ['NNNN'], ['RKDL'], ['NCYH']]
for antigen in antigens:
	pwm_df = potential_df[list(antigen)]
	# Your code here
	fig, ax = plt.subplots(figsize = (10, 8))

	tcrs = [list(tcr) for tcr in list(itertools.product(Alphabet_T, repeat=l))]

	for p in ps:
		Es = []
		for i_t, tcr in enumerate(tqdm(tcrs)):
			E_matrix_df = pwm_df.reindex(tcr)
			W_matrix = np.random.choice([1, 0], size=(l, l), p=[p, 1-p])
			EW_matrix = np.array(E_matrix_df)*W_matrix
			E = np.sum(EW_matrix)
			Es.append(E)

		data_hist = np.histogram(Es, bins = 80, density = False)
		E_hist = data_hist[1][:-1]
		Omega_0 = np.diff(np.cumsum(data_hist[0])/len(Es))/np.diff(E_hist)
		ax.plot(E_hist[:-1], Omega_0, label = '%.2f'%p, ls = '', marker = '^', ms = 5)

	Es = []
	W_matrix = np.identity(l)
	for i_t, tcr in enumerate(tqdm(tcrs)):
		E_matrix_df = pwm_df.reindex(tcr)
		EW_matrix = np.array(E_matrix_df)*W_matrix
		E = np.sum(EW_matrix)
		Es.append(E)

	data_hist = np.histogram(Es, bins = 80, density = False)
	E_hist = data_hist[1][:-1]
	Omega_0 = np.diff(np.cumsum(data_hist[0])/len(Es))/np.diff(E_hist)
	ax.plot(E_hist[:-1], Omega_0, label = r'$\textrm{Identity}$', ls = '', marker = 'o', ms = 5)

	my_plot_layout(ax = ax, yscale = 'log')
	ax.legend(title = r'$p$', fontsize = 20, title_fontsize = 22)
	ax.set_xlim(left = -10, right = 10)
	figname = lambda i: '../../../Figures/tcr_contact_map/densities/Omega_0_E_'+ E_model + '_l-%d_antigen-'%l + antigen + f'_{i}.pdf'
	i=0
	while os.path.exists(figname(i)):
		i+=1
	fig.savefig(figname(i))

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time/60:.3f} minutes")


