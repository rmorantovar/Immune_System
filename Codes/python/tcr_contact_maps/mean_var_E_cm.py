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
N_ens = 100
l = 4
ps = [.1, .2, .5, .8]#
tcrs = [list(tcr) for tcr in list(itertools.product(Alphabet_T, repeat=l))]
antigen_random = np.random.choice(Alphabet_p, l)
antigens = ['ARRR']#['TTMM', 'AAAR', 'WWWR', 'NNNN', 'RKDL', 'NCYH']

# data_mean = dict()
# data_var = dict()
# data_tcr_m = dict()
# data_cm_m = dict()
for antigen in tqdm(antigens):
	fig1, ax1 = plt.subplots(figsize = (10, 8))
	fig2, ax2 = plt.subplots(figsize = (10, 8))
	# data_mean[antigen] = dict()
	# data_var[antigen] = dict()
	# data_tcr_m[antigen] = dict()
	# data_cm_m[antigen] = dict()
	data_mean = dict()
	data_var = dict()
	data_tcr_m = dict()
	data_cm_m = dict()
	pwm_df = potential_df[list(antigen)]
	# Your code here
	for p in ps:
		# data_mean[antigen][p] = []
		# data_var[antigen][p] = []
		# data_tcr_m[antigen][p] = []
		# data_cm_m[antigen][p] = []
		data_mean[p] = []
		data_var[p] = []
		data_tcr_m[p] = []
		data_cm_m[p] = []
		for i_ens in range(N_ens):
			E_ref = 5
			Es = []
			for i_t, tcr in enumerate(tcrs[::1]):
				E_matrix_df = pwm_df.reindex(tcr)
				W_matrix = np.random.choice([1, 0], size=(l, l), p=[p, 1-p])
				EW_matrix = np.array(E_matrix_df)*W_matrix
				E = np.sum(EW_matrix)
				Es.append(E)
				if(E<E_ref):
					E_ref = E
					best_tcr = tcr
					best_cm = W_matrix
			E_mean = np.mean(Es)
			E_var = np.var(Es)
			# data_mean[antigen][p].append(E_mean)
			# data_var[antigen][p].append(E_var)
			# data_tcr_m[antigen][p].append(best_tcr)
			# data_cm_m[antigen][p].append(best_cm)
			data_mean[p].append(E_mean)
			data_var[p].append(E_var)
			data_tcr_m[p].append(best_tcr)
			data_cm_m[p].append(best_cm)

		# ax1.hist(data_mean[antigen][p], label = r'$%.2f$'%p)	
		# ax2.hist(data_var[antigen][p], label = r'$%.2f$'%p)	
		ax1.hist(data_mean[p], label = r'$%.2f$'%p)	
		ax2.hist(data_var[p], label = r'$%.2f$'%p)						
	my_plot_layout(ax = ax1, yscale = 'linear')
	ax1.legend(title = r'$p$', fontsize = 20, title_fontsize = 22)
	#ax1.set_xlim(left = -10, right = 10)
	figname = lambda i: '../../../Figures/tcr_contact_map/mean_E_'+ E_model + '_l-%d_antigen-'%l + antigen + f'_{i}.pdf'
	i=0
	while os.path.exists(figname(i)):
		i+=1
	fig1.savefig(figname(i))

	my_plot_layout(ax = ax2, yscale = 'linear')
	ax2.legend(title = r'$p$', fontsize = 20, title_fontsize = 22)
	#ax2.set_xlim(left = -10, right = 10)
	figname = lambda i: '../../../Figures/tcr_contact_map/var_E_'+ E_model + '_l-%d_antigen-'%l + antigen + f'_{i}.pdf'
	i=0
	while os.path.exists(figname(i)):
		i+=1
	fig2.savefig(figname(i))

	mean_df = pd.DataFrame(data_mean)
	var_df = pd.DataFrame(data_var)
	tcr_m_df = pd.DataFrame(data_tcr_m)
	cm_m_df = pd.DataFrame(data_cm_m)

	mean_df.to_csv('../../out/tcr_contact_maps/mean_df_'+antigen+'.csv')
	var_df.to_csv('../../out/tcr_contact_maps/var_df_'+antigen+'.csv')
	tcr_m_df.to_csv('../../out/tcr_contact_maps/tcr_m_df_'+antigen+'.csv')
	cm_m_df.to_csv('../../out/tcr_contact_maps/cm_m_df_'+antigen+'.csv')

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time/60:.3f} minutes")


