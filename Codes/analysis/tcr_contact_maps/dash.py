import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
import itertools
import ast
import pyrepseq.plotting as rsp
import pyrepseq.distance as rsd
from collections import defaultdict
import scipy.cluster.hierarchy as hc
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#---------- import matrix ----------
E_model = 'TCRen'
#------- TCRen ------
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

#---------- import dataset ----------
epitope_gene = 'BMLF'

database = pd.read_csv('../../in/dash_human.csv')
database = pd.DataFrame(database)
database = database.loc[database['epitope'] == epitope_gene]

cg, linkage, cluster = rsp.similarity_clustermap(database[database['epitope']==epitope_gene],
                                                 alpha_column='cdr3_a_aa',
                                                 beta_column='cdr3_b_aa')

clus1 = database.loc[database['epitope'] == epitope_gene][cluster==6]
background = database.loc[database['epitope'] == epitope_gene][cluster!=6]

tcrs_clus_1 = []
for index, row in clus1.iterrows():
	seq_a = row['cdr3_a_aa']
	first_index_a = seq_a.find('D')
	clus1['cdr3_a_aa'][index] = row['cdr3_a_aa'][first_index_a:first_index_a+4]

	seq_b = row['cdr3_b_aa']
	first_index_b = seq_b.find('R')
	clus1['cdr3_b_aa'][index] = row['cdr3_b_aa'][first_index_b:first_index_b+6]

	tcrs_clus_1.append(clus1['cdr3_a_aa'][index] + clus1['cdr3_b_aa'][index])

tcrs_background = []
for index, row in background.iterrows():
	seq_a = row['cdr3_a_aa']
	# first_index_a = seq_a.find('D')
	background['cdr3_a_aa'][index] = row['cdr3_a_aa'][first_index_a:first_index_a+4]

	seq_b = row['cdr3_b_aa']
	# first_index_b = seq_b.find('R')
	background['cdr3_b_aa'][index] = row['cdr3_b_aa'][first_index_b:first_index_b+6]

	tcrs_background.append(background['cdr3_a_aa'][index] + background['cdr3_b_aa'][index])

# print(clus1[['cdr3_a_aa', 'cdr3_b_aa']])

#---------- import contact map info ----------
c_maps_df = pd.read_csv('../../in/contact_maps/contact_maps_PDB.csv')
miss_df = pd.read_csv('../../in/contact_maps/summary_PDB_structures.csv')

miss_df[['pdb.id', 'complex.species', 'peptide']].to_csv('../../out/epitopes_database.csv')
epitope_db = miss_df.loc[miss_df['pdb.id']=='3o4l']['peptide'][69]

cmap_a = c_maps_df.loc[(c_maps_df['pdb.id']=='3o4l') & (c_maps_df['chain.id.from']=='D') & (c_maps_df['region.type.from']=='CDR3')]
cmap_b = c_maps_df.loc[(c_maps_df['pdb.id']=='3o4l') & (c_maps_df['chain.id.from']=='E') & (c_maps_df['region.type.from']=='CDR3')]

l_t_a = len(cmap_a['residue.index.from'].unique())
l_t_b = len(cmap_b['residue.index.from'].unique())
l_t = l_t_a + l_t_b

cmap_a['residue.index.from'] -= np.min(cmap_a['residue.index.from'])
cmap_b['residue.index.from'] -= np.min(cmap_b['residue.index.from']) 
cmap_b['residue.index.from'] += l_t_a

cmap_merged = pd.concat([cmap_a, cmap_b])[['chain.type.from', 'residue.index.from', 'residue.aa.from', 'residue.index.to', 'residue.aa.to']].sort_values('residue.index.from')

l_p = 9
tcr_db = ''
indices = []
contact_map = np.zeros((l_t, l_p))
for index, row in cmap_merged.iterrows():
	contact_map[row['residue.index.from'], row['residue.index.to']] = 1
	if row['residue.index.from'] not in indices:
		tcr_db+=row['residue.aa.from']
	indices.append(row['residue.index.from'])
tcrs_clus_1.append(tcr_db)

print(tcrs_clus_1)

p = np.sum(contact_map)/(10*9)
ep_spec_e = []
background_e = []
random_e = []
random_cm_e = []

tcrs_random = [''.join(np.random.choice(Alphabet_T, 10)) for i in range(10000)]

pwm_df = potential_df[list(epitope_db)]

for tcr in tcrs_clus_1:
	E_matrix_df = pwm_df.reindex(list(tcr))
	EW_matrix = np.array(E_matrix_df)*contact_map
	E = np.sum(EW_matrix)
	ep_spec_e.append(E)

for tcr in tcrs_background:
	E_matrix_df = pwm_df.reindex(list(tcr))
	EW_matrix = np.array(E_matrix_df)*contact_map
	E = np.sum(EW_matrix)
	background_e.append(E)

for tcr in tcrs_random:
	E_matrix_df = pwm_df.reindex(list(tcr))
	EW_matrix = np.array(E_matrix_df)*contact_map
	E = np.sum(EW_matrix)
	random_e.append(E)

for i in range(10000):
	E_matrix_df = pwm_df.reindex(list(tcr_db))
	EW_matrix = np.array(E_matrix_df)*np.random.choice([1, 0], size=(l_t, l_p), p=[p, 1-p])
	E = np.sum(EW_matrix)
	random_cm_e.append(E)


E_matrix_df_db = pwm_df.reindex(list(tcr_db))
EW_matrix_db = np.array(E_matrix_df_db)*contact_map
E_db = np.sum(EW_matrix_db)

yscale = 'linear'

bins = np.linspace(-6, 9, 20)
plt.close()
plt.hist(ep_spec_e, bins = bins, alpha = .6, density = True, label = 'e-s_c1')
plt.hist(background_e, bins = bins, alpha = .6, density = True, label = 'e-s')
plt.vlines(E_db, 0, .5)
plt.legend()
plt.yscale(yscale)
plt.savefig('../../../Figures/repertoire_entropy/cluster_1_0.png')
plt.close()


plt.hist(ep_spec_e, bins = bins, alpha = .6, density = True, label = 'e-s_c1')
plt.hist(random_e, bins = bins, alpha = .6, density = True, label = 'r')
plt.vlines(E_db, 0, .5)
plt.legend()
plt.yscale(yscale)
plt.savefig('../../../Figures/repertoire_entropy/cluster_1_1.png')
plt.close()


plt.hist(ep_spec_e, bins = bins, alpha = .6, density = True, label = 'e-s_c1')
plt.hist(random_cm_e, bins = bins, alpha = .6, density = True, label = 'r_cm')
plt.vlines(E_db, 0, .5)
plt.legend()
plt.yscale(yscale)
plt.savefig('../../../Figures/repertoire_entropy/cluster_1_2.png')
plt.close()


plt.hist(background_e, bins = bins, alpha = .6, density = True, label = 'e-s')
plt.hist(random_e, bins = bins, alpha = .6, density = True, label = 'r')
plt.vlines(E_db, 0, .5)
plt.legend()
plt.yscale(yscale)
plt.savefig('../../../Figures/repertoire_entropy/cluster_1_3.png')
plt.close()

plt.hist(background_e, bins = bins, alpha = .6, density = True, label = 'e-s')
plt.hist(random_cm_e, bins = bins, alpha = .6, density = True, label = 'r_cm')
plt.vlines(E_db, 0, .5)
plt.legend()
plt.yscale(yscale)
plt.savefig('../../../Figures/repertoire_entropy/cluster_1_4.png')
plt.close()


plt.hist(random_e, bins = bins, alpha = .6, density = True, label = 'r')
plt.hist(random_cm_e, bins = bins, alpha = .6, density = True, label = 'r_cm')
plt.vlines(E_db, 0, .5)
plt.legend()
plt.yscale(yscale)
plt.savefig('../../../Figures/repertoire_entropy/cluster_1_5.png')
plt.close()




epitope_random = [''.join(np.random.choice(Alphabet_p, 9)) for i in range(10000)]

random_e_p = []
for epitope in epitope_random:
	pwm_df = potential_df[list(epitope)]
	E_matrix_df = pwm_df.reindex(list(tcr_db))
	EW_matrix = np.array(E_matrix_df)*contact_map
	E = np.sum(EW_matrix)
	random_e_p.append(E)

plt.hist(random_e_p, bins = 20, alpha = .6, density = True, label = 'random')
plt.vlines(E_db, 0, .2)
plt.legend()
plt.savefig('../../../Figures/repertoire_entropy/cluster_1_peptides.png')
