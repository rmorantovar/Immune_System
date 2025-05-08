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

#---------- import datasets ----------
cm_db = pd.read_csv('../../in/contact_maps/contact_maps_PDB.csv')
#print(cm_db.columns)
cm_db = cm_db.loc[cm_db['region.type.from']=='CDR3'].reset_index(drop=True)
cm_db = cm_db[['pdb.id','residue.index.to',
       'residue.index.from', 'chain.type.from',
       'residue.aa.from',
       'residue.aa.to']]
cm_db = cm_db.rename(columns={'pdb.id':'pdb_id', 'residue.index.to': "index_p", 'residue.index.from': 'index_t', 'chain.type.from':'chain', 'residue.aa.from':'aa_t', 'residue.aa.to':'aa_p'})
miss_cm_db = pd.read_csv('../../in/contact_maps/summary_PDB_structures.csv')
dash_db = pd.read_csv('../../in/dash_human.csv')

organisms = ['EBV', 'CMV', 'Influenza A']
peptides = ['GLCTLVAML', 'NLVPMVATV', 'GILGFVFTL']
genes = ['BMLF', 'pp65', 'M1']
clusters = {'BMLF':[ 5,  5, 22, 34, 34, 35, 41, 12,  6, 16, 13,  6,  6,  6,  7,  6,  6,
        6,  6,  6,  6,  6,  6,  6,  6,  6,  4,  4,  6,  9,  6,  6, 11,  2,
       19, 19, 18, 18, 18, 14, 14,  8,  8,  8,  8, 17, 15,  1, 17, 17,  8,
        8,  8,  8, 24,  8, 21, 21,  6,  6, 32, 37, 28, 28, 28, 40,  6, 23,
       23, 38, 38, 33, 25, 25, 20, 39, 39,  3,  3,  3, 28, 28, 28, 36, 10,
       27,  6, 39, 31, 26, 29, 30], 'pp65' : [21, 13, 20, 29, 26, 50,  2,  2,  1,  1, 19,  8, 38, 43, 28, 31, 48,
       46, 46, 46, 22, 23, 49, 11, 42, 52, 15, 39, 39, 39, 39,  7, 40, 14,
       45, 30, 30, 44, 18, 16, 10, 32, 10, 25, 25, 36, 27, 24, 24, 24, 17,
       35, 35, 33, 34, 12, 12, 47, 32,  5, 41,  5,  6,  4,  3,  5,  9, 51,
       37], 'M1' : [ 41,  80,  57,  57,  57,  39,  57, 103,   1,  57,  57,  15,  34,
        61,  54,  83,  60,  84,  91,  53,  63,  51,  51,  51,  51,  57,
        52,  96,  30,  30,  30,  30,  30,  16,   7,   7,  61,  59,  59,
        54,  54,  59,  59,  61,  61,  61,  56,  61,  54,  54,  54,  59,
        59,  54,  59,  54,  59,  54,  57,  49,  57,  57,  57,  57,  61,
        61,  61,  57,  61,  57,  37,  57,  61,  61,  81,  61,  31,  54,
        61,  61,  61,  61,  61,  57,  61,  61,  61,  61,  61,  61,  59,
        57,  57,  57,  59,  57,  57,  58,  57,  57,  57,  57,  57,  57,
        57,  57,  57,  57,  57,  57,  57,  61,  61,  76,  57,  90,  85,
        54,  57,  57,  62,  64,  57,  89,  43,  32,  57,  46,  46,  54,
        55,  57,  11,  57,  57,  57,  57,  57,  57,  57,  57,  57,  57,
        54,  54,  54,  68,  67,  45,  40,  97,  94,  94,  94,  69,  86,
        78,  87, 101,  33,  79,   4,   4,   4,   6,  35,  22,  54,  21,
        24,  24,  25,   3,   3,  57,  59, 104,  54,  93,  57,  54,  44,
         2,  23,  57,  57,  57,  57,  71,  36,  57,  57,  47,  57,  57,
        26,  57,  88,  10,  70,  74,  57,  56,  98,  92,  18,  61,  57,
        57,  38,  12,  17, 102, 100,  72,  75,  20,  27,  28,  73,  42,
        14,  57,  57,  57,  99,  77,  57,  57,  56,  61,  61,  48,  30,
        95,  57,  57,  30,  30,  30,  65,  50,  54,  82, 105,   9,  66,
        29,  57,   8,   5, 106,  13,  19]}

dash_db['cluster'] = np.ones(len(dash_db))
larger_clusters = dict()
for gene in genes:
       dash_db.loc[dash_db['epitope']==gene, 'cluster'] = [str(int(i)) for i in clusters[gene]]
       cluster_sizes = np.unique(clusters[gene], return_counts = True)[1]
       larger_clusters_gene = [str(int(i)) for i in (np.argsort(cluster_sizes)+1)[-5:]]
       larger_clusters[gene] = larger_clusters_gene

ids_dict = defaultdict(list)

for index, row in miss_cm_db.iterrows():
	p = row['peptide']
	if p in peptides:
		i_p = peptides.index(p)
		ids_dict['peptide'].append(p)
		ids_dict['pdb_id'].append(row['pdb.id'])
		ids_dict['gene'].append(genes[i_p])

ids_df = pd.DataFrame(ids_dict)#.sort_values('peptide').reset_index(drop=True)
cm_db = cm_db.drop(cm_db[~cm_db.pdb_id.isin(ids_df['pdb_id'])].index)

l_p = 9


for i, row in ids_df.iterrows():
       cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRA'), 'index_t'] -= np.min(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRA'), 'index_t'])
       cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRB'), 'index_t'] -= np.min(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRB'), 'index_t'])
       c_t_a = len(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRA'), 'index_t'])
       c_t_b = len(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRB'), 'index_t'])
       c_t = c_t_a + c_t_b
       i_max_t_a = np.max(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRA'), 'index_t'].unique())
       i_max_t_b = np.max(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRB'), 'index_t'].unique())
       cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRB'), 'index_t'] += (i_max_t_a + 1)
       i_max_t = np.max(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRB'), 'index_t'].unique()) + 1
       l_t_a = len(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRA'), 'index_t'].unique())
       l_t_b = len(cm_db.loc[(cm_db['pdb_id'] == row['pdb_id']) & (cm_db['chain']=='TRB'), 'index_t'].unique())
       l_t = l_t_a + l_t_b
       tcr_a = ''
       tcr_b = ''
       tcr = ''
       indices_a = []
       indices_b = []
       indices = []

       contact_map_completed = np.zeros((i_max_t, l_p))
       for j, row2 in cm_db.loc[cm_db['pdb_id'] == row['pdb_id']].sort_values(['chain', 'index_t']).iterrows():
              contact_map_completed[row2['index_t'], row2['index_p']] = 1
              if row2['index_t'] not in indices:
                     if row2['chain'] == 'TRA':
                            tcr_a += row2['aa_t']
                            indices_a.append(row2['index_t'])
                     else:
                            tcr_b += row2['aa_t']
                            indices_b.append(row2['index_t'])
                     tcr = tcr_a  + tcr_b
              indices.append(row2['index_t'])
       contact_map = contact_map_completed[~np.all(contact_map_completed == 0, axis=1)]
       indices_complete = list(range(np.max(indices)+1))
       tcr_completed = ''
       counter =  0
       for k in range(len(indices_complete)):
              if np.isin(indices_complete, indices)[k]:
                     tcr_completed+=tcr[counter]
                     counter+=1
              else:
                     tcr_completed+='X'

       fig, ax = plt.subplots()
       sns.heatmap(contact_map_completed, ax = ax, xticklabels=list(row['peptide']), yticklabels=list(tcr_completed), square = True, cmap = 'coolwarm', cbar = False)
       ax.hlines(l_t_a + 1, ax.get_xlim()[0], ax.get_xlim()[1], color = 'white', ls = 'dashed', lw = 3)
       fig.savefig("../../../Figures/repertoire_entropy/" + row['gene'] + "_" + tcr_completed + "_" + row['pdb_id'] + ".png")

       ids_dict['tcr_a'].append(tcr_a)
       ids_dict['index_a'].append(indices_a)
       ids_dict['tcr_b'].append(tcr_b)
       ids_dict['index_b'].append(indices_b)
       ids_dict['tcr'].append(tcr)
       ids_dict['contact_map'].append(contact_map)
       ids_dict['tcr_completed'].append(tcr_completed)
       ids_dict['contact_map_completed'].append(contact_map_completed)
       #Calciulate consensus energies:
       peptide = row['peptide']
       pwm_df = potential_df[list(peptide)]
       E_matrix_df = pwm_df.reindex(list(tcr))
       EW_matrix = np.array(E_matrix_df)*contact_map
       e = np.sum(EW_matrix)
       ids_dict['consensus_energies'].append(e)

ids_df = pd.DataFrame(ids_dict).sort_values('peptide').reset_index(drop=True)

'''
Until here, the code just curates the contact map database.
It discards all contact maps that does not belong to the three
antigen epitopes present in the epitope-specicic tcr repertoire
database. Now, we proceed to cluster epitope-specific tcr
repertoires and look for correlations between clusters and 
contact maps. 
'''

'''
The first thing to try is to calculate binding energy for all tcr
clusters using all available contact maps.
'''

reduced_tcrs_dict = defaultdict(list)

for i, row in dash_db.iterrows():
       gene = row['epitope']
       cluster = row['cluster']
       # print(gene)
       tcr_dataset_a = row['cdr3_a_aa'] 
       tcr_dataset_b = row['cdr3_b_aa']
       tcr_dataset = tcr_dataset_a + tcr_dataset_b
       tcr_dataset_a_list = np.array(list(tcr_dataset_a))
       tcr_dataset_b_list = np.array(list(tcr_dataset_b))
       
       tcr_test_df = ids_df.loc[ids_df['gene']==gene]
       for k, row_test in tcr_test_df.iterrows():
              tcr_test_a = row_test['tcr_a']
              tcr_test_b = row_test['tcr_b']
              alignment_a = needleman_wunsch(tcr_dataset_a, tcr_test_a)
              alignment_b = needleman_wunsch(tcr_dataset_b, tcr_test_b)
              joint_positions = np.array([k[0] for k in alignment_a if k[1] != None] + [k[0] for k in alignment_b if k[1] != None])
              if (joint_positions==None).any():
                     break
              reduced_tcr_dataset_a = ''.join(tcr_dataset_a_list[[k[0] for k in alignment_a if k[1] != None]])
              reduced_tcr_dataset_b = ''.join(tcr_dataset_b_list[[k[0] for k in alignment_b if k[1] != None]])
              reduced_tcr_dataset = reduced_tcr_dataset_a + reduced_tcr_dataset_b
              # print(reduced_tcr_dataset_a, reduced_tcr_dataset_b)
              reduced_tcrs_dict['db_id'].append(i)
              reduced_tcrs_dict['c'].append(cluster)
              reduced_tcrs_dict['gene'].append(gene)
              reduced_tcrs_dict['pdb_id'].append(row_test['pdb_id'])
              reduced_tcrs_dict['peptide'].append(row_test['peptide'])
              reduced_tcrs_dict['reduced_tcr'].append(reduced_tcr_dataset)
              reduced_tcrs_dict['contact_map'].append(row_test['contact_map'])
              
reduced_tcrs_df = pd.DataFrame(reduced_tcrs_dict)

energies = []
for i, row_dash in reduced_tcrs_df.iterrows():
       pdb_id = row_dash['pdb_id']
       tcr = row_dash['reduced_tcr']
       gene = row_dash['gene']
       peptide = row_dash['peptide']
       contact_map = row_dash['contact_map']
       pwm_df = potential_df[list(peptide)]
       E_matrix_df = pwm_df.reindex(list(tcr))
       EW_matrix = np.array(E_matrix_df)*contact_map
       e = np.sum(EW_matrix)
       energies.append(e)

reduced_tcrs_df['e'] = energies

for g, gene in enumerate(genes):

       cluster = clusters[gene]
       organism = organisms[g]
       
       fig, ax = plt.subplots()
       ax.set_title('Organism:' + organism + " ; Gene:" + gene)
       sns.histplot(reduced_tcrs_df.loc[reduced_tcrs_df['gene']==gene], x = 'e', ax=ax, hue = 'pdb_id', bins = np.linspace(-6, 6, 20), multiple = 'stack')
       fig.savefig("../../../Figures/repertoire_entropy/histogram_energies_" + gene + ".png")

       fig, ax = plt.subplots()
       ax.set_title('Organism:' + organism + " ; Gene:" + gene)
       sns.histplot(reduced_tcrs_df.loc[(reduced_tcrs_df['gene']==gene) & (reduced_tcrs_df['c'].isin(larger_clusters[gene]))], x = 'e', ax=ax, hue = 'c', bins = np.linspace(-6, 6, 20), palette = 'rainbow', multiple = 'stack')
       fig.savefig("../../../Figures/repertoire_entropy/histogram_energies_" + gene + "_c.png")

       pdb_structures = reduced_tcrs_df.loc[reduced_tcrs_df['gene']==gene, 'pdb_id'].unique()
       for pdb_id in pdb_structures:

            fig, ax = plt.subplots()
            ax.set_title('Organism:' + organism + " ; Gene:" + gene + " ; Structure:" + pdb_id)
            sns.histplot(reduced_tcrs_df.loc[(reduced_tcrs_df['gene']==gene) & (reduced_tcrs_df['pdb_id']==pdb_id)], x = 'e', ax=ax, hue = 'c', bins = np.linspace(-6, 6, 20), multiple = 'stack', legend = False)
            ax.vlines(ids_df.loc[ids_df['pdb_id']==pdb_id, 'consensus_energies'], ax.get_ylim()[0], ax.get_ylim()[1], color = 'black', ls = 'dashed', lw = 3)
            fig.savefig("../../../Figures/repertoire_entropy/histogram_energies_" + gene + "_" + pdb_id +".png")





