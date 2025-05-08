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


#---------- import contact map info ----------
cmaps_db = pd.read_csv('../../in/contact_maps/contact_maps_PDB.csv')
miss_cmaps_db = pd.read_csv('../../in/contact_maps/summary_PDB_structures.csv')
dash_db = pd.read_csv('../../in/dash_human.csv')
print(dash_db, cmaps_db)
peptides = ['GLCTLVAML', 'NLVPMVATV', 'GILGFVFTL']
genes = ['BMLF', 'pp65', 'M1']


ids_dict = defaultdict(list)

for index, row in miss_cmaps_db.iterrows():
	p = row['peptide']
	if p in peptides:
		i_p = peptides.index(p)
		ids_dict['peptide'].append(p)
		ids_dict['pdb.id'].append(row['pdb.id'])
		ids_dict['gene'].append(genes[i_p])
ids_df = pd.DataFrame(ids_dict).sort_values('peptide').reset_index(drop=True)
print(ids_df)


# peptide_id = [np.array(miss_df['pdb.id'])[i] for i in range(len(miss_df)) if np.array(miss_df['peptide'])[i] in peptides]
# print(peptide_id)
# epitope_db = miss_df.loc[miss_df['pdb.id']=='3o4l']['peptide'][69]

# cmap_a = c_maps_df.loc[(c_maps_df['pdb.id']=='3o4l') & (c_maps_df['chain.id.from']=='D') & (c_maps_df['region.type.from']=='CDR3')]
# cmap_b = c_maps_df.loc[(c_maps_df['pdb.id']=='3o4l') & (c_maps_df['chain.id.from']=='E') & (c_maps_df['region.type.from']=='CDR3')]

# l_t_a = len(cmap_a['residue.index.from'].unique())
# l_t_b = len(cmap_b['residue.index.from'].unique())
# l_t = l_t_a + l_t_b

# cmap_a['residue.index.from'] -= np.min(cmap_a['residue.index.from'])
# cmap_b['residue.index.from'] -= np.min(cmap_b['residue.index.from']) 
# cmap_b['residue.index.from'] += l_t_a

# cmap_merged = pd.concat([cmap_a, cmap_b])[['chain.type.from', 'residue.index.from', 'residue.aa.from', 'residue.index.to', 'residue.aa.to']].sort_values('residue.index.from')



# epitope_gene = 'M1'

# database = pd.read_csv('../../in/dash_human.csv')
# #database = pd.DataFrame(database)
# #database = database.loc[database['epitope'] == epitope_gene]

# epitope_genes = database['epitope'].unique()

# print(epitope_genes)

# 	