import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

c_maps_df = pd.read_csv('../../in/contact_maps/contact_maps_PDB.csv')
miss_df = pd.read_csv('../../in/contact_maps/summary_PDB_structures.csv')
max_pos_t = np.max(c_maps_df['residue.index.from'])
max_pos_p = np.max(c_maps_df['residue.index.to'])
print(max_pos_t)
c_map_ids = np.unique(c_maps_df['pdb.id'])
c_maps = []
contacs = []
contacs_t = []
contacs_p = []
MHC = miss_df['mhc.class']
for i_map in tqdm(range(len(c_map_ids))):
	if(miss_df.loc[miss_df['pdb.id']==c_map_ids[i_map]]['mhc.class'].values[0]=='MHCI'):
		c_map_df = c_maps_df.loc[c_maps_df['pdb.id']==c_map_ids[i_map]]
		p_pos = list(np.unique(c_map_df['residue.index.to']))
		t_pos = list(np.unique(c_map_df['residue.index.from']))
		c_map = np.zeros((max_pos_p+1, max_pos_t+1))
		for l in range(len(c_map_df['pdb.id'])):
			index_p = int(np.array(c_map_df['residue.index.to'])[l])
			# print(index_p)
			index_t = t_pos.index(np.array(c_map_df['residue.index.from'])[l])
			# print(index_t)
			c_map[np.array(c_map_df['residue.index.to'])[l], np.array(c_map_df['residue.index.from'])[l]] = 1
		c_maps.append(c_map)
		contacs.append(np.sum(c_map))
		contacs_t+= list(np.sum(c_map, axis = 0))
		contacs_p+= list(np.sum(c_map, axis = 1))
		#fig, ax = plt.subplots()
		#map_df_new = pd.DataFrame(c_map)
		#sns.heatmap(map_df_new, ax = ax, cmap = 'binary', cbar = False)
		#fig.savefig('../../../Figures/tcr_contact_map/contact_maps/c_map_' + np.array(c_map_df['pdb.id'])[0] + '.png')

fig, ax = plt.subplots()
data_hist = np.histogram(contacs, bins = np.arange(0, 60, 1), density = False)

P_c = np.diff(np.cumsum(data_hist[0]))/np.diff(data_hist[1][:-1])#/(len(energies_lineages))
P_c /= np.sum(P_c*np.diff(data_hist[1][:-1]))


ns = np.arange(0, 60)
ax.plot(np.arange(0, 60, 1)[:-2], P_c, color = my_red, marker = '*', ls = '--')
ax.plot(ns, np.array([np.mean(contacs)**n/math.factorial(n)*np.exp(-np.mean(contacs)) for n in ns]), color = my_red)
fig.savefig('../../../Figures/tcr_contact_map/contacts.png')

fig, ax = plt.subplots()
ax.hist(contacs_t, bins = np.arange(0, 12, 1), label = 'T', alpha = .6)
ax.hist(contacs_p, bins = np.arange(0, 12, 1), label = 'p', alpha = .6)
ax.legend()
ax.set_yscale('log')
# ax.set_xscale('log')
fig.savefig('../../../Figures/tcr_contact_map/contacts_cols.png')