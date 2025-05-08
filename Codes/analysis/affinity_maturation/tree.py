from Bio import Phylo
import pandas as pd
import ast
import numpy as np

metadf = pd.read_csv('../../out/affinity_maturation/metadata_0.csv')
parents = metadf['parent'].unique()
daughters = metadf['id'].tolist()
parents_dict = {p:[] for p in parents}
#print(metadf)
for row in metadf.itertuples():
	d = row.id
	p = row.parent
	if p in parents:
		parents_dict[p].append(d)


def add_subtree(d):
	subtree = ('node_'+str(d),)
	for dd in parents_dict[d]:
		if dd not in parents:
			subtree+=('node_'+str(dd),)
			daughters_temp.remove(dd)
		else:
			subsubtree = add_subtree(dd)
			subtree+=(subsubtree,)
			daughters_temp.remove(dd)
	return subtree

tree = ()
daughters_temp = [d for d in daughters]
for d in daughters:
	if d in daughters_temp:
		if d not in parents:
			tree+=('node_'+str(d),)
		else:
			subtree = add_subtree(d)
			tree+=(subtree,)
print(tree)
# np.savetxt('../../out/affinity_maturation/tree.newick')
