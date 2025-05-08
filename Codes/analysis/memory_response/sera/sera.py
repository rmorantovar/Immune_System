'''
author=mmeijers
created on 04-10-2022
'''
from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
import flai
import os
import glob
import sqlite3
from tqdm import tqdm
from flai.util.vir.flu.DriverSites import DriverSites
from flai.util.Utils import Utils
from flai.util.Time import Time
import sys
from flai.tree.LeanTree import LeanTree
from flai.tree.LeanCladeTree import LeanCladeTree
    
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
args = parser.parse_args()


args.tree_db = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Temp_Tree_Share/H3N2/2024_07_02/tree.db'
args.output = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Codes/analysis/memory_response/sera/out'
args.lineage = 'H3N2'
DriverSites.setH3N2()

args.c = ['CRICK', 'CDC', 'VIDRL']
args.a = ['FRAVNT', 'HINT', 'HI']
args.sm = ['CELL']

con = sqlite3.connect(args.tree_db, timeout=1000)
cur = con.cursor()
df_tree = pd.read_sql_query("SELECT nid, pid, rcluster, clade, db_sequence, mutations, name,collection_date, multiplicity,epi,region, leaf_count FROM Tree",con)
df_clades = pd.read_sql_query("SELECT * FROM SignificantClades WHERE frequency_threshold = 0.05",con)
df_antigenic = pd.read_sql_query('SELECT strainNode, rstrainNode,'
                                 'alphaClade, rhoClade, epi_dist, titre,'
                                 'centre, assay,serumMedium, train '
                                 'FROM AntigenicData', con)
df_Ear = pd.read_sql_query("SELECT rho, alpha, Ealpharho_ALL_ALL_CELL FROM Ealpharho",con)
con.close()
    
clade2sclade = {a[0]:a[1] for a in zip(list(df_clades.clade),list(df_clades.antigenic_clade))}
clade2label = {a:b for a,b in zip(df_clades.clade, df_clades.clade_label)}
df_tree['antigenic_clade'] = [clade2sclade[a] for a in df_tree.clade]
df_tree['FLAI_time'] = [Time.dateToCoordinate(a) for a in df_tree.collection_date]

tt = LeanTree()
tt.initialize(df_tree)
tt.set_significant_clades(clade2sclade)
tt.compute_node_frequencies()

clade_tree = LeanCladeTree(tt, significantAlphas=True)
clade_tree.compute_clade_frequencies()

# set train-test beforehand
df_antigenic = df_antigenic.loc[[c in args.c for c in df_antigenic.centre]]
df_antigenic = df_antigenic.loc[[a in args.a for a in df_antigenic.assay]]
df_antigenic = df_antigenic.loc[[sm in args.sm
                                 for sm in df_antigenic.serumMedium]]
df_antigenic['alphaClade'] = [tt.nodes[nid].cladeSalpha for nid in df_antigenic.strainNode]
df_antigenic['rhoClade'] = [tt.nodes[nid].cladeSalpha for nid in df_antigenic.rstrainNode]
df_antigenic['mesmul'] = np.ones(len(df_antigenic))

# a = defaultdict(lambda: 0.0)


