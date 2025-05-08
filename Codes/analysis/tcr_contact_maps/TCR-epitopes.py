import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
import itertools
import ast
from collections import defaultdict
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

database = pd.read_table('../../in/epitopes-TCRs.tsv')
database = pd.DataFrame(database)
print(database.columns)
database = database[['CDR3', 'Epitope']]
CDR3s = database['CDR3'].unique()
epitopes = database['Epitope'].unique()
print(len(CDR3s), len(epitopes))

pairs = defaultdict(list)

for cdr3, epitope in zip(database['CDR3'], database['Epitope']):
	pairs['epitope'].append(epitope)
	pairs['cdr3'].append(cdr3)

print([database['CDR3'][i] for i in range(len(database['Epitope'])) if database['Epitope'][i] == 'LLFGYPVYV'])

