import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
import itertools
import ast
# import pyrepseq.plotting as rsp
# import pyrepseq.distance as rsd
from collections import defaultdict
import scipy.cluster.hierarchy as hc
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")



data = pd.read_excel('/Users/robertomorantovar/Dropbox/Research/Immune_system/repertoire_entropy/data/victora_2020/mmc1.xlsx', header=1, sheet_name = 'Photoactivation CGG')
data_early = data.loc[(data['Figure']==1)]
data_early = data_early.dropna()
seqs, counts = np.unique(data_early['CDR3:'], return_counts = True)
lengths = [len(i) for i in seqs]
plt.hist(lengths, bins = range(1, 30))
plt.show()

print(needleman_wunsch(seqs[4], seqs[20])[3][1])