import sys
sys.path.append('../../lib/')
from functions_cm import*
from classes import*
plt.rcParams['text.usetex'] = True
import time
import itertools
import ast
from collections import defaultdict
import re
#import logomaker

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


l = 4
ps = [0.1, 0.2, 0.5, 0.8]

#antigen_random = np.random.choice(Alphabet_p, l)
antigens = ['TTMM', 'AAAR', 'WWWR', 'NNNN', 'RKDL', 'NCYH', 'ARRR']

# mean_df = pd.read_csv('../../out/tcr_contact_maps/mean_df.csv', index_col = 0)
# var_df = pd.read_csv('../../out/tcr_contact_maps/var_df.csv', index_col = 0)
# tcr_m_df = pd.read_csv('../../out/tcr_contact_maps/tcr_m_df.csv', index_col = 0)
# cm_m_df = pd.read_csv('../../out/tcr_contact_maps/cm_m_df.csv', index_col = 0)

tcrs_m = defaultdict(list)
cms_m = defaultdict(list)

fig, ax = plt.subplots(figsize = (10, 8))
for antigen in antigens:
    tcr_m_df = pd.read_csv('../../out/tcr_contact_maps/tcr_m_df_'+antigen+'.csv', index_col = 0)
    cm_m_df = pd.read_csv('../../out/tcr_contact_maps/cm_m_df_'+antigen+'.csv', index_col = 0)
    for p in ps:
        for elem in tcr_m_df[str(p)]:
            tcrs_m["antigen"].append(antigen)
            tcrs_m["p"].append(p)
            # tcrs_m["value"].append([''.join(a) for a in ast.literal_eval(tcr_m_df[str(p)])])
            tcrs_m["value"].append(elem)
        for elem in cm_m_df[str(p)]:
            cms_m["antigen"].append(antigen)
            cms_m["p"].append(p)
            # cms_m["value"].append(ast.literal_eval(cm_m_df[str(p)]))
            cms_m["value"].append(elem)

tcrs_df = pd.DataFrame(tcrs_m)
cms_df = pd.DataFrame(cms_m)

e_motif = defaultdict(list)
for antigen in antigens:
    temp_df1 = tcrs_df.loc[tcrs_df['antigen']==antigen]
    for p in ps:
        temp_df2 = temp_df1.loc[temp_df1['p']==p]
        for elem in temp_df2['value']:
            aa_seq = ast.literal_eval(elem)
            for pos in range(l):
                e_motif["antigen"].append(antigen)
                e_motif["p"].append(p)
                e_motif["pos"].append(pos+1)
                # e_motif["aa"].append(aa_seq[pos])
                e_motif['aa'].append(list(Alphabet_T).index(aa_seq[pos]))
e_motif_df = pd.DataFrame(e_motif)

cm_motif = defaultdict(list)
for antigen in antigens:
    temp_df1 = cms_df.loc[cms_df['antigen']==antigen]
    for p in ps:
        temp_df2 = temp_df1.loc[temp_df1['p']==p]
        for elem in temp_df2['value']:
            cm = ast.literal_eval(','.join(re.sub(r'(?<=\d)(\s+)(?=-?\d)', ',', elem).splitlines()))
            for pos_t in range(l):
                for pos_p in range(l):
                    if cm[pos_t][pos_p]:
                        cm_motif["antigen"].append(antigen)
                        cm_motif["p"].append(p)
                        cm_motif["pos_p"].append(pos_p+1)
                        cm_motif["pos_t"].append(pos_t+1)
cm_motif_df = pd.DataFrame(cm_motif)

for antigen in antigens:
    for p in ps:
        fig, ax = plt.subplots(figsize = (10, 8))
        ax.set_title(antigen + r' ; $p=%.1f$'%p, fontsize = 20)
        sns.histplot(data=e_motif_df.loc[(e_motif_df['antigen']==antigen) & (e_motif_df['p']==p)], ax = ax, x="pos", y="aa", discrete = True, cbar = True, bins = (range(1, l+1), range(1, len(Alphabet_T)+1)))
        ax.tick_params(labelsize = 20)
        ax.set_xlabel('Position peptide', fontsize = 20)
        ax.set_ylabel('aa', fontsize = 20)
        ax.set_xticks(np.arange(l)+1)
        ax.set_yticks(np.arange(len(Alphabet_T))+1)
        ax.set_yticklabels(Alphabet_T)
        fig.savefig('../../../Figures/tcr_contact_map/motifs/e_motif_'+antigen+'_p-%.1f.pdf'%(p))
        plt.close(fig)

for antigen in antigens:
    for p in ps:
        fig, ax = plt.subplots(figsize = (10, 8))
        ax.set_title(antigen + r' ; $p=%.1f$'%p, fontsize = 20)
        sns.histplot(data=cm_motif_df.loc[(cm_motif_df['antigen']==antigen) & (cm_motif_df['p']==p)], ax = ax, x="pos_p", y="pos_t", discrete = True, cbar = True, bins = (range(1, l+1), range(1, l+1)))
        ax.tick_params(labelsize = 20)
        ax.set_xlabel('Position peptide', fontsize = 20)
        ax.set_ylabel('Position tcr', fontsize = 20)
        ax.set_xticks(np.arange(1, l+1))
        ax.set_yticks(np.arange(l)+1)
        fig.savefig('../../../Figures/tcr_contact_map/motifs/cm_motif_'+antigen+'_p-%.1f.pdf'%(p))
        plt.close(fig)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time/60:.3f} minutes")



