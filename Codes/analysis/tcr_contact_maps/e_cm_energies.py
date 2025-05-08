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
energies_m = defaultdict(list)
for antigen in antigens:
    tcr_m_df = pd.read_csv('../../out/tcr_contact_maps/tcr_m_df_'+antigen+'.csv', index_col = 0)
    cm_m_df = pd.read_csv('../../out/tcr_contact_maps/cm_m_df_'+antigen+'.csv', index_col = 0)
    pwm_df = potential_df[list(antigen)]
    for p in ps:
        for tcr_raw, cm_raw in zip(tcr_m_df[str(p)], cm_m_df[str(p)]):
            tcr = ast.literal_eval(tcr_raw)
            E_matrix_df = pwm_df.reindex(list(tcr))
            cm = ast.literal_eval(','.join(re.sub(r'(?<=\d)(\s+)(?=-?\d)', ',', cm_raw).splitlines()))
            EW_matrix = np.array(E_matrix_df)*cm
            E = np.sum(EW_matrix)
            energies_m['antigen'].append(antigen)
            energies_m['p'].append(p)
            energies_m['e_m'].append(E)

energies_m_df = pd.DataFrame(energies_m)    

for antigen in antigens:
    fig, ax = plt.subplots(figsize = (10, 8))
    ax.set_title(antigen, fontsize = 20)
    sns.histplot(data=energies_m_df.loc[(energies_m_df['antigen']==antigen)], ax = ax, x="e_m", hue="p", bins = 80)
    ax.tick_params(labelsize = 20)
    ax.set_xlabel(r'$\textrm{Energy}$', fontsize = 20)
    ax.set_ylabel('counts', fontsize = 20)
    # ax.set_xticks(np.arange(l)+1)
    # ax.set_yticks(np.arange(len(Alphabet_T))+1)
    # ax.set_yticklabels(Alphabet_T)
    fig.savefig('../../../Figures/tcr_contact_map/e_m_'+antigen+'_p-%.1f.pdf'%(p))
    plt.close(fig)


end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time/60:.3f} minutes")



