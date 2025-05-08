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

E_model = 'TCRen'
#------- TCRen ------
potential_csv = pd.read_csv('../../in/TCRen_potential.csv')
Alphabet_T = np.unique(potential_csv['residue.aa.from'])
Alphabet_p = np.unique(potential_csv['residue.aa.to'])
#print(Alphabet_T, Alphabet_p)

l = 4
ps = [0.1, 0.2, 0.5, 0.8]

#antigen_random = np.random.choice(Alphabet_p, l)
antigens = ['TTMM', 'AAAR', 'WWWR', 'NNNN', 'RKDL', 'NCYH', 'ARRR']
makers_antigen = ['o', '^', '*', 's', 'd', '.', 'v']
colors_antigen = [my_blue, my_red, my_green, my_brown, my_blue2, my_purple2, my_cyan]

# mean_df = pd.read_csv('../../out/tcr_contact_maps/mean_df.csv', index_col = 0)
# var_df = pd.read_csv('../../out/tcr_contact_maps/var_df.csv', index_col = 0)
# tcr_m_df = pd.read_csv('../../out/tcr_contact_maps/tcr_m_df.csv', index_col = 0)
# cm_m_df = pd.read_csv('../../out/tcr_contact_maps/cm_m_df.csv', index_col = 0)

means = defaultdict(list)
varss = defaultdict(list)
means_summ = defaultdict(list)
varss_summ = defaultdict(list)
tcrs_m = defaultdict(list)
cms_m = defaultdict(list)

fig, ax = plt.subplots(figsize = (10, 8))
for antigen in antigens:
    mean_df = pd.read_csv('../../out/tcr_contact_maps/mean_df_'+antigen+'.csv', index_col = 0)
    var_df = pd.read_csv('../../out/tcr_contact_maps/var_df_'+antigen+'.csv', index_col = 0)
    for p in ps:
        means_summ["antigen"].append(antigen)
        means_summ["p"].append(p)
        means_summ["value"].append(np.mean(mean_df[str(p)]))
        varss_summ["antigen"].append(antigen)
        varss_summ["p"].append(p)
        varss_summ["value"].append(np.mean(var_df[str(p)]))
        #for elem in ast.literal_eval(mean_df[str(p)]):
        for elem in mean_df[str(p)]:
            means["antigen"].append(antigen)
            means["p"].append(p)
            means["value"].append(elem)
        # for elem in ast.literal_eval(var_df[str(p)]):
        for elem in var_df[str(p)]:
            varss["antigen"].append(antigen)
            varss["p"].append(p)
            varss["value"].append(elem)
        
means_df = pd.DataFrame(means)
varss_df = pd.DataFrame(varss)
means_summ_df = pd.DataFrame(means_summ)
varss_summ_df = pd.DataFrame(varss_summ)

fig, ax = plt.subplots(figsize = (10, 8))
sns.lineplot(data=means_df, ax = ax, x="p", y="value", hue="antigen", err_style="bars", errorbar=("sd", 4),)
fig.savefig('../../../Figures/tcr_contact_map/means.pdf')

p_array = np.linspace(0, 1 , 20)
fig, ax = plt.subplots(figsize = (10, 8))
for i, antigen in enumerate(antigens):
    df_temp = means_summ_df.loc[means_summ_df['antigen']==antigen]
    sns.lineplot(ax = ax, x=ps, y= df_temp['value']/np.array(df_temp['value'])[0], ls = '', marker = makers_antigen[i], ms = 10, label = antigen)
    sns.lineplot(ax = ax, x=p_array, y= p_array/0.1, ls = '--', color = 'gray')
    ax.tick_params(labelsize = 20)
    ax.set_xlabel(r'$p$', fontsize = 20)
    ax.set_ylabel(r'$\bar{\epsilon}$', fontsize = 20)
    fig.savefig('../../../Figures/tcr_contact_map/means_summ_colapse.pdf')
fig, ax = plt.subplots(figsize = (10, 8))
for i, antigen in enumerate(antigens):
    df_temp = means_summ_df.loc[means_summ_df['antigen']==antigen]
    sns.lineplot(ax = ax, x=ps, y= df_temp['value'], ls = '', marker = makers_antigen[i], ms = 10, color = colors_antigen[i], label = antigen)
    sns.lineplot(ax = ax, x=p_array, y= (1/0.8)*np.array(df_temp['value'])[-1]*np.array(p_array), ls = '--', color = colors_antigen[i])
    print((1/0.8)*np.array(df_temp['value'])[-1], (1/0.1)*np.array(df_temp['value'])[0])
    ax.tick_params(labelsize = 20)
    ax.set_xlabel(r'$p$', fontsize = 20)
    ax.set_ylabel(r'$\bar{\epsilon}$', fontsize = 20)
    fig.savefig('../../../Figures/tcr_contact_map/means_summ.pdf')


fig, ax = plt.subplots(figsize = (10, 8))
for i, antigen in enumerate(antigens):
    df_temp = varss_summ_df.loc[varss_summ_df['antigen']==antigen]
    sns.lineplot(ax = ax, x=ps, y= df_temp['value'], ls = '', marker = makers_antigen[i], ms = 10, label = antigen, color = colors_antigen[i])
    popt, pcov = curve_fit(my_linear_func, np.log(ps), np.log(df_temp['value']/np.array(df_temp['value'])[0]))
    sns.lineplot(ax = ax, x=p_array, y= (1/0.8)**popt[1]*np.array(df_temp['value'])[-1]*np.array(p_array)**(popt[1]), ls = '--', color = colors_antigen[i])

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(0.09, 1.005)
ax.set_ylim(bottom = 0.2, top = 10)
ax.tick_params(labelsize = 20)
ax.set_xlabel(r'$p$', fontsize = 20)
ax.set_ylabel(r'$\sigma^2_{\epsilon}$', fontsize = 20)
fig.savefig('../../../Figures/tcr_contact_map/vars.pdf')

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time/60:.3f} minutes")



