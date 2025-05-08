import sys
sys.path.append('../../lib/')
from funcs import*
plt.rcParams['text.usetex'] = True

project = 'epitope_complexity'
subproject = 'data'
experiment = 0
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}"
output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
os.makedirs(output_plot, exist_ok=True)

# Parameters
n_ens = 10000
gs = [2]  # Number of Poisson processes
mu = 1.0  # Poisson rate
T = 15  # Total simulation time
theta = 1.5  # Values of theta to compare
gamma = 0.4
my_colors = [my_red, my_cyan, my_green, my_green, my_brown]
my_colors2 = ['darkred', my_blue2, my_green2, my_green2]
alpha = 1e-10
depth = 6
anti_mut_epi = 5/4

data = pd.read_csv(root_dir + "/df_potency.txt", sep="\t")
dataID = pd.read_csv(root_dir + "/df_immuno_dominance.txt", sep="\t")

print(data['muts_A'])

# data['Z'] = np.exp(data['dE_A']) + np.exp(data['dE_B']) + np.exp(data['dE_C']) + np.exp(data['dE_D']) + np.exp(data['dE_E'])

chi = []

for index, row in data.iterrows():
    DDGs = row[['dE_A', 'dE_B', 'dE_C', 'dE_D', 'dE_E']].astype(float).values
    ferret = row['ferret']
    
    # Get the corresponding weights row
    ID_row = dataID[dataID['ferret'] == ferret][['A', 'B', 'C', 'D', 'E']]
    
    # Make sure we only proceed if we found a match
    if not ID_row.empty:
        weights = ID_row.iloc[0].astype(float).values  # Extract the row as a 1D NumPy array
        chi.append(np.sum(np.exp(-2*DDGs) * weights))
    else:
        print(f"Warning: No match for ferret '{ferret}' in dataID")

data['log_chi'] = np.log(chi)

data.to_csv(root_dir + "/df_potency_chi.txt")

# sns.violinplot(data = data, x = 'chi', hue = 'ferret', legend = False)
# plt.show()
# # data = data[data['epi_dist']<10]
# data['muts'] = data['muts'].apply(literal_eval)
# avg_data = data.groupby(['ferret', 'epi_dist'], as_index=False)['dh'].mean()

# grand_avg = avg_data.groupby('epi_dist', as_index=False)['dh'].mean()

# # Get unique ferrets
# ferrets = data['ferret'].unique()


# fig0, ax0 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.1, 'right':.95, 'bottom':.1, 'top': 0.96})
# for f, ferret in enumerate(ferrets):
#     fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.15, 'top': 0.96})
#     # Raw points
#     sub_df = data[data['ferret'] == ferret]
#     ax.scatter(sub_df['epi_dist']/anti_mut_epi, -np.log(2**sub_df['dh']), alpha=0.4, label=f'{ferret} (raw)')

#     # Average line
#     sub_avg = avg_data[avg_data['ferret'] == ferret]
#     ax.step(sub_avg['epi_dist']/anti_mut_epi, -np.log(2**sub_avg['dh']), marker='', linestyle='-', label=f'{ferret} (avg)')
#     ax0.step(sub_avg['epi_dist']/anti_mut_epi, -np.log(2**sub_avg['dh']), marker='', linestyle='-', label=f'{ferret} (avg)', alpha = 0.4)

#     ax.plot(sub_avg['epi_dist']/anti_mut_epi, -sub_avg['epi_dist']/anti_mut_epi, marker='', linestyle='--', label=f'{ferret} (avg)')

#     # Final plot touches
#     my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
#     ax.set_xlabel('epi_dist')
#     ax.set_ylabel('dh')
#     # ax.legend()
#     fig.savefig(output_plot + '/T_epi_dist_'+str(f)+'.pdf')

# # ax0.plot(grand_avg['epi_dist']/anti_mut_epi, -grand_avg['epi_dist']/anti_mut_epi, color='black', marker='', linestyle='--', label=f'{ferret} (avg)')
# ax0.plot(grand_avg['epi_dist']/anti_mut_epi, -np.log(2**grand_avg['dh']), color='black', linestyle='-', linewidth=2, label='Grand average')

# # Final plot touches
# my_plot_layout(ax = ax0, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
# # ax0.set_xlabel(r'$\textrm{distance}/5$')
# # ax0.set_ylabel(r'$\textrm{dh}$')
# # ax0.legend()
# fig0.savefig(output_plot + '/T_epi_dist.pdf')

