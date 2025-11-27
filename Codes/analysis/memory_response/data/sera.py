import sys
sys.path.append('../../../my_lib/')
from funcs import*
plt.rcParams['text.usetex'] = True

project = 'memory_response'
subproject = 'data'
experiment = 0
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/schiepers2023"
output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+ subproject + '/schiepers2023/sera'
os.makedirs(output_plot, exist_ok=True)

# Parameters
n_ens = 10000
gs = [2]  # Number of Poisson processes
mu = 1.0  # Poisson rate
T = 15  # Total simulation time
theta = 1.5  # Values of theta to compare
gamma = 0.4
alpha = 1e-10
depth = 6
anti_mut_epi = 5/4
n_ensemble = 1000

zeta_min = 0.3
zeta_max = 1.4
color_vals = np.linspace(0, 1, int((zeta_max - zeta_min)*100))
cmap = plt.get_cmap('rainbow_r')
my_colors_alpha = [cmap(val) for val in color_vals] 
my_colors_alpha_proposal = [my_blue, my_red]
my_colors = [my_blue2, my_purple, my_purple, my_purple, my_cyan]

boosts = [0, 32, 64]

fig, ax = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
# color per top-level header
color_map = {
    'Anti-TNP IgG titer': 'grey',
    'Anti-TNP IgK titer': 'black',
    'Anti-TNP FLAG titer': my_blue,
    'Anti-TNP Strep titer': my_red
}
#------------ Experiment 1 (Figure 1D) ------------
data_sera = pd.read_excel(root_dir + "/Fig_2.xlsx", sheet_name = 'Fig 2c-d', header = [0,1])
# Drop the first column you don't want
data_sera = data_sera.drop(data_sera.columns[0], axis=1)
data_sera = data_sera.set_index(data_sera.columns[0])
data_sera = data_sera[['Anti-TNP FLAG titer', 'Anti-TNP Strep titer']]
print(data_sera.columns)
for (top, sub) in data_sera.columns:
    ax.plot(
        data_sera.index.to_numpy()[1:],    # x-axis = time
        (data_sera[(top, sub)].to_numpy()[1:]-data_sera[(top, sub)].to_numpy()[:-1])/data_sera[(top, sub)].to_numpy()[2],    # y-axis = the column
        label=f"{sub} ({top})",
        color=color_map[top], ls = '', marker='o', markersize=8
    )
ax.vlines(boosts, ymin=1e-1, ymax=1e3, colors='grey', linestyles='dashed', label='Boosts')
ax.hlines(1, xmin=0, xmax=np.max(data_sera.index.to_numpy()[:]), colors='grey', linestyles='dashed', label='Baseline')
ax.set_yscale('log')
plt.show()

