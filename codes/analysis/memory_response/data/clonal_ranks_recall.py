import sys
sys.path.append('../../../my_lib/')
from funcs import*
plt.rcParams['text.usetex'] = True

project = 'memory_response'
subproject = 'data'
experiment = 0
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/mesin2020"
output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+ subproject + '/mesin2020/recall'
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
markers_ph = ['o', 's', '^', 'D', 'v', 'P', '*']

def model(x, m):
    return m * x 

fig_r2, ax_r2 = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
fig_zeta2, ax_zeta2 = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})

#------------ Experiment 2 and 3 (Figure 4A and 4C) ------------

data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
phs = data_recall['Phenotype'].unique()
print(phs)
figures = ['4A', '4C-H']

max_rank = 100
max_rank_fit = 20

for i_ph, ph in enumerate(phs):
    zetas = []
    data_ph = data_recall[(data_recall['Phenotype']==ph)]
    print(ph)
    for rep in tqdm(range(n_ensemble)):
        x_avg = np.zeros(max_rank)
        counts_per_ranking = np.zeros(max_rank)
        min_max_rank_mouse = max_rank
        max_max_rank_mouse = 0
        len_mice = 0

        for i_fig, fig in enumerate(figures):
            data_ph_fig = data_ph[(data_ph['Figure']==fig)]
            data_ph_grouped = data_ph_fig.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
            # data_recall_grouped = data_recall_fig.groupby(['Mouse', 'CDR3:']).size().reset_index(name='count')
            # data_recall_grouped = data_recall_fig.groupby(['Mouse', 'Sequence']).size().reset_index(name='count')
            mice = data_ph_grouped['Mouse'].unique()
            if rep == n_ensemble - 1:
                mice_rep = mice
            else:
                mice_rep = np.random.choice(mice, len(mice), replace = True)

            len_mice+=len(mice)
            for mouse in mice_rep:
                data_mouse = data_ph_grouped[data_ph_grouped['Mouse']==mouse]
                # CDR3, counts = np.unique(np.array((list(data_mouse['CDR3:']))), return_counts = True)
                counts = data_mouse['count'].to_numpy()
                sort_index = counts.argsort()
                largest = np.max(counts)
                x = np.flip(counts[sort_index])
                max_rank_mouse = len(x)
                if rep == n_ensemble - 1:
                    ax_r2.step(range(1, max_rank_mouse+1), x/largest, color = my_colors[0], alpha = .5, lw = 0.5)
                
                if max_rank_mouse>max_rank:
                    x = x[:max_rank]
                else:
                    x = np.pad(x, (0, max_rank - max_rank_mouse), mode='constant')
            
                for k in range(max_rank):
                    if(x[k]>0):
                        counts_per_ranking[k]+=1
                        x_avg[k]+=x[k]/largest


        max_rank_eff = len(counts_per_ranking[counts_per_ranking>2])

        x_avg = x_avg[:max_rank_eff]/counts_per_ranking[:max_rank_eff]

        params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1))[:max_rank_fit], np.log(x_avg)[:max_rank_fit])
        # print(np.sqrt(pcov))
        slope = params[0]
        zetas.append(-slope)
        zeta = 3*3.5/(4.5*2.1)
        
    print(np.mean(zetas), int((np.mean(zetas)-zeta_min)*100))
    for j in range(len_mice):
        ax_r2.lines[-(j+1)].set_color(my_colors_alpha[int((np.mean(zetas)-zeta_min)*100)])

    ax_r2.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int((np.mean(zetas)-zeta_min)*100)], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = markers_ph[i_ph], label = r'$%.2f$'%(np.mean(zetas)) + ' ; ' + ph)
    ax_r2.plot(np.arange(1, max_rank_eff + 1), np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int((np.mean(zetas)-zeta_min)*100)], alpha = .8, lw = 3)

    ax_zeta2.hist(zetas, bins = np.linspace(0.2, 1.6, 30), alpha = .7, label = r'$\mathrm{GC+m}$', color = my_colors_alpha[int((np.mean(zetas)-zeta_min)*100)], density = True, histtype = 'stepfilled', edgecolor = 'k')

my_plot_layout(ax =ax_r2, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r2.set_ylim(bottom = 2e-2, top = 1.1)
ax_r2.set_xlim(right = 5e1)
ax_r2.legend(title = r'$\mathrm{exponent }\,\zeta$', fontsize = 24, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r2.savefig(output_plot + '/ranking_B_cells2.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta2, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
# ax_zeta2.set_ylim(bottom = 2e-2, top = 1.1)
ax_zeta2.set_xlim(left = 0.15, right = 1.6)
# ax_zeta2.legend(title = r'$\mathrm{sub-pop}$', fontsize = 22, title_fontsize = 30, loc = (1, 0))
fig_zeta2.savefig(output_plot + '/zetas.pdf', transparent=.5)
