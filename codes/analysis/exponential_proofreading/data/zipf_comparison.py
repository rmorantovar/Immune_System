import sys
sys.path.append('../../../library/')
from funcs import*
plt.rcParams['text.usetex'] = True

project = 'exponential_proofreading'
subproject = 'data'
experiment = 0
root_dir = f"/Users/robertomorantovar/Dropbox/Research_files/Immune_system/{project}/{subproject}/mesin2020"
output_plot = '/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/'+project+'/'+ subproject + '/mesin2020/zipf/comparison'
os.makedirs(output_plot, exist_ok=True)

# Parameters
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

def model(x, m):
    return m * x 

fig_r2, ax_r2 = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
fig_zeta, ax_zeta = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})

#------------ Experiment 1 (Figure 1D) ------------
print('Experiment 1 (Figure 1D)')
data_primary = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Photoactivation CGG', header = 1)
data_primary = data_primary[(data_primary['Figure']==1)]
data_primary_grouped = data_primary.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
# data_primary_grouped = data_primary.groupby(['Mouse', 'CDR3:']).size().reset_index(name='count')
# data_primary_grouped = data_primary.groupby(['Mouse', 'Sequence']).size().reset_index(name='count')
mice = data_primary_grouped['Mouse'].unique()
print(len(mice))
max_rank = 100
max_rank_fit = 20

zetas = []
zetas_mice = []
for rep in tqdm(range(n_ensemble)):
	if rep == n_ensemble - 1:
		mice_rep = mice
	else:
		mice_rep = np.random.choice(mice, len(mice), replace = True)
	
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	for mouse in mice_rep:
		data_mouse = data_primary_grouped[data_primary_grouped['Mouse']==mouse]
		# CDR3, counts = np.unique(np.array((list(data_mouse['CDR3:']))), return_counts = True)
		counts = data_mouse['count'].to_numpy()
		# print(counts)
		N = np.sum(counts)
		S_i = -np.sum((counts/N)*np.log((counts/N)))

		sort_index = counts.argsort()
		largest = np.max(counts)
		x = np.flip(counts[sort_index])
		max_rank_mouse = len(x)
		if rep == n_ensemble - 1:
			ax_r2.step(range(1, max_rank_mouse+1), x/largest, color = my_colors[0], alpha = .5, lw = 0.5)
			params_mouse, pcov_mouse = curve_fit(model, np.log(range(1, len(x)+1)), np.log(x/largest))
			zetas_mice.append(-params_mouse[0])
		if len(x)>max_rank:
			x = x[:max_rank]
		else:
			x = np.pad(x, (0, max_rank - len(x)), mode='constant')
		
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

print(zetas[-1])

ax_r2.plot(range(1, max_rank_eff+1), x_avg, color = my_red, markerfacecolor="None", ms = 18, alpha = 1, ls = '', marker = '*', label = r'$%.2f$'%(np.mean(zetas)) + ' ; ' + 'GC')

ax_r2.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_red, alpha = .8, lw = 3)

parts = ax_zeta.violinplot([zetas], positions = [0], showmeans = True, showextrema = False)
# ax_zeta.scatter(np.random.normal(0, 0.04, len(zetas)), (zetas), color = my_red, alpha = 0.5, edgecolor = 'None', s = 100)
for i, body in enumerate(parts['bodies']):
    body.set_facecolor(my_red)
    body.set_edgecolor('black')
    body.set_alpha(0.5)
# parts['cbars'].set_color('black')
# parts['cmins'].set_color('black')
# parts['cmaxes'].set_color('black')
parts['cmeans'].set_color(my_red)
ax_zeta.scatter(np.random.normal(0, 0.04, len(zetas_mice)), (zetas_mice), color = my_red, edgecolor = 'k', s = 80, alpha = .5)
ax_zeta.scatter(0, (np.mean(zetas_mice)), color = my_red, edgecolor = 'k', s = 150)


my_plot_layout(ax =ax_r2, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r2.set_ylim(bottom = 2e-2, top = 1.1)
ax_r2.set_xlim(right = 5e1)
ax_r2.legend(title = r'$\mathrm{exponent }\,\zeta$', fontsize = 24, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r2.savefig(output_plot + '/ranking_B_cells1.pdf', transparent=.5)


#------------ Experiment 2 and 3 (Figure 4A and 4C) ------------
print('Experiment 2 and 3 (Figure 4A and 4C)')
data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
# data_recall = data_recall[(data_recall['Phenotype']=='GC + fm')]

figures = ['4A', '4C-H']

max_rank = 100
zetas = []
zetas_mice = []

for rep in tqdm(range(n_ensemble)):
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	x_avg_fig = np.zeros(max_rank)
	counts_per_ranking_fig = np.zeros(max_rank)


	for i_fig, fig in enumerate(figures):
		data_recall_fig = data_recall[(data_recall['Figure']==fig)]
		data_recall_grouped = data_recall_fig.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
		# data_recall_grouped = data_recall_fig.groupby(['Mouse', 'CDR3:']).size().reset_index(name='count')
		# data_recall_grouped = data_recall_fig.groupby(['Mouse', 'Sequence']).size().reset_index(name='count')
		mice = data_recall_grouped['Mouse'].unique()
		if rep == n_ensemble - 1:
			mice_rep = mice
		else:
			mice_rep = np.random.choice(mice, len(mice), replace = True)

		for mouse in mice_rep:
			data_mouse = data_recall_grouped[data_recall_grouped['Mouse']==mouse]
			counts = data_mouse['count'].to_numpy()
			# print(counts)
			N = np.sum(counts)

			sort_index = counts.argsort()
			largest = np.max(counts)
			x = np.flip(counts[sort_index])
			max_rank_mouse = len(x)
			if rep == n_ensemble - 1:
				ax_r2.step(range(1, max_rank_mouse+1), x/largest, color = my_colors[0], alpha = .5, lw = 0.5)
				params_mouse, pcov_mouse = curve_fit(model, np.log(range(1, len(x)+1)), np.log(x/largest))
				zetas_mice.append(-params_mouse[0])
			
			if max_rank_mouse>max_rank:
				x = x[:max_rank]
			else:
				x = np.pad(x, (0, max_rank - max_rank_mouse), mode='constant')

			for k in range(max_rank):
				if(x[k]>0):
					counts_per_ranking[k]+=1
					x_avg[k]+=x[k]/largest
					if i_fig==1:
						counts_per_ranking_fig[k]+=1
						x_avg_fig[k]+=x[k]/largest


	max_rank_eff = len(counts_per_ranking[counts_per_ranking>2])
	max_rank_eff_fig = len(counts_per_ranking_fig[counts_per_ranking_fig>2])

	x_avg = x_avg[:max_rank_eff]/counts_per_ranking[:max_rank_eff]
	x_avg_fig = x_avg_fig[:max_rank_eff_fig]/counts_per_ranking_fig[:max_rank_eff_fig]

	params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1))[:max_rank_fit], np.log(x_avg)[:max_rank_fit])
	params_fig, pcov_fig = curve_fit(model, np.log(range(1, max_rank_eff_fig+1))[:max_rank_fit], np.log(x_avg_fig)[:max_rank_fit])
	
	slope = params[0]
	zetas.append(-slope)
	
print(zetas[-1])


ax_r2.plot(range(1, max_rank_eff+1), x_avg, color = my_blue, markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = 'o', label = r'$%.2f$'%(np.mean(zetas)) + ' ; ' + 'fm')
ax_r2.plot(np.arange(1, max_rank_eff + 1), np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_blue, alpha = .8, lw = 3)

parts = ax_zeta.violinplot([zetas], positions = [1], showmeans = True, showextrema = False)
# ax_zeta.scatter(np.random.normal(1, 0.04, len(zetas)), (zetas), color = my_blue, alpha = 0.5, edgecolor = 'None', s = 100)
for i, body in enumerate(parts['bodies']):
    body.set_facecolor(my_blue)
    body.set_edgecolor('black')
    body.set_alpha(0.5)
# parts['cbars'].set_color('black')
# parts['cmins'].set_color('black')
# parts['cmaxes'].set_color('black')
parts['cmeans'].set_color(my_blue)
ax_zeta.scatter(np.random.normal(1, 0.04, len(zetas_mice)), (zetas_mice), color = my_blue, edgecolor = 'k', s = 80, alpha = .5)
ax_zeta.scatter(1, (np.mean(zetas_mice)), color = my_blue, edgecolor = 'k', s = 150)

my_plot_layout(ax =ax_r2, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r2.set_ylim(bottom = 2e-2, top = 1.1)
ax_r2.set_xlim(right = 5e1)
ax_r2.legend(title = r'$\mathrm{exponent }\,\zeta$', fontsize = 24, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r2.savefig(output_plot + '/ranking_B_cells2.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
# ax_zeta.set_ylim(bottom = 2e-2, top = 1.1)
# ax_zeta.set_xlim(left = 0.15, right = 1.6)
ax_zeta.set_xticks([0, 1], ['', ''])
ax_zeta.set_ylabel(r'$\zeta$', fontsize = 30)
# ax_zeta.legend(title = r'$\mathrm{sub-pop}$', fontsize = 22, title_fontsize = 30, loc = (1, 0))
fig_zeta.savefig(output_plot + '/zetas.pdf', transparent=.5)

