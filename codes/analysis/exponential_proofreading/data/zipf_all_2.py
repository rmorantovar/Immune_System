import sys
sys.path.append('../../../library/')
from funcs import*
plt.rcParams['text.usetex'] = True

project = 'exponential_proofreading'
subproject = 'data'
experiment = 0
root_dir = f"/Users/robertomorantovar/Dropbox/Research_files/Immune_system/{project}/{subproject}/mesin2020"
output_plot = '/Users/robertomorantovar/Dropbox/_Documents/Research/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject+'/mesin2020/zipf/all'
os.makedirs(output_plot, exist_ok=True)

# Parameters
gs = [2]  # Number of Poisson processes
mu = 1.0  # Poisson rate
T = 15  # Total simulation time
theta = 1.5  # Values of theta to compare
gamma = 0.4
my_colors = [my_blue2, my_purple, my_purple, my_purple, my_cyan]
my_colors2 = [my_purple, my_purple, my_purple, my_cyan, my_purple, my_blue2]
my_colors3 = [my_blue2, my_purple, my_purple, my_blue, my_blue2, my_purple, my_purple, my_blue2, my_blue2, my_blue2]
my_colors4 = [my_blue2, my_purple, my_purple, my_blue, my_blue2, my_purple, my_purple, my_blue2, my_blue2, my_blue2]
alpha = 1e-10
depth = 6
anti_mut_epi = 5/4
n_ensemble = 100

color_vals = np.linspace(0, 2, 200)
cmap = plt.get_cmap('autumn_r')
my_colors_alpha = [cmap(val) for val in color_vals] 
my_colors_alpha_proposal = [my_blue, my_red]

def model(x, m):
    return m * x 

# fig_r, ax_r = plt.subplots(figsize=(10*1.62,8), gridspec_kw={'left':0.12, 'right':.8, 'bottom':.15, 'top': 0.94})
fig_r, ax_r = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 4e1)
# ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_0.pdf', transparent=.5)

fig_zeta, ax_zeta = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
fig_scaling1, ax_scaling1 = plt.subplots(figsize = (8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
fig_scaling2, ax_scaling2 = plt.subplots(figsize = (8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})

scaling_dict = defaultdict(list)

#------------ Experiment 1 (Figure 1D) ------------
print("Experiment 1")
data_primary = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Photoactivation CGG', header = 1)
data_primary = data_primary[(data_primary['Figure']==1)]
data_primary_grouped = data_primary.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
# data_primary_grouped = data_primary.groupby(['Mouse', 'AA JUNCTION']).size().reset_index(name='count')
# data_primary_grouped = data_primary.groupby(['Mouse', 'Sequence']).size().reset_index(name='count')
mice = data_primary_grouped['Mouse'].unique()
# phenotypes = data_primary_grouped['Phenotype'].unique()

# for i_ph, ph in enumerate(phenotypes):
# 	max_rank = max_ranks[i_ph]
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
	min_max_rank_mouse = max_rank
	max_max_rank_mouse = 0
	for mouse in mice_rep:
		data_mouse = data_primary_grouped[data_primary_grouped['Mouse']==mouse]
		counts = data_mouse['count'].to_numpy()
		N = np.sum(counts)
		S_i = -np.sum((counts/N)*np.log((counts/N)))
		sort_index = counts.argsort()
		largest = np.max(counts)
		x = np.flip(counts[sort_index])
		max_rank_mouse = len(x)
		if rep == n_ensemble - 1:
			ax_r.step(range(1, max_rank_mouse+1), x/largest, color = my_red, alpha = .5, lw = 0.5)
			params_mouse, pcov_mouse = curve_fit(model, np.log(range(1, len(x)+1)), np.log(x/largest))
			zetas_mice.append(-params_mouse[0])
			scaling_dict['experiment'].append('1')
			scaling_dict['response'].append('naive')
			scaling_dict['phenotype'].append('GC')
			scaling_dict['N1'].append(np.max(counts))
			scaling_dict['L_act'].append(np.size(counts))
			scaling_dict['barN'].append(np.sum(counts))
			scaling_dict['S'].append(S_i)
		if len(x)>max_rank:
			x = x[:max_rank]
		else:
			x = np.pad(x, (0, max_rank - len(x)), mode='constant')

		if max_rank_mouse < min_max_rank_mouse:
			min_max_rank_mouse = max_rank_mouse
		if max_rank_mouse > max_max_rank_mouse:
			max_max_rank_mouse = max_rank_mouse
		
		for k in range(max_rank):
			if(x[k]>0):
				counts_per_ranking[k]+=1
				x_avg[k]+=x[k]/largest

	max_rank_eff = len(counts_per_ranking[counts_per_ranking>2])

	x_avg = x_avg[:max_rank_eff]/counts_per_ranking[:max_rank_eff]

	params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1))[:max_rank_fit], np.log(x_avg)[:max_rank_fit])
	slope = params[0]

	zetas.append(-slope)

print(np.mean(zetas), np.std(zetas))

for j in range(len(mice)):
	ax_r.lines[-(j+1)].set_color(my_red)

ax_r.plot(range(1, max_rank_eff+1), x_avg, color = my_red, markerfacecolor="None", ms = 18, alpha = 1, ls = '', marker = '*', label = r'$%.2f$'%(np.mean(zetas)) + ' ; ' + 'GC')
ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_red, alpha = .8, lw = 3)

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

ax_scaling1.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), (np.linspace(1, 400, 100)*1.0)**(np.mean(zetas)), color = my_red, linestyle = '--')
ax_scaling2.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), np.linspace(1, 400, 100), color = my_red, linestyle = '--')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
# ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_1.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_zeta.set_xticks([0, 1], ['', ''])
ax_zeta.set_ylabel(r'$\zeta$', fontsize = 30)
ax_zeta.tick_params(axis='y', labelsize=30)
ax_zeta.tick_params(axis='x', labelsize=14, rotation = 45)
# ax_zeta.legend(title = r'$\mathrm{sub-pop}$', fontsize = 30, title_fontsize = 30, loc = (1, 0))
fig_zeta.savefig(output_plot + '/zetas_1.pdf', transparent=.5)


#------------ Experiment 2 (Figure 4A) ------------
print("Experiment 2 (Figure 4A)")
data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
data_recall = data_recall[(data_recall['Figure']=='4A')]
data_recall_grouped = data_recall.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_recall.groupby(['Mouse', 'AA JUNCTION']).size().reset_index(name='count')
# data_recall_grouped = data_recall.groupby(['Mouse', 'Sequence']).size().reset_index(name='count')
mice = data_recall_grouped['Mouse'].unique()

max_rank = 100
zetas = []
zetas_mice = []
for rep in tqdm(range(n_ensemble)):

	if rep == n_ensemble - 1:
		mice_rep = mice
	else:
		mice_rep = np.random.choice(mice, len(mice), replace = True)
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	# data_ph = data_recall_grouped[(data_recall_grouped['Phenotype']==ph)]
	min_max_rank_mouse = max_rank
	max_max_rank_mouse = 0
	for mouse in mice_rep:
		data_mouse = data_recall_grouped[data_recall_grouped['Mouse']==mouse]
		counts = data_mouse['count'].to_numpy()
		N = np.sum(counts)
		S_i = -np.sum((counts/N)*np.log((counts/N)))
		sort_index = counts.argsort()
		largest = np.max(counts)
		x = np.flip(counts[sort_index])
		max_rank_mouse = len(x)
		if rep == n_ensemble - 1:
			ax_r.step(range(1, max_rank_mouse+1), x/largest, color = my_blue, alpha = .5, lw = 0.5)
			params_mouse, pcov_mouse = curve_fit(model, np.log(range(1, len(x)+1)), np.log(x/largest))
			zetas_mice.append(-params_mouse[0])
			scaling_dict['experiment'].append('2')
			scaling_dict['response'].append('recall')
			scaling_dict['phenotype'].append('GC + fm')
			scaling_dict['N1'].append(np.max(counts))
			scaling_dict['L_act'].append(np.size(counts))
			scaling_dict['barN'].append(np.sum(counts))
			scaling_dict['S'].append(S_i)
    
		if max_rank_mouse>max_rank:
			x = x[:max_rank]
		else:
			x = np.pad(x, (0, max_rank - max_rank_mouse), mode='constant')
		 
		if max_rank_mouse < min_max_rank_mouse:
			min_max_rank_mouse = max_rank_mouse
		if max_rank_mouse > max_max_rank_mouse:
			max_max_rank_mouse = max_rank_mouse

		for k in range(max_rank):
			if(x[k]>0):
				counts_per_ranking[k]+=1
				x_avg[k]+=x[k]/largest

	max_rank_eff = len(counts_per_ranking[counts_per_ranking>2])

	x_avg = x_avg[:max_rank_eff]/counts_per_ranking[:max_rank_eff]

	params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1))[:max_rank_fit], np.log(x_avg)[:max_rank_fit])
	slope = params[0]
	zetas.append(-slope)
	zeta = 3*3.5/(4.5*2.1)
		
print(np.mean(zetas), np.std(zetas))

for j in range(len(mice)):
	ax_r.lines[-(j+1)].set_color(my_blue)

ax_r.plot(range(1, max_rank_eff+1), x_avg, color = my_blue, markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = 'o', label = r'$%.2f$'%(np.mean(zetas)) + ' ; ' + 'GC + fm')
ax_r.plot(np.arange(1, max_rank_eff + 1), np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_blue, alpha = .8, lw = 3)

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

ax_scaling1.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), (np.linspace(1, 400, 100)*1.0)**(np.mean(zetas)), color = my_blue, linestyle = '--')
ax_scaling2.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), np.linspace(1, 400, 100), color = my_blue, linestyle = '--')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
# ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_2.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_zeta.set_xticks([0, 1], ['', ''])
ax_zeta.set_ylabel(r'$\zeta$', fontsize = 30)
ax_zeta.tick_params(axis='y', labelsize=30)
ax_zeta.tick_params(axis='x', labelsize=14, rotation = 45)
fig_zeta.savefig(output_plot + '/zetas_2.pdf', transparent=.5)

#------------ Experiment 3 (Figure 4C) ------------
print('Experiment 3 (Figure 4C)')
data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
data_recall = data_recall[(data_recall['Figure']=='4C-H')]
data_recall_grouped = data_recall.groupby(['Mouse', 'Phenotype', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_recall.groupby(['Mouse', 'Phenotype', 'AA JUNCTION']).size().reset_index(name='count')
# data_recall_grouped = data_recall.groupby(['Mouse', 'Phenotype', 'Sequence']).size().reset_index(name='count')
mice = data_recall_grouped['Mouse'].unique()
phenotypes = data_recall_grouped['Phenotype'].unique()
print(phenotypes)

for i_ph, ph in enumerate(phenotypes):
	max_rank = 100
	zetas = []
	zetas_mice = []
	for rep in tqdm(range(n_ensemble)):

		if rep == n_ensemble - 1:
			mice_rep = mice
		else:
			mice_rep = np.random.choice(mice, len(mice), replace = True)

		x_avg = np.zeros(max_rank)
		counts_per_ranking = np.zeros(max_rank)
		data_ph = data_recall_grouped[(data_recall_grouped['Phenotype']==ph)]
		for mouse in mice_rep:
			data_mouse = data_ph[data_ph['Mouse']==mouse]
			counts = data_mouse['count'].to_numpy()
			N = np.sum(counts)
			S_i = -np.sum((counts/N)*np.log((counts/N)))
			sort_index = counts.argsort()
			largest = np.max(counts)
			x = np.flip(counts[sort_index])
			max_rank_mouse = len(x)
			if rep == n_ensemble - 1:
				ax_r.step(range(1, max_rank_mouse+1), x/largest, color = my_blue, alpha = .5, lw = 0.5)
				params_mouse, pcov_mouse = curve_fit(model, np.log(range(1, len(x)+1)), np.log(x/largest))
				zetas_mice.append(-params_mouse[0])
				scaling_dict['experiment'].append('3')
				scaling_dict['response'].append('recall')
				scaling_dict['phenotype'].append(ph)
				scaling_dict['N1'].append(np.max(counts))
				scaling_dict['L_act'].append(np.size(counts))
				scaling_dict['barN'].append(np.sum(counts))
				scaling_dict['S'].append(S_i)
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
		slope = params[0]
		zetas.append(-slope)
		zeta = 3*3.5/(4.5*2.1)
					
	for j in range(len(mice)):
		ax_r.lines[-(j+1)].set_color(my_blue)

	ax_r.plot(range(1, max_rank_eff+1), x_avg, color = my_blue, markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = '^', label = r'$%.2f$'%(np.mean(zetas))+ ' ; ' + ph)
	ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_blue, alpha = .8, lw = 3)

	parts = ax_zeta.violinplot([zetas], positions = [i_ph+2], showmeans = True, showextrema = False)
	for i, body in enumerate(parts['bodies']):
		body.set_facecolor(my_blue)
		body.set_edgecolor('black')
		body.set_alpha(0.5)
	parts['cmeans'].set_color(my_blue)
	# parts['cbars'].set_color('black')
	# parts['cmins'].set_color('black')
	# parts['cmaxes'].set_color('black')
	parts['cmeans'].set_color(my_blue)
	ax_zeta.scatter(np.random.normal(i_ph+2, 0.04, len(zetas_mice)), (zetas_mice), color = my_blue, edgecolor = 'k', s = 80, alpha = .5)
	ax_zeta.scatter(i_ph+2, (np.mean(zetas_mice)), color = my_blue, edgecolor = 'k', s = 150)

	ax_scaling1.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), (np.linspace(1, 400, 100)*1.0)**(np.mean(zetas)), color = my_blue, linestyle = '--')
	ax_scaling2.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), np.linspace(1, 400, 100), color = my_blue, linestyle = '--')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
# ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_3.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_zeta.set_xticks([0, 1, 2, 3], ['', '', 'GC + fm', 'PB + fm'])
ax_zeta.set_ylabel(r'$\zeta$', fontsize = 30)
ax_zeta.tick_params(axis='y', labelsize=30)
ax_zeta.tick_params(axis='x', labelsize=14, rotation = 45)
fig_zeta.savefig(output_plot + '/zetas_3.pdf', transparent=.5)

#------------ Experiment 3 (Figure 4C) ------------ 2
print('Experiment 3 (Figure 4C) - 2')
data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
data_recall = data_recall[(data_recall['Figure']=='4C-H')]
data_recall_grouped = data_recall.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_recall.groupby(['Mouse', 'AA JUNCTION']).size().reset_index(name='count')
# data_recall_grouped = data_recall.groupby(['Mouse', 'Phenotype', 'Sequence']).size().reset_index(name='count')
mice = data_recall_grouped['Mouse'].unique()

max_rank = 100
zetas = []
zetas_mice = []
for rep in tqdm(range(n_ensemble)):
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	if rep == n_ensemble - 1:
		mice_rep = mice
	else:
		mice_rep = np.random.choice(mice, len(mice), replace = True)
	
	for mouse in mice_rep:
		data_mouse = data_recall_grouped[data_recall_grouped['Mouse']==mouse]
		counts = data_mouse['count'].to_numpy()
		N = np.sum(counts)
		S_i = -np.sum((counts/N)*np.log((counts/N)))
		sort_index = counts.argsort()
		largest = np.max(counts)
		x = np.flip(counts[sort_index])
		max_rank_mouse = len(x)
		if rep == n_ensemble - 1:
			ax_r.step(range(1, max_rank_mouse+1), x/largest, color = my_colors2[0], alpha = .5, lw = 0.5)
			params_mouse, pcov_mouse = curve_fit(model, np.log(range(1, len(x)+1)), np.log(x/largest))
			zetas_mice.append(-params_mouse[0])
			scaling_dict['experiment'].append('3')
			scaling_dict['response'].append('recall')
			scaling_dict['phenotype'].append('combined')
			scaling_dict['N1'].append(np.max(counts))
			scaling_dict['L_act'].append(np.size(counts))
			scaling_dict['barN'].append(np.sum(counts))
			scaling_dict['S'].append(S_i)

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
	slope = params[0]
	zetas.append(-slope)
	zeta = 3*3.5/(4.5*2.1)

for j in range(len(mice)):
	ax_r.lines[-(j+1)].set_color(my_colors_alpha[int(np.mean(zetas)*100)])

print(np.mean(zetas), np.std(zetas))
# ax_r.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int(np.mean(zetas)*100)], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = '.', label = r'$%.2f$'%(np.mean(zetas))+ ' ; ' + 'combined')
# ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int(np.mean(zetas)*100)], alpha = .8, lw = 3)

parts = ax_zeta.violinplot([zetas], positions = [4], showmeans = True, showextrema = False)
for i, body in enumerate(parts['bodies']):
	body.set_facecolor(my_blue)
	body.set_edgecolor('black')
	body.set_alpha(0.5)
# parts['cbars'].set_color('black')
# parts['cmins'].set_color('black')
# parts['cmaxes'].set_color('black')
parts['cmeans'].set_color(my_blue)
ax_zeta.scatter(np.random.normal(4, 0.04, len(zetas_mice)), (zetas_mice), color = my_blue, edgecolor = 'k', s = 80, alpha = .5)
ax_zeta.scatter(4, (np.mean(zetas_mice)), color = my_blue, edgecolor = 'k', s = 150)

ax_scaling1.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), (np.linspace(1, 400, 100)*1.0)**(np.mean(zetas)), color = my_blue2, linestyle = '--')
ax_scaling2.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), np.linspace(1, 400, 100), color = my_blue2, linestyle = '--')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
# ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_3b.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_zeta.set_xticks([0, 1, 2, 3, 4], ['', '', 'GC + fm', 'GC + fm + recall', 'combined'])
ax_zeta.set_ylabel(r'$\zeta$', fontsize = 30)
ax_zeta.tick_params(axis='y', labelsize=30)
ax_zeta.tick_params(axis='x', labelsize=14, rotation = 45)
fig_zeta.savefig(output_plot + '/zetas_3b.pdf', transparent=.5)

#------------ Experiment 4 (Figure 5) ------------
colors_ph = [my_red, my_blue, my_red, my_blue, my_red]
colors_ph = [my_red, my_blue, my_purple, my_blue]
print('Experiment 4 (Figure 5)')

data_infection = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Influenza_IGH', header = 2)
data_infection = data_infection[(data_infection['Sort']=='GC')  | (data_infection['Sort']=='GC + fm') | (data_infection['Sort']=='PB + fm') | (data_infection['Sort']=='M')]
data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort', 'CDR3:']).size().reset_index(name='count')
# data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort2', 'Sequence']).size().reset_index(name='count')

mice = data_recall_grouped['Experiment / Mouse'].unique()
print(mice)
phenotypes = data_recall_grouped['Sort'].unique()
print(phenotypes)

max_ranks = [100, 100, 100, 100, 100]
# max_ranks = [10, 10, 10, 10, 10]
for i_ph, ph in enumerate(phenotypes):
	max_rank = 100
	zetas = []
	zetas_mice = []

	for rep in tqdm(range(n_ensemble)):

		if rep == n_ensemble - 1:
			mice_rep = mice
		else:
			mice_rep = np.random.choice(mice, len(mice), replace = True)

		x_avg = np.zeros(max_rank)
		counts_per_ranking = np.zeros(max_rank)
		data_ph = data_recall_grouped[(data_recall_grouped['Sort']==ph)]
		min_max_rank_mouse = max_rank
		max_max_rank_mouse = 0
		for mouse in mice_rep:
			data_mouse = data_ph[data_ph['Experiment / Mouse']==mouse]
			counts = data_mouse['count'].to_numpy()
			N = np.sum(counts)
			S_i = -np.sum((counts/N)*np.log((counts/N)))
			sort_index = counts.argsort()
			largest = np.max(counts)
			x = np.flip(counts[sort_index])
			max_rank_mouse = len(x)
			if rep == n_ensemble - 1:
				ax_r.step(range(1, max_rank_mouse+1), x/largest, color = colors_ph[i_ph], alpha = .5, lw = 0.5)
				params_mouse, pcov_mouse = curve_fit(model, np.log(range(1, len(x)+1)), np.log(x/largest))
				zetas_mice.append(-params_mouse[0])
				scaling_dict['experiment'].append('4')
				scaling_dict['response'].append('recall')
				scaling_dict['phenotype'].append(ph)
				scaling_dict['N1'].append(np.max(counts))
				scaling_dict['L_act'].append(np.size(counts))
				scaling_dict['barN'].append(np.sum(counts))
				scaling_dict['S'].append(S_i)

			if len(x)>max_rank:
				x = x[:max_rank]
			else:
				x = np.pad(x, (0, max_rank - len(x)), mode='constant')
			
			if max_rank_mouse < min_max_rank_mouse:
				min_max_rank_mouse = max_rank_mouse
			if max_rank_mouse > max_max_rank_mouse:
				max_max_rank_mouse = max_rank_mouse

			for k in range(max_rank):
				if(x[k]>0):
					counts_per_ranking[k]+=1
					x_avg[k]+=x[k]/largest

		max_rank_eff = len(counts_per_ranking[counts_per_ranking>2])
		x_avg = x_avg[:max_rank_eff]/counts_per_ranking[:max_rank_eff]
    
		params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1))[:max_rank_fit], np.log(x_avg[:])[:max_rank_fit])
		slope = params[0]

		zetas.append(-slope)

		zeta = 3*3.5/(4.5*2.1)
		

	for j in range(len(mice)):
		ax_r.lines[-(j+1)].set_color(colors_ph[i_ph])

	print(np.mean(zetas), np.std(zetas))
	
	ax_r.plot(range(1, max_rank_eff+1), x_avg, color = colors_ph[i_ph], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = 'D', label = r'$%.2f$'%(np.mean(zetas)) + ' ; ' + ph)
	ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = colors_ph[i_ph], alpha = .8, lw = 3)
	
	parts = ax_zeta.violinplot([zetas], positions = [i_ph+5], showmeans = True, showextrema = False)
	for i, body in enumerate(parts['bodies']):
		body.set_facecolor(colors_ph[i_ph])
		body.set_edgecolor('black')
		body.set_alpha(0.5)
	parts['cmeans'].set_color(colors_ph[i_ph])
	# parts['cbars'].set_color('black')
	# parts['cmins'].set_color('black')
	# parts['cmaxes'].set_color('black')
	parts['cmeans'].set_color(colors_ph[i_ph])
	ax_zeta.scatter(np.random.normal(i_ph+5, 0.04, len(zetas_mice)), (zetas_mice), color = colors_ph[i_ph], edgecolor = 'k', s = 80, alpha = .5)
	ax_zeta.scatter(i_ph+5, (np.mean(zetas_mice)), color = colors_ph[i_ph], edgecolor = 'k', s = 150)

	ax_scaling1.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), (np.linspace(1, 400, 100)*1.0)**(np.mean(zetas)), color = colors_ph[i_ph], linestyle = '--')
	ax_scaling2.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), np.linspace(1, 400, 100), color = colors_ph[i_ph], linestyle = '--')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
# ax_r.legend(title = r'$\zeta$', fontsize = 24, title_fontsize = 28, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_4.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 20, x_fontsize=30, y_fontsize=30 )
ax_zeta.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8], ['exp1/naive', 'exp2/recall', 'exp3/recall', 'exp3/recall', 'exp3/recall', 'exp4/naive', 'exp4/recall', 'exp4/M', 'exp4/recall'])
ax_zeta.set_ylabel(r'$\zeta$', fontsize = 30)
ax_zeta.tick_params(axis='y', labelsize=30)
ax_zeta.tick_params(axis='x', labelsize=14, rotation = 45)
fig_zeta.savefig(output_plot + '/zetas_4.pdf', transparent=.5)


#------------ Experiment 4 (Figure 5) ------------ 2
colors_ph = [my_red, my_blue, my_red, my_blue, my_red]
colors_ph = [my_blue, my_blue, my_blue, my_blue, my_red]
print('Experiment 4 (Figure 5) - 2')

data_infection = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Influenza_IGH', header = 2)
data_infection = data_infection[(data_infection['Sort2']=='fm')]
data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort2', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort2', 'CDR3:']).size().reset_index(name='count')
# data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort2', 'Sequence']).size().reset_index(name='count')

mice = data_recall_grouped['Experiment / Mouse'].unique()
print(mice)
phenotypes = data_recall_grouped['Sort2'].unique()
print(phenotypes)

max_ranks = [100, 100, 100, 100, 100]
# max_ranks = [10, 10, 10, 10, 10]
for i_ph, ph in enumerate(phenotypes):
	max_rank = 100
	zetas = []
	zetas_mice = []

	for rep in tqdm(range(n_ensemble)):

		if rep == n_ensemble - 1:
			mice_rep = mice
		else:
			mice_rep = np.random.choice(mice, len(mice), replace = True)

		x_avg = np.zeros(max_rank)
		counts_per_ranking = np.zeros(max_rank)
		data_ph = data_recall_grouped[(data_recall_grouped['Sort2']==ph)]
		min_max_rank_mouse = max_rank
		max_max_rank_mouse = 0
		for mouse in mice_rep:
			data_mouse = data_ph[data_ph['Experiment / Mouse']==mouse]
			counts = data_mouse['count'].to_numpy()
			N = np.sum(counts)
			S_i = -np.sum((counts/N)*np.log((counts/N)))
			sort_index = counts.argsort()
			largest = np.max(counts)
			x = np.flip(counts[sort_index])
			max_rank_mouse = len(x)
			if rep == n_ensemble - 1:
				ax_r.step(range(1, max_rank_mouse+1), x/largest, color = colors_ph[i_ph], alpha = .5, lw = 0.5)
				params_mouse, pcov_mouse = curve_fit(model, np.log(range(1, len(x)+1)), np.log(x/largest))
				zetas_mice.append(-params_mouse[0])
				scaling_dict['experiment'].append('4')
				scaling_dict['response'].append('recall')
				scaling_dict['phenotype'].append('combined')
				scaling_dict['N1'].append(np.max(counts))
				scaling_dict['L_act'].append(np.size(counts))
				scaling_dict['barN'].append(np.sum(counts))
				scaling_dict['S'].append(S_i)

			if len(x)>max_rank:
				x = x[:max_rank]
			else:
				x = np.pad(x, (0, max_rank - len(x)), mode='constant')
			
			if max_rank_mouse < min_max_rank_mouse:
				min_max_rank_mouse = max_rank_mouse
			if max_rank_mouse > max_max_rank_mouse:
				max_max_rank_mouse = max_rank_mouse

			for k in range(max_rank):
				if(x[k]>0):
					counts_per_ranking[k]+=1
					x_avg[k]+=x[k]/largest

		max_rank_eff = len(counts_per_ranking[counts_per_ranking>2])
		x_avg = x_avg[:max_rank_eff]/counts_per_ranking[:max_rank_eff]
    
		params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1))[:max_rank_fit], np.log(x_avg[:])[:max_rank_fit])
		slope = params[0]

		zetas.append(-slope)

		zeta = 3*3.5/(4.5*2.1)
		

	for j in range(len(mice)):
		ax_r.lines[-(j+1)].set_color(colors_ph[i_ph])

	print(np.mean(zetas), np.std(zetas))
	
	# ax_r.plot(range(1, max_rank_eff+1), x_avg, color = colors_ph[i_ph], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = 'D', label = r'$%.2f$'%(np.mean(zetas)) + ' ; ' + ph)
	# ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = colors_ph[i_ph], alpha = .8, lw = 3)
	
	parts = ax_zeta.violinplot([zetas], positions = [i_ph+9], showmeans = True, showextrema = False)
	for i, body in enumerate(parts['bodies']):
		body.set_facecolor(colors_ph[i_ph])
		body.set_edgecolor('black')
		body.set_alpha(0.5)
	parts['cmeans'].set_color(colors_ph[i_ph])
	# parts['cbars'].set_color('black')
	# parts['cmins'].set_color('black')
	# parts['cmaxes'].set_color('black')
	parts['cmeans'].set_color(colors_ph[i_ph])
	ax_zeta.scatter(np.random.normal(i_ph+9, 0.04, len(zetas_mice)), (zetas_mice), color = colors_ph[i_ph], edgecolor = 'k', s = 80, alpha = .5)
	ax_zeta.scatter(i_ph+9, (np.mean(zetas_mice)), color = colors_ph[i_ph], edgecolor = 'k', s = 150)

	ax_scaling1.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), (np.linspace(1, 400, 100)*1.0)**(np.mean(zetas)), color = my_blue2, linestyle = '--')
	ax_scaling1.plot((np.linspace(1, 400, 100))**(np.mean(zetas))/(1-np.mean(zetas))*((np.linspace(1, 400, 100))**(1-np.mean(zetas))-1), (np.linspace(1, 400, 100)*1.0)**(np.mean(zetas)), color = my_blue2, linestyle = '--')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
# ax_r.legend(title = r'$\zeta$', fontsize = 24, title_fontsize = 28, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_4b.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 20, x_fontsize=30, y_fontsize=30 )
ax_zeta.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], ['exp1/naive', 'exp2/recall', 'exp3/recall', 'exp3/recall', 'exp3/recall', 'exp4/naive', 'exp4/recall', 'exp4/M', 'exp4/recall', 'exp4/recall'])
ax_zeta.set_ylabel(r'$\zeta$', fontsize = 30)
ax_zeta.tick_params(axis='y', labelsize=30)
ax_zeta.tick_params(axis='x', labelsize=14, rotation = 45)
fig_zeta.savefig(output_plot + '/zetas_4b.pdf', transparent=.5)


# #------------ Scaling plot ------------

scaling_results = pd.DataFrame(scaling_dict)

sns.scatterplot(data = scaling_results, x = 'barN', y = 'N1', hue = 'phenotype', style = 'experiment', ax = ax_scaling1, s = 100, palette = [my_red, my_blue, my_blue, my_blue2, my_purple, 'darkorange', 'purple', 'brown', 'pink'], edgecolors = 'black', alpha = 0.8)
# ax_scaling1.plot(np.linspace(1, 130, 100), (np.linspace(1, 130, 100)*1.0)**(1), color = 'k', linestyle = '--', alpha = .8)
ax_scaling1.set_xlabel(r'$N_B^{\mathrm{tot}}$', fontsize = 30)
ax_scaling1.set_ylabel(r'$N_1$', fontsize = 30)
ax_scaling1.set_ylim(bottom=0, top = 80)
ax_scaling1.set_xlim(left=1, right = 400)
ax_scaling1.tick_params(labelsize = 20)
# ax_scaling1.set_yscale('log')
# ax_scaling1.set_xscale('log')
ax_scaling1.legend(title = 'Response', title_fontsize = 20, fontsize = 15, loc = 4)
fig_scaling1.savefig(output_plot + '/size_scaling_1_linear.pdf', bbox_inches = 'tight')

sns.scatterplot(data = scaling_results, x = 'barN', y = 'L_act', hue = 'phenotype', style = 'experiment', ax = ax_scaling2, s = 100, palette = [my_red, my_blue, my_blue, my_blue2, my_purple, 'darkorange', 'purple', 'brown', 'pink'], edgecolors = 'black', alpha = 0.8)
# ax_scaling2.plot(np.linspace(1, 130, 100), (np.linspace(1, 130, 100)*1.0)**(1), color = 'k', linestyle = '--', alpha = .8)
ax_scaling2.set_xlabel(r'$N_B^{\mathrm{tot}}$', fontsize = 30)
ax_scaling2.set_ylabel(r'$L_{act}$', fontsize = 30)
ax_scaling2.set_ylim(bottom=0, top = 200)
ax_scaling2.set_xlim(left=1, right = 400)
ax_scaling2.tick_params(labelsize = 20)
# ax_scaling2.set_yscale('log')
# ax_scaling2.set_xscale('log')
ax_scaling2.legend(title = 'Response', title_fontsize = 20, fontsize = 15, loc = 2)
fig_scaling2.savefig(output_plot + '/size_scaling_2_linear.pdf', bbox_inches = 'tight')


# #------------ Experiment 4 (day 9) ------------

# data_infection = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Influenza9', header = 2)
# # data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort', 'CDR3:']).size().reset_index(name='count')
# mice = data_recall_grouped['Experiment / Mouse'].unique()
# phenotypes = data_recall_grouped['Sort'].unique()

# max_ranks = [7, 7, 7, 20, 7]
# for i_ph, ph in enumerate(phenotypes):
# 	max_rank = max_ranks[i_ph]
# 	x_avg = np.zeros(max_rank)
# 	counts_per_ranking = np.zeros(max_rank)
# 	data_ph = data_recall_grouped[(data_recall_grouped['Sort']==ph)]
# 	for mouse in mice:
# 		data_mouse = data_ph[data_ph['Experiment / Mouse']==mouse]
# 		# CDR3, counts = np.unique(np.array((list(data_mouse['CDR3:']))), return_counts = True)
# 		counts = data_mouse['count'].to_numpy()
# 		N = np.sum(counts)

# 		if(N>0):
# 			sort_index = counts.argsort()
# 			largest = np.max(counts)
# 			x = np.flip(counts[sort_index])
# 			if len(x)>max_rank:
# 				x = x[:max_rank]
# 			else:
# 				x = np.pad(x, (0, max_rank - len(x)), mode='constant')
# 			ax_r.step(range(1, max_rank+1), x/largest, color = my_colors4[i_ph], alpha = .5, lw = 0.5)
# 			for k in range(max_rank):
# 				if(x[k]>0):
# 					counts_per_ranking[k]+=1
# 					x_avg[k]+=x[k]/largest

# 	x_avg/=counts_per_ranking

# 	params, pcov = curve_fit(model, np.log(range(1, max_rank+1)), np.log(x_avg))
# 	slope = params[0]

# 	zeta = 3*3.5/(4.5*2.1)
# 	ax_r.plot(np.arange(1, 30), np.exp(0)*np.arange(1, 30)**(slope), color = my_colors4[i_ph], alpha = 1)
# 	ax_r.plot(range(1, max_rank+1), x_avg, color = my_colors4[i_ph], alpha = .8, ls = '', marker = '.', label = r'$%.2f\pm %.3f$'%(-slope, np.sqrt(pcov[0][0])) + ' ; ' + ph)

