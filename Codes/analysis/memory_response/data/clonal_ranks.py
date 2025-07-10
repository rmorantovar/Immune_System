import sys
sys.path.append('../../../my_lib/')
from funcs import*
plt.rcParams['text.usetex'] = True

project = 'memory_response'
subproject = 'data'
experiment = 0
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/mesin2020"
output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject+'/'+str(experiment) + '/mesin2020'
os.makedirs(output_plot, exist_ok=True)

# Parameters
n_ens = 10000
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
n_ensemble = 1000

color_vals = np.linspace(0, 2, 200)
cmap = plt.get_cmap('managua_r')
my_colors_alpha = [cmap(val) for val in color_vals] 

def model(x, m):
    return m * x 

# fig_r, ax_r = plt.subplots(figsize=(10*1.62,8), gridspec_kw={'left':0.12, 'right':.8, 'bottom':.15, 'top': 0.94})
fig_r, ax_r = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
fig_r2, ax_r2 = plt.subplots(figsize=(8,8), gridspec_kw={'left':0.15, 'right':.98, 'bottom':.15, 'top': 0.98})

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 4e1)
# ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_0.pdf', transparent=.5)

fig_zeta, ax_zeta = plt.subplots(figsize=(10*1.62,8), gridspec_kw={'left':0.12, 'right':.8, 'bottom':.15, 'top': 0.94})
fig_zeta2, ax_zeta2 = plt.subplots(figsize=(10*1.62,8), gridspec_kw={'left':0.12, 'right':.8, 'bottom':.15, 'top': 0.94})

#------------ Experiment 1 (Figure 1D) ------------
data_primary = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Photoactivation CGG', header = 1)
data_primary = data_primary[(data_primary['Figure']==1)]
data_primary_grouped = data_primary.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
# data_primary_grouped = data_primary.groupby(['Mouse', 'CDR3:']).size().reset_index(name='count')
mice = data_primary_grouped['Mouse'].unique()
# phenotypes = data_primary_grouped['Phenotype'].unique()

# for i_ph, ph in enumerate(phenotypes):
# 	max_rank = max_ranks[i_ph]
max_rank = 100
zetas = []
for rep in tqdm(range(n_ensemble)):
	if rep == n_ensemble - 1:
		mice_rep = mice
		# print(mice_rep)
	else:
		mice_rep = np.random.choice(mice, len(mice), replace = True)
	
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	min_max_rank_mouse = max_rank
	max_max_rank_mouse = 0
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
			ax_r.step(range(1, max_rank_mouse+1), x/largest, color = my_colors[0], alpha = .5, lw = 0.5)
			ax_r2.step(range(1, max_rank_mouse+1), x/largest, color = my_colors[0], alpha = .5, lw = 0.5)
		
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

	params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1)), np.log(x_avg))
	# print(np.sqrt(pcov))
	slope = params[0]

	zeta = 3*3.5/(4.5*2.1)
	zetas.append(-slope)
	# print(-slope)

print(len(mice))

for j in range(len(mice)):
	ax_r.lines[-(j+1)].set_color(my_colors_alpha[int(np.mean(zetas)*100)])
	ax_r2.lines[-(j+1)].set_color(my_colors_alpha[int(np.mean(zetas)*100)])

ax_r.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int(np.mean(zetas)*100)], markerfacecolor="None", ms = 18, alpha = 1, ls = '', marker = '*', label = r'$%.2f$'%(np.mean(zetas)))
ax_r2.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int(np.mean(zetas)*100)], markerfacecolor="None", ms = 18, alpha = 1, ls = '', marker = '*', label = r'$%.2f$'%(np.mean(zetas)))

ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int(np.mean(zetas)*100)], alpha = .8, lw = 3)
ax_r2.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int(np.mean(zetas)*100)], alpha = .8, lw = 3)

ax_zeta.hist(zetas, bins = np.linspace(0.2, 1.6, 20), alpha = .7, label = r'$\mathrm{GC}$', color = my_colors_alpha[int(np.mean(zetas)*100)], density = True, histtype = 'stepfilled', edgecolor = 'k')
ax_zeta2.hist(zetas, bins = np.linspace(0.2, 1.6, 20), alpha = .7, label = r'$\mathrm{GC}$', color = my_colors_alpha[int(np.mean(zetas)*100)], density = True, histtype = 'stepfilled', edgecolor = 'k')


my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_1.pdf', transparent=.5)


my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
# ax_zeta.set_ylim(bottom = 2e-2, top = 1.1)
ax_zeta.set_xlim(left = 0.2, right = 1.6)
ax_zeta.legend(title = r'$\mathrm{sub-pop}$', fontsize = 30, title_fontsize = 30, loc = (1, 0))
fig_zeta.savefig(output_plot + '/zetas_1.pdf', transparent=.5)


#------------ Experiment 2 (Figure 4A) ------------

data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
data_recall = data_recall[(data_recall['Figure']=='4A')]
data_recall_grouped = data_recall.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_recall.groupby(['Mouse', 'CDR3:']).size().reset_index(name='count')
mice = data_recall_grouped['Mouse'].unique()
# phenotypes = data_recall_grouped['Phenotype'].unique()

# for i_ph, ph in enumerate(phenotypes):
max_rank = 100
zetas = []
for rep in tqdm(range(n_ensemble)):

	if rep == n_ensemble - 1:
		mice_rep = mice
		# print(mice_rep)
	else:
		mice_rep = np.random.choice(mice, len(mice), replace = True)
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	# data_ph = data_recall_grouped[(data_recall_grouped['Phenotype']==ph)]
	min_max_rank_mouse = max_rank
	max_max_rank_mouse = 0
	for mouse in mice_rep:
		data_mouse = data_recall_grouped[data_recall_grouped['Mouse']==mouse]
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
			ax_r.step(range(1, max_rank_mouse+1), x/largest, color = my_colors2[0], alpha = .5, lw = 0.5)
			ax_r2.step(range(1, max_rank_mouse+1), x/largest, color = my_colors2[0], alpha = .5, lw = 0.5)
		
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

	max_rank_eff = len(counts_per_ranking[counts_per_ranking>3])

	x_avg = x_avg[:max_rank_eff]/counts_per_ranking[:max_rank_eff]

	params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1)), np.log(x_avg))
	# print(np.sqrt(pcov))
	slope = params[0]
	zetas.append(-slope)
	zeta = 3*3.5/(4.5*2.1)
		

print(len(mice))
for j in range(len(mice)):
	ax_r.lines[-(j+1)].set_color(my_colors_alpha[int(np.mean(zetas)*100)])

ax_r.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int(np.mean(zetas)*100)], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = 'o', label = r'$%.2f$'%(np.mean(zetas)))

ax_r.plot(np.arange(1, max_rank_eff + 1), np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int(np.mean(zetas)*100)], alpha = .8, lw = 3)
ax_zeta.hist(zetas, bins = np.linspace(0.2, 1.6, 20), alpha = .7, label = r'$\mathrm{GC+m}$', color = my_colors_alpha[int(np.mean(zetas)*100)], density = True, histtype = 'stepfilled', edgecolor = 'k')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_2.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
# ax_zeta.set_ylim(bottom = 2e-2, top = 1.1)
ax_zeta.set_xlim(left = 0.2, right = 1.6)
ax_zeta.legend(title = r'$\mathrm{sub-pop}$', fontsize = 30, title_fontsize = 30, loc = (1, 0))
fig_zeta.savefig(output_plot + '/zetas_2.pdf', transparent=.5)

#------------ Experiment 3 (Figure 4C) ------------

data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
data_recall = data_recall[(data_recall['Figure']=='4C-H')]
data_recall_grouped = data_recall.groupby(['Mouse', 'Phenotype', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_recall.groupby(['Mouse', 'Phenotype', 'CDR3:']).size().reset_index(name='count')
# print(data_recall_grouped)
mice = data_recall_grouped['Mouse'].unique()
phenotypes = data_recall_grouped['Phenotype'].unique()
print(phenotypes)

max_ranks = [100, 100]
for i_ph, ph in enumerate(phenotypes):
	max_rank = max_ranks[i_ph]
	zetas = []
	for rep in tqdm(range(n_ensemble)):

		if rep == n_ensemble - 1:
			mice_rep = mice
			# print(mice_rep)
		else:
			mice_rep = np.random.choice(mice, len(mice), replace = True)

		x_avg = np.zeros(max_rank)
		counts_per_ranking = np.zeros(max_rank)
		data_ph = data_recall_grouped[(data_recall_grouped['Phenotype']==ph)]
		min_max_rank_mouse = max_rank
		max_max_rank_mouse
		for mouse in mice_rep:
			data_mouse = data_ph[data_ph['Mouse']==mouse]
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
				ax_r.step(range(1, max_rank_mouse+1), x/largest, color = my_colors2[i_ph], alpha = .5, lw = 0.5)
			
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

		max_rank_eff = len(counts_per_ranking[counts_per_ranking>3])

		x_avg = x_avg[:max_rank_eff]/counts_per_ranking[:max_rank_eff]

		params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1)), np.log(x_avg))
		# print(np.sqrt(pcov))
		slope = params[0]
		zetas.append(-slope)
		zeta = 3*3.5/(4.5*2.1)
		# print(-slope)
					
	print(len(mice))
	for j in range(len(mice)):
		ax_r.lines[-(j+1)].set_color(my_colors_alpha[int(np.mean(zetas)*100)])

	ax_r.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int(np.mean(zetas)*100)], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = '^', label = r'$%.2f$'%(np.mean(zetas)))

	ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int(np.mean(zetas)*100)], alpha = .8, lw = 3)
	ax_zeta.hist(zetas, bins = np.linspace(0.2, 1.6, 20), alpha = .7, label = r"$\mathrm{" + ph + "}$", color = my_colors_alpha[int(np.mean(zetas)*100)], density = True, histtype = 'stepfilled', edgecolor = 'k')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_3.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
# ax_zeta.set_ylim(bottom = 2e-2, top = 1.1)
ax_zeta.set_xlim(left = 0.2, right = 1.6)
ax_zeta.legend(title = r'$\mathrm{sub-pop}$', fontsize = 30, title_fontsize = 30, loc = (1, 0))
fig_zeta.savefig(output_plot + '/zetas_3.pdf', transparent=.5)

#------------ Experiment 4 (Figure 5) ------------

fig_r_Flu, ax_r_Flu = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
fig_zeta_Flu, ax_zeta_Flu = plt.subplots(figsize=(10*1.62,8), gridspec_kw={'left':0.12, 'right':.8, 'bottom':.15, 'top': 0.94})

data_infection = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Influenza', header = 2)
data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort2', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort2', 'CDR3:']).size().reset_index(name='count')

mice = data_recall_grouped['Experiment / Mouse'].unique()
phenotypes = data_recall_grouped['Sort2'].unique()
print(phenotypes)

max_ranks = [100, 100, 100, 100, 100]
# max_ranks = [10, 10, 10, 10, 10]
for i_ph, ph in enumerate(phenotypes):
	max_rank = max_ranks[i_ph]
	zetas = []

	for rep in tqdm(range(n_ensemble)):

		if rep == n_ensemble - 1:
			mice_rep = mice
			# print(mice_rep)
		else:
			mice_rep = np.random.choice(mice, len(mice), replace = True)

		x_avg = np.zeros(max_rank)
		counts_per_ranking = np.zeros(max_rank)
		data_ph = data_recall_grouped[(data_recall_grouped['Sort2']==ph)]
		min_max_rank_mouse = max_rank
		max_max_rank_mouse = 0
		for mouse in mice_rep:
			data_mouse = data_ph[data_ph['Experiment / Mouse']==mouse]
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
				ax_r.step(range(1, max_rank_mouse+1), x/largest, color = my_colors3[i_ph], alpha = .5, lw = 0.5)
				ax_r_Flu.step(range(1, max_rank_mouse+1), x/largest, color = my_colors3[i_ph], alpha = .5, lw = 0.5)
			
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

		params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1)), np.log(x_avg))
		slope = params[0]

		zetas.append(-slope)

		zeta = 3*3.5/(4.5*2.1)
		
		if rep == n_ensemble - 1:
			ax_r.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int(np.mean(zetas)*100)], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = 'D', label = r'$%.2f$'%(np.mean(zetas)))
			ax_r_Flu.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int(np.mean(zetas)*100)], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = 'D', label = r'$%.2f$'%(np.mean(zetas)))

	print(len(mice))
	for j in range(len(mice)):
		ax_r.lines[-(j+1)].set_color(my_colors_alpha[int(np.mean(zetas)*100)])
		ax_r_Flu.lines[-(j+1)].set_color(my_colors_alpha[int(np.mean(zetas)*100)])

	ax_r.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int(np.mean(zetas)*100)], alpha = .8, lw = 3)
	ax_r_Flu.plot(np.arange(1, max_rank_eff + 1), np.exp(0)*np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int(np.mean(zetas)*100)], alpha = .8, lw = 3)
	
	ax_zeta.hist(zetas, bins = np.linspace(0.2, 1.6, 20), alpha = .7, label = r"$\mathrm{" + ph + "}$", color = my_colors_alpha[int(np.mean(zetas)*100)], density = True, histtype = 'stepfilled', edgecolor = 'k')
	ax_zeta_Flu.hist(zetas, bins = np.linspace(0.2, 1.6, 20), alpha = .7, label = r"$\mathrm{" + ph + "}$", color = my_colors_alpha[int(np.mean(zetas)*100)], density = True, histtype = 'stepfilled', edgecolor = 'k')

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 2e-2, top = 1.1)
ax_r.set_xlim(right = 5e1)
ax_r.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells_4.pdf', transparent=.5)

my_plot_layout(ax =ax_r_Flu, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r_Flu.set_ylim(bottom = 2e-2, top = 1.1)
ax_r_Flu.set_xlim(right = 5e1)
ax_r_Flu.legend(title = r'$\zeta$', fontsize = 30, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r_Flu.savefig(output_plot + '/ranking_B_cells_Flu.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
# ax_zeta.set_ylim(bottom = 2e-2, top = 1.1)
ax_zeta.set_xlim(left = 0.2, right = 1.6)
ax_zeta.legend(title = r'$\mathrm{sub-pop}$', fontsize = 30, title_fontsize = 30, loc = (1, 0))
fig_zeta.savefig(output_plot + '/zetas_4.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta_Flu, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
# ax_zeta_Flu.set_ylim(bottom = 2e-2, top = 1.1)
ax_zeta_Flu.set_xlim(left = 0.2, right = 1.6)
ax_zeta_Flu.legend(title = r'$\mathrm{sub-pop}$', fontsize = 30, title_fontsize = 30, loc = (1, 0))
fig_zeta_Flu.savefig(output_plot + '/zetas_Flu.pdf', transparent=.5)


# #------------ Experiment 4 (day 9) ------------

# data_infection = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Influenza9', header = 2)
# # data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort', 'V', 'J', 'D']).size().reset_index(name='count')
# data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort', 'CDR3:']).size().reset_index(name='count')
# # print(data_recall_grouped)
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
# 		# print(counts)
# 		N = np.sum(counts)
# 		S_i = -np.sum((counts/N)*np.log((counts/N)))

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
# 	print(np.sqrt(pcov))
# 	slope = params[0]

# 	zeta = 3*3.5/(4.5*2.1)
# 	print(-slope)
# 	ax_r.plot(np.arange(1, 30), np.exp(0)*np.arange(1, 30)**(slope), color = my_colors4[i_ph], alpha = 1)
# 	ax_r.plot(range(1, max_rank+1), x_avg, color = my_colors4[i_ph], alpha = .8, ls = '', marker = '.', label = r'$%.2f\pm %.3f$'%(-slope, np.sqrt(pcov[0][0])) + ' ; ' + ph)

#------------ Experiment 2 and 3 (Figure 4A and 4C) ------------

data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
data_recall = data_recall[(data_recall['Phenotype']=='GC + m')]

figures = ['4A', '4C-H']


max_rank = 100
zetas = []
for rep in tqdm(range(n_ensemble)):
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	# data_ph = data_recall_grouped[(data_recall_grouped['Phenotype']==ph)]
	min_max_rank_mouse = max_rank
	max_max_rank_mouse = 0

	len_mice = 0

	for i_fig, fig in enumerate(figures):
		data_recall_fig = data_recall[(data_recall['Figure']==fig)]
		data_recall_grouped = data_recall_fig.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
		# data_recall_grouped = data_recall.groupby(['Mouse', 'CDR3:']).size().reset_index(name='count')
		mice = data_recall_grouped['Mouse'].unique()

		if rep == n_ensemble - 1:
			mice_rep = mice
			# print(mice_rep)
		else:
			mice_rep = np.random.choice(mice, len(mice), replace = True)

		len_mice+=len(mice)
		for mouse in mice_rep:
			data_mouse = data_recall_grouped[data_recall_grouped['Mouse']==mouse]
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
				ax_r2.step(range(1, max_rank_mouse+1), x/largest, color = my_colors2[0], alpha = .5, lw = 0.5)
			
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

	params, pcov = curve_fit(model, np.log(range(1, max_rank_eff+1)), np.log(x_avg))
	# print(np.sqrt(pcov))
	slope = params[0]
	zetas.append(-slope)
	zeta = 3*3.5/(4.5*2.1)
	

print(len_mice)
for j in range(len_mice):
	ax_r2.lines[-(j+1)].set_color(my_colors_alpha[int(np.mean(zetas)*100)])

ax_r2.plot(range(1, max_rank_eff+1), x_avg, color = my_colors_alpha[int(np.mean(zetas)*100)], markerfacecolor="None", ms = 12, alpha = 1, ls = '', marker = 'o', label = r'$%.2f$'%(np.mean(zetas)))

ax_r2.plot(np.arange(1, max_rank_eff + 1), np.arange(1, max_rank_eff + 1)**(-np.mean(zetas)), color = my_colors_alpha[int(np.mean(zetas)*100)], alpha = .8, lw = 3)

ax_zeta2.hist(zetas, bins = np.linspace(0.2, 1.6, 20), alpha = .7, label = r'$\mathrm{GC+m}$', color = my_colors_alpha[int(np.mean(zetas)*100)], density = True, histtype = 'stepfilled', edgecolor = 'k')

my_plot_layout(ax =ax_r2, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r2.set_ylim(bottom = 2e-2, top = 1.1)
ax_r2.set_xlim(right = 5e1)
# ax_r2.legend(title = r'$\zeta$', fontsize = 24, title_fontsize = 30, loc = 3)#, loc = (1, 0))
fig_r2.savefig(output_plot + '/ranking_B_cells_proposal.pdf', transparent=.5)

my_plot_layout(ax =ax_zeta2, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
# ax_zeta2.set_ylim(bottom = 2e-2, top = 1.1)
ax_zeta2.set_xlim(left = 0.2, right = 1.6)
ax_zeta2.legend(title = r'$\mathrm{sub-pop}$', fontsize = 22, title_fontsize = 30, loc = (1, 0))
fig_zeta2.savefig(output_plot + '/zetas_proposal.pdf', transparent=.5)


my_plot_layout(ax =ax_r, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 1e-2, top = 1.1)
ax_r.legend(title = r'$\zeta$', fontsize = 11, title_fontsize = 13)
fig_r.savefig(output_plot + '/ranking_B_cells_linear.pdf', transparent=.2)

