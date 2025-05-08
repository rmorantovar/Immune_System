import sys
sys.path.append('../../../lib/')
from funcs import*
plt.rcParams['text.usetex'] = True

project = 'memory_response'
subproject = 'data'
experiment = 0
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/mesin2020"
output_plot = '../../../../Figures/'+project+'/'+subproject+'/'+str(experiment) + '/mesin2020'
os.makedirs(output_plot, exist_ok=True)

# Parameters
n_ens = 10000
gs = [2]  # Number of Poisson processes
mu = 1.0  # Poisson rate
T = 15  # Total simulation time
theta = 1.5  # Values of theta to compare
gamma = 0.4
my_colors = [my_red, my_blue2, my_blue2, my_blue2, my_cyan]
my_colors2 = [my_blue2, my_blue2, my_blue2, my_cyan, my_blue2, my_red]
my_colors3 = [my_red, my_blue2, my_blue2, my_red, my_red, my_blue2, my_blue2, my_red, my_red, my_red]
my_colors4 = [my_red, my_blue2, my_blue2, my_red, my_red, my_blue2, my_blue2, my_red, my_red, my_red]
alpha = 1e-10
depth = 6
anti_mut_epi = 5/4

def model(x, m):
    return m * x 

fig_r, ax_r = plt.subplots(figsize=(10*1.62,8), gridspec_kw={'left':0.12, 'right':.8, 'bottom':.15, 'top': 0.94})

#------------ Experiment 1 (Figure 1D) ------------

data_primary = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Photoactivation CGG', header = 1)
data_primary = data_primary[(data_primary['Figure']==1)]
# data_primary_grouped = data_primary.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
data_primary_grouped = data_primary.groupby(['Mouse', 'CDR3:']).size().reset_index(name='count')
mice = data_primary_grouped['Mouse'].unique()
# phenotypes = data_primary_grouped['Phenotype'].unique()

# for i_ph, ph in enumerate(phenotypes):
# 	max_rank = max_ranks[i_ph]
max_rank = 30
x_avg = np.zeros(max_rank)
counts_per_ranking = np.zeros(max_rank)
for mouse in mice:
	data_mouse = data_primary_grouped[data_primary_grouped['Mouse']==mouse]
	# CDR3, counts = np.unique(np.array((list(data_mouse['CDR3:']))), return_counts = True)
	counts = data_mouse['count'].to_numpy()
	# print(counts)
	N = np.sum(counts)
	S_i = -np.sum((counts/N)*np.log((counts/N)))

	if(N>0):
		sort_index = counts.argsort()
		largest = np.max(counts)
		x = np.flip(counts[sort_index])
		if len(x)>max_rank:
			x = x[:max_rank]
		else:
			x = np.pad(x, (0, max_rank - len(x)), mode='constant')
		ax_r.step(range(1, max_rank+1), x/largest, color = my_colors[0], alpha = .5, lw = 0.5)
		for k in range(max_rank):
			if(x[k]>0):
				counts_per_ranking[k]+=1
				x_avg[k]+=x[k]/largest

x_avg/=counts_per_ranking

# Linear fit
coeffs = np.polyfit(np.log(range(1, max_rank+1)), np.log(x_avg), 1)  # 1 = degree of the polynomial (linear)
slope, intercept = coeffs

params, pcov = curve_fit(model, np.log(range(1, max_rank+1)), np.log(x_avg))
print(np.sqrt(pcov))
slope = params[0]

zeta = 3*3.5/(4.5*2.1)
print(-slope)
ax_r.plot(np.arange(1, 30), np.exp(0)*np.arange(1, 30)**(slope), color = my_colors[0], alpha = 1)
ax_r.plot(range(1, max_rank+1), x_avg, color = my_colors[0], alpha = .8, ls = '', marker = '^', label = r'$%.2f$'%(-slope) + ' ; ' + 'GC')

#------------ Experiment 2 (Figure 4A) ------------

data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
data_recall = data_recall[(data_recall['Figure']=='4A')]
# data_recall_grouped = data_recall.groupby(['Mouse', 'V', 'J', 'D']).size().reset_index(name='count')
data_recall_grouped = data_recall.groupby(['Mouse', 'CDR3:']).size().reset_index(name='count')
mice = data_recall_grouped['Mouse'].unique()
# phenotypes = data_recall_grouped['Phenotype'].unique()

# for i_ph, ph in enumerate(phenotypes):
max_rank = 18
x_avg = np.zeros(max_rank)
counts_per_ranking = np.zeros(max_rank)
# data_ph = data_recall_grouped[(data_recall_grouped['Phenotype']==ph)]
for mouse in mice:
	data_mouse = data_recall_grouped[data_recall_grouped['Mouse']==mouse]
	# CDR3, counts = np.unique(np.array((list(data_mouse['CDR3:']))), return_counts = True)
	counts = data_mouse['count'].to_numpy()
	# print(counts)
	N = np.sum(counts)
	S_i = -np.sum((counts/N)*np.log((counts/N)))

	if(N>0):
		sort_index = counts.argsort()
		largest = np.max(counts)
		x = np.flip(counts[sort_index])
		if len(x)>max_rank:
			x = x[:max_rank]
		else:
			x = np.pad(x, (0, max_rank - len(x)), mode='constant')
		ax_r.step(range(1, max_rank+1), x/largest, color = my_colors2[0], alpha = .5, lw = 0.5)
		for k in range(max_rank):
			if(x[k]>0):
				counts_per_ranking[k]+=1
				x_avg[k]+=x[k]/largest

x_avg/=counts_per_ranking

# Linear fit
coeffs = np.polyfit(np.log(range(1, max_rank+1)), np.log(x_avg), 1)  # 1 = degree of the polynomial (linear)
slope, intercept = coeffs

params, pcov = curve_fit(model, np.log(range(1, max_rank+1)), np.log(x_avg))
print(np.sqrt(pcov))
slope = params[0]

zeta = 3*3.5/(4.5*2.1)
print(-slope)
ax_r.plot(np.arange(1, 30), np.exp(intercept)*np.arange(1, 30)**(slope), color = my_colors2[0], alpha = 1)
ax_r.plot(range(1, max_rank+1), x_avg, color = my_colors2[0], alpha = .8, ls = '', marker = '*', label = r'$%.2f$'%(-slope) + ' ; ' + 'GC + m')

#------------ Experiment 3 (Figure 4C) ------------

data_recall = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Fate-mapping CGG', header = 1)
data_recall = data_recall[(data_recall['Figure']=='4C-H')]
# data_recall_grouped = data_recall.groupby(['Mouse', 'Phenotype', 'V', 'J', 'D']).size().reset_index(name='count')
data_recall_grouped = data_recall.groupby(['Mouse', 'Phenotype', 'CDR3:']).size().reset_index(name='count')
# print(data_recall_grouped)
mice = data_recall_grouped['Mouse'].unique()
phenotypes = data_recall_grouped['Phenotype'].unique()

max_ranks = [20, 20]
for i_ph, ph in enumerate(phenotypes):
	max_rank = max_ranks[i_ph]
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	data_ph = data_recall_grouped[(data_recall_grouped['Phenotype']==ph)]
	for mouse in mice:
		data_mouse = data_ph[data_ph['Mouse']==mouse]
		# CDR3, counts = np.unique(np.array((list(data_mouse['CDR3:']))), return_counts = True)
		counts = data_mouse['count'].to_numpy()
		# print(counts)
		N = np.sum(counts)
		S_i = -np.sum((counts/N)*np.log((counts/N)))

		if(N>0):
			sort_index = counts.argsort()
			largest = np.max(counts)
			x = np.flip(counts[sort_index])
			if len(x)>max_rank:
				x = x[:max_rank]
			else:
				x = np.pad(x, (0, max_rank - len(x)), mode='constant')
			ax_r.step(range(1, max_rank+1), x/largest, color = my_colors2[i_ph], alpha = .5, lw = 0.5)
			for k in range(max_rank):
				if(x[k]>0):
					counts_per_ranking[k]+=1
					x_avg[k]+=x[k]/largest

	x_avg/=counts_per_ranking

	# Linear fit
	coeffs = np.polyfit(np.log(range(1, max_rank+1)), np.log(x_avg), 1)  # 1 = degree of the polynomial (linear)
	slope, intercept = coeffs

	params, pcov = curve_fit(model, np.log(range(1, max_rank+1))[:-5], np.log(x_avg)[:-5])
	print(np.sqrt(pcov))
	slope = params[0]

	zeta = 3*3.5/(4.5*2.1)
	print(-slope)
	ax_r.plot(np.arange(1, 30), np.exp(0)*np.arange(1, 30)**(slope), color = my_colors2[i_ph], alpha = 1)
	ax_r.plot(range(1, max_rank+1), x_avg, color = my_colors2[i_ph], alpha = .8, ls = '', marker = 'D', label = r'$%.2f$'%(-slope) + ' ; ' + ph)

#------------ Experiment 4 (Figure 5) ------------

data_infection = pd.read_excel(root_dir + "/1-s2.0-S0092867419313170-mmc1.xlsx", sheet_name = 'Influenza', header = 2)
# data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort2', 'V', 'J', 'D']).size().reset_index(name='count')
data_recall_grouped = data_infection.groupby(['Experiment / Mouse', 'Sort2', 'CDR3:']).size().reset_index(name='count')

mice = data_recall_grouped['Experiment / Mouse'].unique()
phenotypes = data_recall_grouped['Sort2'].unique()
print(phenotypes)

max_ranks = [30, 22, 25, 25, 25]
for i_ph, ph in enumerate(phenotypes):
	max_rank = max_ranks[i_ph]
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	data_ph = data_recall_grouped[(data_recall_grouped['Sort2']==ph)]
	for mouse in mice:
		data_mouse = data_ph[data_ph['Experiment / Mouse']==mouse]
		# CDR3, counts = np.unique(np.array((list(data_mouse['CDR3:']))), return_counts = True)
		counts = data_mouse['count'].to_numpy()
		# print(counts)
		N = np.sum(counts)
		S_i = -np.sum((counts/N)*np.log((counts/N)))

		if(N>0):
			sort_index = counts.argsort()
			largest = np.max(counts)
			x = np.flip(counts[sort_index])
			if len(x)>max_rank:
				x = x[:max_rank]
			else:
				x = np.pad(x, (0, max_rank - len(x)), mode='constant')
			ax_r.step(range(1, max_rank+1), x/largest, color = my_colors3[i_ph], alpha = .5, lw = 0.5)
			for k in range(max_rank):
				if(x[k]>0):
					counts_per_ranking[k]+=1
					x_avg[k]+=x[k]/largest

	x_avg/=counts_per_ranking

	# Linear fit
	coeffs = np.polyfit(np.log(range(1, max_rank+1)), np.log(x_avg), 1)  # 1 = degree of the polynomial (linear)
	slope, intercept = coeffs
	params, pcov = curve_fit(model, np.log(range(1, max_rank+1))[:-10], np.log(x_avg)[:-10])
	print(np.sqrt(pcov))
	slope = params[0]

	zeta = 3*3.5/(4.5*2.1)
	print(-slope)
	ax_r.plot(np.arange(1, 30), np.exp(0)*np.arange(1, 30)**(slope), color = my_colors3[i_ph], alpha = 1)
	ax_r.plot(range(1, max_rank+1), x_avg, color = my_colors3[i_ph], alpha = .8, ls = '', marker = 'o', label = r'$%.2f$'%(-slope) + ' ; ' + ph)

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

# 	# Linear fit
# 	coeffs = np.polyfit(np.log(range(1, max_rank+1)), np.log(x_avg), 1)  # 1 = degree of the polynomial (linear)
# 	slope, intercept = coeffs

# 	params, pcov = curve_fit(model, np.log(range(1, max_rank+1)), np.log(x_avg))
# 	print(np.sqrt(pcov))
# 	slope = params[0]

# 	zeta = 3*3.5/(4.5*2.1)
# 	print(-slope)
# 	ax_r.plot(np.arange(1, 30), np.exp(0)*np.arange(1, 30)**(slope), color = my_colors4[i_ph], alpha = 1)
# 	ax_r.plot(range(1, max_rank+1), x_avg, color = my_colors4[i_ph], alpha = .8, ls = '', marker = '.', label = r'$%.2f\pm %.3f$'%(-slope, np.sqrt(pcov[0][0])) + ' ; ' + ph)

my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 1e-2, top = 1.1)
ax_r.legend(title = r'$\zeta$', fontsize = 24, title_fontsize = 26, loc = (1, 0))
fig_r.savefig(output_plot + '/ranking_B_cells.pdf', transparent=.2)

my_plot_layout(ax =ax_r, yscale = 'linear', xscale = 'linear', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30 )
ax_r.set_ylim(bottom = 1e-2, top = 1.1)
ax_r.legend(title = r'$\zeta$', fontsize = 11, title_fontsize = 13)
fig_r.savefig(output_plot + '/ranking_B_cells_2.pdf', transparent=.2)

