import sys
sys.path.append('../../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/immune_system/primary_immune_response/'


data = pd.read_excel(Text_files_path + 'data/Victora_2016/aad3439-databases1.xlsx', header=1, sheet_name = 'Photoactivation')
data_early = data.loc[data['Figures']=='1/S2'][['Antigen', 'LN', 'GC', 'CDR3:']]

antigens = ['CGG']#, 'HA', 'NP-OVA'] #, 'OVA')]
colors_antigen = [my_blue, my_red, my_green, my_purple]
N_ensemble = 1000
max_rank = 24

antigen_slopes = []
antigen_err = []
fig_R, ax_R = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
for i_a, antigen in enumerate(tqdm(antigens)):
	data_antigen = data_early.loc[data['Antigen']==antigen]
	lns = np.unique(data_antigen['LN'])
	n_ln = len(lns)
	fig_S, ax_S = plt.subplots(1, n_ln, figsize = (5*n_ln, 6), gridspec_kw={'left':0.02, 'right':.95, 'bottom':.05, 'top': 0.95}, edgecolor = 'black')
	fig_r, ax_r = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	S = []
	L_act = []
	x_avg = np.zeros(max_rank)
	counts_per_ranking = np.zeros(max_rank)
	for i, ln in enumerate(lns):
		#print(str(ln))
		data_ln = data_antigen.loc[data_antigen['LN']==ln]
		CDR3, counts = np.unique(np.array((list(data_ln['CDR3:']))), return_counts = True)
		N = np.sum(counts)
		S_i = -np.sum((counts/N)*np.log((counts/N)))
		S.append(S_i)
		L_act.append(np.size(counts))
		if S_i != 0:
			ax_S[i].pie(counts)
			ax_S[i].set_title(r'$S_c=%.1f$'%(S_i), fontsize = 30)
		else:
			ax_S[i].axis('off')
			ax_S[i].text(x = .1, y = .6, s = r'$\bar S_c = %.1f \pm %.1f$'%(np.mean(S), np.sqrt(np.var(S))), fontsize = 30)


		if(N>0):
			sort_index = counts.argsort()
			largest = np.max(counts)
			x = np.flip(counts[sort_index])[:max_rank]
			for k in range(max_rank):
				if(x[k]>0):
					counts_per_ranking[k]+=1
					x_avg[k]+=x[k]/np.max(x)
			ax_r.plot(range(1,max_rank+1), x/np.max(x) , color = 'tab:red', alpha = .4)

	x_avg/=counts_per_ranking
	ax_r.plot(range(1,max_rank+1), x_avg, lw = 2, marker = '*', ls = '', color = 'darkred', ms = 12)
	ax_R.plot(range(1,max_rank+1), x_avg , alpha = 1, color = colors_antigen[i_a])

	fig_S.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2016/clonal_entropy_'+antigen+'.pdf')

	my_plot_layout(ax =ax_r, yscale = 'log', xscale = 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	# ax_r.set_ylim(bottom = 4e-2, top = 1.1)
	fig_r.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2016/ranking_'+antigen+'.pdf', transparent=.2)

	slopes_p = []

	for n in range(N_ensemble):
		counts_per_ranking = np.zeros(max_rank)
		x_p_avg = np.zeros(max_rank)
		#x_g_avg = np.zeros(max_rank)
		for i, ln in enumerate(lns):
			#print(str(ln))
			data_ln = data_antigen.loc[data_antigen['LN']==ln]
			CDR3, counts = np.unique(np.array((list(data_ln['CDR3:']))), return_counts = True)
			N = np.sum(counts)
			#print(counts)
			if(N>0):
				largest = np.max(counts)
				counts_p = np.random.poisson(lam = counts)
				sort_index = counts_p.argsort()
				#cs =  np.flip(counts[sort_index])[:max_rank]
				#cs_p = np.random.poisson(lam = cs)
				cs_p = np.flip(counts_p[sort_index])[:max_rank]
				#x_g = x + np.random.normal(0, np.sqrt(x))
				for k in range(max_rank):
					#if(cs_p[k]>0):
					#counts_per_ranking[k]+=1
					x_p_avg[k]+=np.max([cs_p[k], 1])/cs_p[0]
					# x_p_avg[k]+=np.max([cs_p[k], 1])/np.max(cs_p)
				#x_g_avg[k]+=x_g[k]/np.max(x_g)
			    #Y.append(np.log(x) + dy)
			    #X.append(np.log(max_rank+1))
		x_p_avg/=n_ln
		popt_p, pcov_p = curve_fit(my_linear_func, np.log(range(1,max_rank+1)), np.log(x_p_avg))
		#popt_g, pcov = curve_fit(my_linear_func, np.log(range(1,max_rank+1)), np.log(x_g_avg))
		slopes_p.append(popt_p[1])
		#slopes_g.append(popt_g[1])

	ax_R.plot(np.linspace(1,max_rank+1, 100), np.linspace(1,max_rank+1, 100)**(np.mean(slopes_p)), lw = 5, marker = '',
		ls = '-', color = colors_antigen[i_a], alpha = .6, label = antigen + r' $(%.2f\pm%.2f)$'%(np.mean(slopes_p), np.std(slopes_p)))
	antigen_slopes.append(np.mean(slopes_p))
	antigen_err.append(np.std(slopes_p))

print(antigen_slopes, np.mean(antigen_slopes))
my_plot_layout(ax =ax_R, yscale = 'log', xscale = 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
ax_R.set_ylim(bottom = 4e-2, top = 1.1)
ax_R.legend(title = r'$\mathrm{antigen}$', fontsize = 20, title_fontsize = 24)
fig_R.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2016/ranking.pdf', transparent=.2)


antigen_slopes = np.array(antigen_slopes)
lameff = -1/antigen_slopes
w = np.array([160])#, 65, 50])#, 45])
x = w**(-1/3)*np.log(1/w)**2
w_array = np.linspace(np.min(w), np.max(w), 100)
x_array = np.linspace(np.min(x), np.max(x), 100)
popt, pcov = curve_fit(my_linear_func, np.log(w), np.log(lameff))
popt2, pcov2 = curve_fit(my_linear_func, np.log(x), np.log(lameff))
print(popt)
print(popt2)
print(popt2[1], np.exp(popt[0]))

fig_a, ax_a = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
ax_a.scatter(w, (antigen_slopes), color = my_red, lw = 5)
ax_a.errorbar(w, (antigen_slopes), yerr=antigen_err, color = my_red, ls = '')
#ax_a.plot(w_array, np.exp(popt[1]*w_array), label = r'$w^{$.1f}$'%(popt[1]))
my_plot_layout(ax =ax_a, yscale = 'linear', xscale = 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30)
#ax_a.set_ylim(bottom = 4e-2, top = 1.1)
fig_a.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2016/zeta_vs_w.pdf', transparent=.2)

fig_a, ax_a = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
ax_a.scatter(w, lameff, color = my_red, lw = 5)
ax_a.errorbar(w, lameff, yerr=antigen_err/(-antigen_slopes)**2, color = my_red, ls = '')
ax_a.plot(w_array, np.exp(popt[0])*w_array**popt[1], label = r'$w^{%.1f}$'%(popt[1]), color = my_red, lw = 5, alpha = .4)
my_plot_layout(ax =ax_a, yscale = 'linear', xscale = 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30)
#ax_a.set_ylim(bottom = 4e-2, top = 1.1)
ax_a.legend(fontsize = 20)
fig_a.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2016/lamAeff_vs_w.pdf', transparent=.2)

fig_a, ax_a = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
ax_a.scatter(x, lameff, color = my_red, lw = 5)
ax_a.errorbar(x, lameff, yerr=antigen_err/(-antigen_slopes)**2, color = my_red, ls = '')
ax_a.plot(x_array, np.exp(popt2[0])*x_array**popt2[1], label = r'$x^{%.1f}$'%(popt2[1]), color = my_red, lw = 5, alpha = .4)
my_plot_layout(ax =ax_a, yscale = 'linear', xscale = 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30)
#ax_a.set_ylim(bottom = 4e-2, top = 1.1)
ax_a.legend(fontsize = 20)
fig_a.savefig('../../../../Figures/primary_immune_response/11_data_Victora/2016/lamAeff_vs_x.pdf', transparent=.2)













    