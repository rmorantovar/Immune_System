import sys
sys.path.append('../../lib/')
from functions_memory import*
# from classes import*
#from functions_2 import*

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 1)
	parser.add_argument('--N_epi', type=int, default = 1)
	args = parser.parse_args()

	# Parameters -----------------------------------------------------

	N_ens = args.N_ens
	N_inf = args.N_inf
	N_evo = args.N_evo
	N_epi = args.N_epi

	L0 = int(1e6)  # Number of random sequences
	E_lim = -11.  # Threshold for the sum of entries
	t_lim = 8.  # Threshold for the sum of entries
	chunk_size = 1e6  # Size of each chunk
	p = 3
	k_step = 720
	n_jobs = -1

	lambda_A = 6.0
	lambda_B = 3 * np.log(2) #(days)^-1
	dT = 0.05
	C = 1e4
	time_array = np.linspace(0, 15, int((10-0)/dT))
	colors_inf = plt.cm.jet(np.linspace(0,1,N_inf))
	colors_epi = [my_blue, my_red]

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	#energy_model = 'MJ2'
	antigen = args.antigen
	epitopes = antigen.split('-')
	l=len(epitopes[0])

	data_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/memory/multi-epitope/output_N_ens_{N_ens}_L0_{L0}_p_{p}_k_step_{k_step}_E_lim_{E_lim}_t_lim_{t_lim}_lamA_{lambda_A}_n_evo_{N_evo}/"+ antigen
	
	fig, ax = plt.subplots(figsize=(4*5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

	fig_times, ax_times = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	
	ts_next = []
	ts_next = [5, 2, 10, 7, 32, 10, 2, 5, 9, 2, 7, 11, 1, 0, 20, 20, 29, 14, 20, 6, 8, 8, 15, 3, 16, 3, 35, 22, 32, 10, 9, 15, 16, 2, 5, 4, 23, 17, 0, 0, 4, 26, 34, 12, 20, 21, 16, 11, 10, 18, 23, 9, 17, 23, 16, 4, 4, 12, 2, 5, 0, 0, 33, 10, 1, 13, 14, 8, 3, 8, 4, 11, 36, 7, 13, 25, 40, 14, 1, 3, 16, 21, 10, 10, 13, 9, 31, 8, 1, 1]
	counter_t_next = 0
	for infection in tqdm(range(1, N_inf + 1)):
		data = pd.read_csv(data_dir + '/act_pop_'+ str(infection) + '.csv')
		if infection>1:
			mut_data = json.load(open(data_dir + '/mutation_'+ str(infection) + '.json'))
			for j in range(N_evo):
				epi_mut = mut_data['epi_mut_%d'%j]
				ax.scatter((infection-1)*time_array[-1], 2e4*(j+1), marker = 'D', color = colors_epi[epi_mut-1])
		# Find an input directory 

		Z = np.zeros_like(time_array)
		Z_naive = np.zeros_like(time_array)
		Z_memory = np.zeros_like(time_array)
		counter_0 = 0
		counter_1 = 0
		for i in range(N_ens):
			data_i = data.loc[data['ens_id']==i]
			data_i.reset_index(drop=True, inplace=True)			
			data_i = ensemble_of_expansions_time(data_i, N_ens, p, time_array, lambda_B, C, dT)
			
			Z_i = 0
			for epi in range(N_epi):
				data_epi = data_i.loc[data_i['epi']==epi+1]
				data_epi.reset_index(drop=True, inplace=True)	
				data_epi_0 = data_epi.loc[data_epi['m']==0]
				data_epi_1 = data_epi.loc[data_epi['m']==1]

				if(len(data_epi)>0):
					Z_epi = (data_epi['n_t']/np.exp(data_epi['E']))[data_epi.index].sum(axis = 0)
					# print((data_epi['n_t']/np.exp(data_epi['E']))[data_epi.index].max(axis = 0))
					# arg_max = np.argmax((data_epi['n_t']/np.exp(data_epi['E']))[data_epi.index][:, -1])
					# Z_epi_max = (data_epi['n_t'][-1]/np.exp(data_epi['E']))[data_epi.index][arg_max, :]
					Z_i += Z_epi
					
					# ax.plot(time_array + (infection-1)*time_array[-1], Z_epi_max, color = colors_epi[epi], alpha = .5, lw = 3)
					final_Z_epi_0 = 0
					final_Z_epi_1 = 0
					if len(data_epi_0.index)>0:
						counter_0+=1
						Z_epi_0 = (data_epi['n_t']/np.exp(data_epi['E']))[data_epi_0.index].sum(axis = 0)
						final_Z_epi_0 = Z_epi_0[-1]
						# Z_naive += Z_epi_0
						ax.plot(time_array + (infection-1)*time_array[-1], Z_epi_0, color = colors_epi[epi], alpha = .8, lw = 3, ls = '--')
					if len(data_epi_1.index)>0:
						counter_1+=1
						Z_epi_1 = (data_epi['n_t']/np.exp(data_epi['E']))[data_epi_1.index].sum(axis = 0)
						final_Z_epi_1 = Z_epi_1[-1]
						# Z_memory += Z_epi_1
						ax.plot(time_array + (infection-1)*time_array[-1], Z_epi_1, color = colors_epi[epi], alpha = .8, lw = 3, ls = ':')

					if final_Z_epi_1>final_Z_epi_0:
						counter_t_next+=1
					else:
						ts_next.append(counter_t_next)
						counter_t_next = 0

			ax.hlines(Z_i[-1], (time_array + (infection-1)*time_array[-1])[0], (time_array + (infection-1)*time_array[-1])[-1], color = 'grey', alpha = .5, lw = 3)
					

		# ax.plot(time_array + (infection-1)*time_array[-1], Z/N_ens, marker = '', color = 'k', alpha = .8, ls = '--')
		# ax.plot(time_array, Z_naive/counter_0 + Z_memory/counter_1, marker = '', color = colors_inf[infection], alpha = 1, label = '%d'%(infection+1), lw = 2, ls = '--')
		# ax.plot(time_array, Z_naive/counter_0, marker = '', color = colors_inf[infection], alpha = 1, ls = '--')
		# ax.plot(time_array, Z_memory/counter_1, marker = '', color = colors_inf[infection], alpha = 1, ls = ':')

		# ax_0.plot(time_array, Z_naive/counter_0, marker = '', color = colors_inf[infection], alpha = 1, ls = '--', label = '%d'%(infection+1))
		# ax_1.plot(time_array, Z_memory/counter_1, marker = '', color = colors_inf[infection], alpha = 1, ls = ':', label = '%d'%(infection+1))

	# ts_next = [1, 5, 2, 10, 7, 32, 10, 2, 5, 9, 2, 7, 11, 1, 0, 20, 20, 29, 14, 20, 6, 8, 8, 15, 3, 16, 3, 35, 22, 32, 10, 9, 15, 16, 2, 5, 4,
	# 23, 17, 0, 0, 4, 26, 34, 12, 20, 21, 16, 11, 10, 18, 23, 9, 17, 23, 16, 4, 4, 12, 2, 5, 0, 0,
	# 33, 10, 1, 13, 14, 8, 3, 8, 4, 11, 36, 7, 13, 25, 40, 14, 1, 3, 16, 21, 10, 10, 13, 9, 31, 8, 1, 1]
	ax_times.hist(ts_next[1:], bins = np.logspace(0, 2, 14))

	my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30, xlabel = 'infection' )
	# ax.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# ax.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
	# ax.set_ylim(bottom = 1e6, top = 1.1e11)
	ax.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	ax.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig.savefig('../../../Figures/memory/multi-epitope/Z_t_'+energy_model + '_' + antigen +'_'+str(N_evo)+'.pdf')

	my_plot_layout(ax = ax_times, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30, xlabel = 'infection' )
	# ax_times.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# ax_times.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
	# ax_times.set_ylim(bottom = 1e6, top = 1.1e11)
	# ax_times.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax_times.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig_times.savefig('../../../Figures/memory/multi-epitope/times_'+energy_model + '_' + antigen +'_'+str(N_evo)+'.pdf')

	# 	my_plot_layout(ax = ax_0, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	# 	ax_0.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# 	# ax_0.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
	# 	ax_0.set_ylim(bottom = 1e6, top = 6e11)
	# 	#ax_0.set_yticks([1, 0.1, 0.01, 0.001])
	# 	#ax_0.set_yticklabels([1, 0.1, 0.01])
	# 	fig_0.savefig('../../../Figures/memory/statistics/Z_stats'+energy_model+'_'+str(N_evo)+'_0.pdf')

	# 	my_plot_layout(ax = ax_1, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	# 	ax_1.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# 	# ax_1.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
	# 	ax_1.set_ylim(bottom = 1e6, top = 6e11)
	# 	#ax_1.set_yticks([1, 0.1, 0.01, 0.001])
	# 	#ax_1.set_yticklabels([1, 0.1, 0.01])
	# 	fig_1.savefig('../../../Figures/memory/statistics/Z_stats'+energy_model+'_'+str(N_evo)+'_1.pdf')


if __name__ == "__main__":
    main()

