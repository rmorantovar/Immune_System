import sys
sys.path.append('../../lib/')
from functions_1 import*
from classes import*
#from functions_2 import*

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=2, help="Number of infections.")
	parser.add_argument('--L0', type=int, default=10**8, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=20, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default=8., help="Threshold for activation time.")
	parser.add_argument('--E_lim', type=float, default=-11., help="Threshold for the sum of entries.")
	parser.add_argument('--E_m', type=float, default=-24, help="Threshold for the sum of entries.")
	parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=3, help="# steps.")
	parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
	parser.add_argument('--lamA', type=float, default=6., help="Antigen growth rate.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
	parser.add_argument('--random_antigen', type=int, default=0)
	# parser.add_argument('--antigen', type=str, default='TACNSEYPNTTRAKCGRWYR')
	parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR')
	parser.add_argument('--energy_model', type=str, default = 'TCRen')
	parser.add_argument('--n_evo', type=int, default = 0)
	parser.add_argument('--seqs', type=int, default = 1)
	args = parser.parse_args()


	# Parameters -----------------------------------------------------

	N_ens = args.N_ens  # Number of times to execute the process
	N_inf = args.N_inf
	L0 = args.L0  # Number of random sequences
	E_lim = args.E_lim  # Threshold for the sum of entries
	t_lim = args.t_lim  # Threshold for the sum of entries
	E_m = args.E_m  #
	chunk_size = args.chunk_size  # Size of each chunk
	p = args.p
	k_step = args.k_step
	n_jobs = args.n_jobs
	n_evo = args.n_evo

	lambda_A = 6.0
	lambda_B = 3 * np.log(2) #(days)^-1
	dT = 0.05
	C = 1e4
	time_array = np.linspace(0, 10, int((10-0)/dT))
	colors_inf = plt.cm.jet(np.linspace(0,1,N_inf))
	colors_inf = [my_green, my_blue, my_red]

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	#energy_model = 'MJ2'
	antigens = ['TACNSEYPNTTRAKCGRWYC', 'TACNSEYPNTTFDKCGRWYC', 'MACNSEYPNTTRAKCGRWYC', 'MACNSEYPNTTRCKCLRWYC', 'YACNSEYPNTTFDKCGRWYC', 'TACNSTYPNTERAKCGRWYC', 'MACNSEYPNTTRCRKLRWYC',
	'TACNSKYPNDTTDKCGRWSC', 'MACNSECPNTTRCKWLRWYC', 'MACNSEYPNTTRCKEFRWYC'] #L=20
	antigens = ['TACNSEYPNTTRAKCGRWYR']

	n_evos = [0, 2, 4, 6, 8]
	for n_evo in n_evos:
		# Find an input directory 
		data_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/memory/exploration/output_N_ens_{N_ens}_L0_{L0}_p_{p}_k_step_{k_step}_E_lim_{E_lim}_t_lim_{t_lim}_lamA_{lambda_A}_n_evo_{n_evo}"

		fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		fig_0, ax_0= plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		fig_1, ax_1= plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		for antigen in antigens:    
			antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
			l=len(antigen)

			for infection in range(N_inf):
				data = pd.read_csv(data_dir + '/act_pop_'+ str(infection+1) + '.csv')
				Z = np.zeros_like(time_array)
				Z_naive = np.zeros_like(time_array)
				Z_memory = np.zeros_like(time_array)
				counter_0 = 0
				counter_1 = 0
				for i in range(N_ens):
					data_i = data.loc[data['ens_id']==i]
					data_i.reset_index(drop=True, inplace=True)
					data_i_0 = data_i.loc[data_i['m']==0]
					data_i_1 = data_i.loc[data_i['m']==1]
					
					# print(data_i.index, data_i_0.index, data_i_1.index)

					data_i = ensemble_of_expansions_time(data_i, N_ens, p, time_array, lambda_B, C, dT)
					
					Z_i = (data_i['n_t']/np.exp(data_i['E']))[data_i.index].sum(axis = 0)
					
					Z += Z_i
					
					ax.plot(time_array, Z_i, color = colors_inf[infection], alpha = .1, lw = 1)

					if len(data_i_0.index)>0:
						counter_0+=1
						Z_i_0 = (data_i['n_t']/np.exp(data_i['E']))[data_i_0.index].sum(axis = 0)
						Z_naive += Z_i_0
						ax_0.plot(time_array, Z_i_0, color = colors_inf[infection], alpha = .1, lw = 1)
					if len(data_i_1.index)>0:
						counter_1+=1
						Z_i_1 = (data_i['n_t']/np.exp(data_i['E']))[data_i_1.index].sum(axis = 0)
						Z_memory += Z_i_1
						ax_1.plot(time_array, Z_i_1, color = colors_inf[infection], alpha = .1, lw = 1)
					

				ax.plot(time_array, Z/N_ens, marker = '', color = colors_inf[infection], alpha = .8, label = '%d'%(infection+1), lw = 2)
				# ax.plot(time_array, Z_naive/counter_0 + Z_memory/counter_1, marker = '', color = colors_inf[infection], alpha = 1, label = '%d'%(infection+1), lw = 2, ls = '--')
				# ax.plot(time_array, Z_naive/counter_0, marker = '', color = colors_inf[infection], alpha = 1, ls = '--')
				# ax.plot(time_array, Z_memory/counter_1, marker = '', color = colors_inf[infection], alpha = 1, ls = ':')

				ax_0.plot(time_array, Z_naive/counter_0, marker = '', color = colors_inf[infection], alpha = 1, ls = '--', label = '%d'%(infection+1))
				ax_1.plot(time_array, Z_memory/counter_1, marker = '', color = colors_inf[infection], alpha = 1, ls = ':', label = '%d'%(infection+1))


		my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		ax.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
		# ax.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
		ax.set_ylim(bottom = 1e6, top = 6e11)
		#ax.set_yticks([1, 0.1, 0.01, 0.001])
		#ax.set_yticklabels([1, 0.1, 0.01])
		fig.savefig('../../../Figures/memory/statistics/Z_stats'+energy_model+'_'+str(n_evo)+'.pdf')

		my_plot_layout(ax = ax_0, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		ax_0.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
		# ax_0.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
		ax_0.set_ylim(bottom = 1e6, top = 6e11)
		#ax_0.set_yticks([1, 0.1, 0.01, 0.001])
		#ax_0.set_yticklabels([1, 0.1, 0.01])
		fig_0.savefig('../../../Figures/memory/statistics/Z_stats'+energy_model+'_'+str(n_evo)+'_0.pdf')

		my_plot_layout(ax = ax_1, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		ax_1.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
		# ax_1.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
		ax_1.set_ylim(bottom = 1e6, top = 6e11)
		#ax_1.set_yticks([1, 0.1, 0.01, 0.001])
		#ax_1.set_yticklabels([1, 0.1, 0.01])
		fig_1.savefig('../../../Figures/memory/statistics/Z_stats'+energy_model+'_'+str(n_evo)+'_1.pdf')


if __name__ == "__main__":
    main()

