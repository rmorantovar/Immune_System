import sys
sys.path.append('../../lib/')
from functions_1 import*
from classes import*
#from functions_2 import*

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=4, help="Number of infections.")
	parser.add_argument('--L0', type=int, default=10**8, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=20, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default= 5., help="Threshold for the sum of entries.")
	parser.add_argument('--E_lim', type=float, default= -8.0, help="Threshold for the sum of entries.")
	parser.add_argument('--E_m', type=float, default= -24, help="Threshold for the sum of entries.")
	parser.add_argument('--chunk_size', type=int, default=100000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=4, help="# steps.")
	parser.add_argument('--k_step', type=float, default=720, help="step rate.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
	parser.add_argument('--n_evo', type=int, default=0)
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

	n_evos = [0, 2, 4, 6]
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
				average = np.zeros_like(np.linspace(-20, -14+2, 50)[:-1])
				average_0 = np.zeros_like(np.linspace(-20, -14+2, 50)[:-1])
				average_1 = np.zeros_like(np.linspace(-20, -14+2, 50)[:-1])
				counter_0 = 0
				counter_1 = 0
				for i in range(N_ens):
					data_i = data.loc[data['ens_id']==i]
					hist= np.histogram(data_i['E'], bins = np.linspace(-20, -14+2, 50), density = False, weights = data_i['n'])
					average +=np.cumsum(hist[0])/np.sum(data_i['n'])
					ax.plot(np.exp(hist[1][:-1]), np.cumsum(hist[0])/np.sum(data_i['n']), marker = '', color = colors_inf[infection], alpha = .1, lw = .5)
					
					data_i_0 = data_i.loc[data_i['m']==0]
					if np.sum(data_i_0['n']) != 0:
						hist_0= np.histogram(data_i_0['E'], bins = np.linspace(-20, -14+2, 50), density = False, weights = data_i_0['n'])
						average_0 +=np.cumsum(hist_0[0])/np.sum(data_i_0['n'])
						ax_0.plot(np.exp(hist_0[1][:-1]), np.cumsum(hist_0[0])/np.sum(data_i_0['n']), marker = '', color = colors_inf[infection], alpha = .1, lw = .5)
						counter_0+=1


					data_i_1 = data_i.loc[data_i['m']==1]
					if np.sum(data_i_1['n']) != 0:
						hist_1= np.histogram(data_i_1['E'], bins = np.linspace(-20, -14+2, 50), density = False, weights = data_i_1['n'])
						average_1 +=np.cumsum(hist_1[0])/np.sum(data_i_1['n'])
						ax_1.plot(np.exp(hist_1[1][:-1]), np.cumsum(hist_1[0])/np.sum(data_i_1['n']), marker = '', color = colors_inf[infection], alpha = .1, lw = .5)
						counter_1+=1

				ax.plot(np.exp(hist[1][:-1]), average/N_ens, marker = '', color = colors_inf[infection], alpha = 1, label = '%d'%(infection+1))
				ax_0.plot(np.exp(hist_0[1][:-1]), average_0/counter_0, marker = '', color = colors_inf[infection], alpha = 1, label = '%d'%(infection+1))
				ax_1.plot(np.exp(hist_0[1][:-1]), average_1/counter_1, marker = '', color = colors_inf[infection], alpha = 1, label = '%d'%(infection+1))


		my_plot_layout(ax = ax, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		ax.legend(fontsize = 18, title_fontsize = 22, loc = 0, title = r'$\mathrm{Infection}$')
		ax.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
		#ax.set_ylim(bottom = 0.25)
		#ax.set_yticks([1, 0.1, 0.01, 0.001])
		#ax.set_yticklabels([1, 0.1, 0.01])
		fig.savefig('../../../Figures/memory/statistics/statistics'+energy_model+'_'+str(n_evo)+'.pdf')

		my_plot_layout(ax = ax_0, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		ax_0.legend(fontsize = 18, title_fontsize = 22, loc = 0, title = r'$\mathrm{Infection}$')
		ax_0.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
		#ax_0.set_ylim(bottom = 0.25)
		#ax_0.set_yticks([1, 0.1, 0.01, 0.001])
		#ax_0.set_yticklabels([1, 0.1, 0.01])
		fig_0.savefig('../../../Figures/memory/statistics/statistics'+energy_model+'_'+str(n_evo)+'_0.pdf')

		my_plot_layout(ax = ax_1, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		ax_1.legend(fontsize = 18, title_fontsize = 22, loc = 0, title = r'$\mathrm{Infection}$')
		ax_1.set_xlim(left = np.exp(-20.5), right = np.exp(-13+1.5))
		#ax_1.set_ylim(bottom = 0.25)
		#ax_1.set_yticks([1, 0.1, 0.01, 0.001])
		#ax_1.set_yticklabels([1, 0.1, 0.01])
		fig_1.savefig('../../../Figures/memory/statistics/statistics'+energy_model+'_'+str(n_evo)+'_1.pdf')


if __name__ == "__main__":
    main()

