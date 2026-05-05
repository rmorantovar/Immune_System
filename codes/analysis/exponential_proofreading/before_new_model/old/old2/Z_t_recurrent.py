import sys
sys.path.append('../../lib/')
from functions_memory import*
# from classes import*
#from functions_2 import*
from matplotlib.colors import LogNorm
from scipy.interpolate import splrep, BSpline

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=40, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=2, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 0)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--lamB', type=float, default=2., help="Antigen growth rate.")
	parser.add_argument('--new', type=int, default=0, help="run Z values again.")
	parser.add_argument('--exp', type=int, default=1, help="experiment.")
	args = parser.parse_args()

	# Parameters -----------------------------------------------------

	N_ant = args.N_ant
	N_ens = args.N_ens
	N_inf = args.N_inf
	N_evo = args.N_evo
	N_epi = args.N_epi
	L0 = args.L0
	l = args.l
	lamB = args.lamB
	new = args.new

	if N_evo == -1:
		N_evo = 'R'


	E_lim = -11.  # Threshold for the sum of entries
	t_lim = 8.  # Threshold for the sum of entries
	chunk_size = 1e6  # Size of each chunk
	p = 3
	k_step = 720
	n_jobs = -1

	lamA = 6.0
	T = 12
	dT = 0.01
	time_array = np.linspace(0, T, int((T-0)/dT))
	C = 1e4
	# time_array = np.linspace(0, 15, 1500)
	colors_inf = plt.cm.jet(np.linspace(0,1,N_inf))
	colors_mut = [my_blue, my_red]

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	#energy_model = 'MJ2'
	# antigen = args.antigen
	# epitopes = antigen.split('-')
	# l=len(epitopes[0])
	
	project = 'memory_response'
	subproject = 'multi-epitope'
	subproject = 'Recurrent'
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
	antigens = antigens_data['antigen']
	input_dir = root_dir + pars_dir_1 + pars_dir_2

	N_ant = len(antigens)
	for a1, antigen1 in enumerate(antigens):
		antigen_seq = from_aa_to_i(antigen1, energy_model, '../../')
		#--------------------------Energy Motif--------------------------
		motif = get_motif(antigen_seq, energy_model, '../../')*1.2

		#Change values by the minimum
		E_m = -3
		for i in np.arange(l):
			E_m+=np.min(motif[:,i], axis=0)
			motif[:,i]-=np.min(motif[:,i], axis=0)
		# if E_m<-22:
	
		# print('Em:', E_m)
		#--------------------------Entropy function--------------------------
		Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
		Kds = np.exp(Es[:-1])
		#--------------------------Repertoire properties--------------------------
		beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
		print("b* = ", beta_r)
		print("sigma = ", beta_r - 1 - lamB*p/lamA)
		print("sigma' = ", beta_r - 1 - 2*lamB*p/lamA)
		if new:
			sera = defaultdict(list)
			for inf in tqdm(range(N_inf)):
				input_dir = input_dir + "/%d"%(a1+1)
				if new:
					data1 = pd.read_csv(input_dir + '/activated_repertoire.csv')
					# for a11, antigen11 in enumerate(antigens):
					sera = run_essay_recurrent_time(data1, inf, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, beta_r, time_array, sera)	
			sera_df = pd.DataFrame(sera)
			sera_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Z_time.csv', index = False)
		

		sera_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Z_time.csv')
		
		# sera_df_grouped = group_by_clones(sera_df)
		fig1, ax1 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		fig2, ax2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
		fig3, ax3 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

		t_12s = []
		Z_12s = []
		exp_Zs = []
		v = lamA/p
		sigma = beta_r - 1 - 1*lamB/v
		for inf in range(N_inf):
			deltaT = 5
			Z = np.array(literal_eval(sera_df['Z'].iloc[inf]))
			lamZ = np.array(literal_eval(sera_df['lamZ'].iloc[inf]))
			Lact = np.array(literal_eval(sera_df['Lact'].iloc[inf]))
			if inf%3==0:
				ax1.plot(time_array, np.exp(np.sum(np.log(Z), axis = 0)/40), lw = 3, alpha = 1,  label = r'$%d$'%(inf+1))
				ax2.plot(time_array[::80][:-1], np.sum(lamZ, axis = 0)/40, lw = 3, alpha = 1, label = r'$%d$'%(inf+1))
				ax3.plot(time_array, np.exp(np.sum(np.log(Lact), axis = 0)/40), lw = 3, alpha = 1, label = r'$%d$'%(inf+1))
				for i in range(0, 40, 5):
					ax1.plot(time_array, Z[i, :], alpha = 0.3, color = ax1.get_lines()[-1].get_color())
					ax2.plot(time_array[::80][:-1], lamZ[i, :], alpha = 0.3, color = ax2.get_lines()[-1].get_color())
					ax3.plot(time_array, Lact[i, :], alpha = 0.3, color = ax2.get_lines()[-1].get_color())
					ax2.vlines(time_array[Z[i, :]>Z[i, :][-1]*0.5][0], 0, 4, ls = ':', color = ax2.get_lines()[-1].get_color(), alpha = 0.4)


		ax1.hlines(1/Kd_r, 0, 12, lw = 3, ls = '--', color = 'k')
		ax1.hlines(1e4/Kd_r, 0, 12, lw = 3, ls = '--', color = 'k')
		ax1.plot(time_array[:1000], (1/Kd_r)*np.exp(lamB*(time_array[:1000]-6)), lw = 3, ls = ':', color = 'k', label = r'$\exp[\lambda_B]$')
		ax1.plot(time_array[:800], (1/Kd_r)*np.exp(v*(beta_r - 1)*(time_array[:800]-4)), lw = 3, ls = ':', color = my_red, label = r'$\exp[v(\beta^* - 1)]$')

		ax2.hlines(v*(beta_r - 1), 0, T, lw = 3, ls = '--', color = my_red, label = r'$v(\beta^* - 1)$')
		ax2.hlines(lamB, 0, T, ls = '--', color = 'k', label = r'$\lambda_B$')

		ax3.plot(time_array[:900], 1*np.exp(v*beta_r*0.9*(time_array[:900]-6)), lw = 3, ls = ':', color = 'k')
		ax3.plot(time_array[:900], 1*np.exp(v*(beta_r - lamB/v)*0.9*(time_array[:900]-3)), lw = 3, ls = ':', color = my_green)

		my_plot_layout(ax = ax1, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax1.legend(fontsize = 14, title_fontsize = 16, loc = 0, title = r'$\mathrm{Infection}$')
		ax1.set_xlim(left = 2, right = 11)
		ax1.set_ylim(bottom = 2e5, top = 5e10)
		# ax1.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax1.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig1.savefig('../../../Figures/memory_response/Z_dynamics/Z1_t_rec_'+energy_model+'_L0-'+str(int(L0/10**int(np.log10(L0)))) + "e" + str(int(np.log10(L0))) + '_exp-'+str(experiment)+'.pdf')

		my_plot_layout(ax = ax2, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax2.legend(fontsize = 14, title_fontsize = 16, loc = 0, title = r'$\mathrm{Infection}$')
		ax2.set_xlim(left = 2, right = 11)
		ax2.set_ylim(bottom = -0.5, top = 10)
		# ax2.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax2.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig2.savefig('../../../Figures/memory_response/Z_dynamics/lamZ1_t_rec_'+energy_model+'_L0-'+str(int(L0/10**int(np.log10(L0)))) + "e" + str(int(np.log10(L0))) + '_exp-'+str(experiment)+'.pdf')

		my_plot_layout(ax = ax3, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		ax3.legend(fontsize = 14, title_fontsize = 20, loc = 0, )
		ax3.set_xlim(left = 2, right = 11)
		ax3.set_ylim(bottom = 0.5, top = 8e3)
		# ax3.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
		# ax3.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
		fig3.savefig('../../../Figures/memory_response/Z_dynamics/Lact_t_rec_'+energy_model+'_L0-'+str(int(L0/10**int(np.log10(L0)))) + "e" + str(int(np.log10(L0))) + '_exp-'+str(experiment)+'.pdf')


if __name__ == "__main__":
    main()

