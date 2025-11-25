import sys
sys.path.append('../../lib/')
from funcs import*
from matplotlib.colors import LogNorm
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.signal import medfilt

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=40, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = -1)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**8, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
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
	new = args.new

	if N_evo == -1:
		N_evo = 'R'


	E_lim = -11.  # Threshold for the sum of entries
	t_lim = 8.  # Threshold for the sum of entries
	chunk_size = 1e6  # Size of each chunk
	p = 3
	k_step = 1/(60*2) #s^-1
	k_step = k_step*3600 # hour^-1
	k_step = k_step*24 #days^-1
	k_on = 1e6*24*3600; #(M*days)^-1
	n_jobs = -1

	T = 12.
	lamA = 6.0
	lamB = 3 * np.log(2) #(days)^-1
	lamB = 2.
	dT = 0.002
	C = 1e5
	v = lamA/p
	b0 = 1e5
	t0 = np.log(lamA/(b0*k_on/N_A))/lamA
	Kstep = k_step/k_on
	ab = 1e-9
	time_array = np.linspace(0, T, int((T-0)/dT))#[::100]

	cmap = plt.get_cmap('cividis_r')
	# Create array of n evenly spaced values between 0 and 1
	color_vals = np.linspace(0, 1, 17)
	my_affinity_colors = [cmap(val) for val in color_vals] 
	color_vals2 = np.linspace(0, 1, 5)
	cmap2 = plt.get_cmap('plasma')
	my_pmem_colors = [cmap2(val) for val in color_vals2] 

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	#energy_model = 'MJ2'
	# antigen = args.antigen
	# epitopes = antigen.split('-')
	# l=len(epitopes[0])
	figEvo, axEvo = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
	figEvo2, axEvo2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})

	DDEs = np.linspace(0, 8, 17)
	# DDEs = np.linspace(0, 4.5, 10)
	pmems = [1., 1.5, 2, 2.5, 3.]
	# for i_exp, exp in enumerate([3, 4, 5, 6, 7]):
	for i_exp, exp in enumerate([0, 2, 3, 4, 5]):
		N_ens = args.N_ens
		pmem = pmems[i_exp]
		project = 'memory_response'
		subproject = 'multi-epitope'
		subproject = 'Z_dynamics'
		experiment = exp
		root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
		pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{int(k_step)}_lamA-{lamA}_lamB-{lamB}"
		pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}_N_evo-{N_evo}"
		antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv", converters={"antigen": literal_eval})
		antigens = antigens_data['antigen']
		output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
		os.makedirs(output_plot, exist_ok=True)

		N_ant = len(antigens[:])
		figDNA, axDNA = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
		figZ, axZ = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
		figDNAZ, axDNAZ = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
		figt, axt = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})

		DNAs_averages = []
		DZs_averages = []

		WT = antigens.iloc[0]
		dataNaive = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(1) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})
		#--------------------------Energy Motif--------------------------
		motif = get_motif(WT, energy_model, '../../')*1.2

		#Change values by the minimum
		E_m = -3
		for i in np.arange(l):
			E_m+=np.min(motif[:,i], axis=0)
			motif[:,i]-=np.min(motif[:,i], axis=0)
		# print('Em:', E_m)
		#--------------------------Entropy function--------------------------
		Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
		Kds = np.exp(Es[:-1])
		#--------------------------Repertoire properties--------------------------
		beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
		print(beta_r)
		tNaives = []
		N_A_peaks_Naive = []
		ZNaive = np.zeros_like(time_array)
		E_Naive = []
		# N_ens = 1
		for i in range(N_ens):
			#----------NAIVE----------
			dataNaive_i = dataNaive[dataNaive['ens_id']==i].reset_index()
			E_Naive.append(dataNaive_i['E'].min())
			tNaive = dataNaive_i['t'].min()
			dataNaive_i_t = ensemble_of_expansions_time(dataNaive_i, N_ens, p, time_array, lamB, C, dT)
			ZNaive_i = dataNaive_i_t.apply(lambda row: np.array(row['N_t']) / np.exp(row['E']), axis=1).sum()
			with np.errstate(divide='ignore', invalid='ignore'):
				pbNaive = (1+(1/(ab *ZNaive_i)))**(-1)
			# Interpolate f(t) for use in ODE
			pbNaive_interp = interp1d(time_array, pbNaive, kind='linear', fill_value='extrapolate')
			# Define the ODE: dN/dt = lambda * (1 - f(t)) * N
			def dNdtNaive(t, N):
			    return (lamA * (1 - pbNaive_interp(t)) - 1.) * N
			# Initial condition
			N0 = 1.0
			# Solve the ODE over the time span of your data
			solNaive = solve_ivp(dNdtNaive, t_span=(time_array[0], time_array[-1]), y0=[N0], t_eval=time_array)
			# Result
			N_ANaive_i = solNaive.y[0]  # solution N(t) evaluated at t_vals
			N_A_peak_Naive_i = np.max(N_ANaive_i)
			tNaives.append(tNaive)
			N_A_peaks_Naive.append(N_A_peak_Naive_i)
			with np.errstate(divide='ignore', invalid='ignore'):
				ZNaive+=np.log(ZNaive_i)

		ZNaive = np.exp(ZNaive/N_ens)
		axt.plot(time_array, ZNaive, color = my_blue, lw = 3, alpha = .8, ls = '-', label = r'$\textrm{Naive}$')
		axt.hlines(1/Kd_r, 0, T, color = 'grey', ls = '--')
		# axt.hlines(1/np.exp(np.mean(E_Naive)), 0, T, color = 'k', ls = '--')

		for i_DDE, DDE in enumerate(tqdm(DDEs)):
			ZMemory = np.zeros_like(time_array)
			DNAs = []
			DZs = []
					
			dataMemory = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(1) + '/%d'%(i_DDE) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})
			E_ms = []

			for i in range(N_ens):
				#----------MEMORY----------
				dataMemory_i = dataMemory[dataMemory['ens_id']==i].reset_index()
				tMemory = dataMemory_i['t'].min()
				dataMemory_i_t = ensemble_of_expansions_time(dataMemory_i, N_ens, p, time_array, lamB, C, dT)
				ZMemory_i = dataMemory_i_t.apply(lambda row: np.array(row['N_t']) / np.exp(row['E']), axis=1).sum()
				with np.errstate(divide='ignore', invalid='ignore'):
					pbMemory = (1+(1/(ab *ZMemory_i)))**(-1)
				# Interpolate f(t) for use in ODE
				pbMemory_interp = interp1d(time_array, pbMemory, kind='linear', fill_value='extrapolate')
				# Define the ODE: dN/dt = lambda * (1 - f(t)) * N
				def dNdtMemory(t, N):
				    return (lamA * (1 - pbMemory_interp(t)) - 1.) * N
				# Initial condition
				N0 = 1.0
				# Solve the ODE over the time span of your data
				solMemory = solve_ivp(dNdtMemory, t_span=(time_array[0], time_array[-1]), y0=[N0], t_eval=time_array)
				# Result
				N_AMemory_i = solMemory.y[0]  # solution N(t) evaluated at t_vals
				N_A_peak_Memory_i = np.max(N_AMemory_i)

				DNAs.append(np.log10(N_A_peak_Memory_i))#/N_A_peaks_Naive[i]))
				DZs.append(lamB*(tNaives[i] - tMemory))
				with np.errstate(divide='ignore', invalid='ignore'):
					ZMemory+=np.log(ZMemory_i)

			axDNA.hist(DNAs, color = my_affinity_colors[i_DDE], alpha = .8)
			axZ.hist(DZs, color = my_affinity_colors[i_DDE], alpha = .8)

			axDNAZ.scatter(DZs, DNAs, edgecolor = my_affinity_colors[i_DDE], facecolor="None", alpha = .2)
			axDNAZ.scatter(np.mean(DZs), np.mean(DNAs), edgecolor = my_affinity_colors[i_DDE], facecolor = my_affinity_colors[i_DDE], s = 30)

			ZMemory = np.exp(ZMemory/N_ens)
			
			if i_DDE%2 == 0:
				axt.plot(time_array, ZMemory, color = my_affinity_colors[i_DDE], lw = 2, alpha = 1, ls = '-')#, label = r'$\textrm{Memory}$')

			DNAs_averages.append(np.mean(DNAs))
			DZs_averages.append(np.mean(DZs) + np.log(100*(ZMemory[ZMemory>0][0])/(ZNaive[ZNaive>0][0])))

		my_plot_layout(ax = axt, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		axt.set_ylim(bottom = 1e6, top = 5e11)
		axt.set_xlim(left = 1, right = 8)
		axt.legend(fontsize = 18, loc = 2)
		figt.savefig(output_plot + '/Z_t_Evo.pdf')
				

		my_plot_layout(ax = axDNA, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# axDNA.set_xticks([])
		# axDNA.set_yticks([])
		# axDNA.set_ylim(bottom = 1e6, top = 2e12)
		# axDNA.set_xlim(left = 1, right = 8)
		# axDNA.legend(fontsize = 18, loc = 2)
		figDNA.savefig(output_plot + '/DNA.pdf')

		my_plot_layout(ax = axZ, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# axZ.set_xticks([])
		# axZ.set_yticks([])
		# axZ.set_ylim(bottom = 1e6, top = 2e12)
		# axZ.set_xlim(left = 1, right = 8)
		# axZ.legend(fontsize = 18, loc = 2)
		figZ.savefig(output_plot + '/DZ.pdf')

		my_plot_layout(ax = axDNAZ, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
		# axDNAZ.set_xticks([])
		# axDNAZ.set_yticks([])
		# axDNAZ.set_ylim(bottom = -7.5, top = -1)
		# axDNAZ.set_xlim(left = 3.5, right = 8.5)
		# axDNAZ.legend(fontsize = 18, loc = 2)
		figDNAZ.savefig(output_plot + '/DNAZ.pdf')

		axEvo.plot(DDEs, DZs_averages, color = my_pmem_colors[i_exp], ls = '-', lw = 2, marker = '', markeredgecolor = my_red, markerfacecolor=my_red, label = r'$%.1f$'%(pmem - beta_r))

		axEvo2.plot(DDEs, DNAs_averages, color = my_pmem_colors[i_exp], ls = '-', lw = 2, marker = '', markeredgecolor = my_red, markerfacecolor=my_red, label = r'$%.1f$'%(pmem - beta_r))


	subproject = 'Z_Evo'
	output_plot = '../../../Figures/'+project+'/'+subproject
	os.makedirs(output_plot, exist_ok=True)

	my_plot_layout(ax = axEvo, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# axEvo.set_xticks([])
	# axEvo.set_yticks([])
	# axEvo.set_ylim(bottom = -7.5, top = -1)
	# axEvo.set_xlim(left = 3.5, right = 8.5)
	axEvo.legend(fontsize = 14, loc = 0, title = r'$p_{\textrm{m}} - \beta^*$', title_fontsize = 16)
	figEvo.savefig(output_plot + '/DZ_DDE.pdf')

	my_plot_layout(ax = axEvo2, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# axEvo2.set_xticks([])
	# axEvo2.set_yticks([])
	# axEvo2.set_ylim(bottom = -7.5, top = -1)
	# axEvo2.set_xlim(left = 3.5, right = 8.5)
	axEvo2.legend(fontsize = 14, loc = 0, title = r'$p_{\textrm{m}} - \beta^*$', title_fontsize = 16)
	figEvo2.savefig(output_plot + '/DNA_DDE.pdf')

	
					

if __name__ == "__main__":
    main()

