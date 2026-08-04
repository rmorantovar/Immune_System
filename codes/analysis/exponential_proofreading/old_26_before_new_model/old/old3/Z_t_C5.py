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
	parser.add_argument('--N_ens', type=int, default=80, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = -1)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**7, help="Number of random sequences.")
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
	p = 4.0
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
	C = 2e4
	v = lamA/p
	b0 = 1e5
	t0 = np.log(lamA/(b0*k_on/N_A))/lamA
	Kstep = k_step/k_on
	ab = 1e-9
	time_array = np.linspace(0, T, int((T-0)/dT))#[::100]

	cmap = plt.get_cmap('cividis_r')
	# Create array of n evenly spaced values between 0 and 1
	color_vals = np.linspace(0, 1, 19)
	my_affinity_colors = [cmap(val) for val in color_vals] 
	color_vals2 = np.linspace(0, 1, 7)
	cmap2 = plt.get_cmap('managua')
	my_pmem_colors = [cmap2(val) for val in color_vals2] 

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	DDEs = [0]
	pmems = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
	exps = [1, 2, 3, 4, 5, 6, 7]
	# pmems = [4.0]
	# exps = [7]
	# figZ, axZ = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
	figP, axP = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
	figexp, axexp = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
	exponents = []

	for i_exp, exp in enumerate(tqdm(exps)):
		# print(i_exp)
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
		
		# subproject = 'Z_dynamics_C5'
		# output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
		# os.makedirs(output_plot, exist_ok=True)

		N_ant = len(antigens[:])
		# figt, axt = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
		# figt2, axt2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
		
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
		# print(beta_r)
		tNaives = []
		N_A_peaks_Naive = []
		ZNaive = np.zeros_like(time_array)
		larger_Naive = []
		Kd_r_sim = []
		max_rank = 5
		# N_ens = 5
		# print('%.1e ; %.1e'%(Kd_r, np.exp(np.mean(Kd_r_sim))))
		ZMemory = np.zeros_like(time_array)
		DNAs = []
		Zs = []
				
		dataMemory = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(1) + '/%d'%(0) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})
		E_ms = []
		j=0
		E_range = (-17, -13.5)

		U_eq = []
		SlogZeq = []

		U_p = []
		SlogZp = []

		U_q = []
		SlogZq = []
		# N_ens = 10
		for i in range(N_ens):
			# print(i)
			#----------NAIVE----------
			dataNaive_i = dataNaive[dataNaive['ens_id']==i].reset_index()
				
			# dataNaive_i_t = ensemble_of_expansions_time(dataNaive_i, N_ens, p, time_array, lamB, C, dT)
			# ZNaive_i = dataNaive_i_t.apply(lambda row: np.array(row['N_t']) / np.exp(row['E']), axis=1).sum()
			ZNaive_i_final = np.sum(dataNaive_i['N']/np.exp(dataNaive_i['E']))
			Zs.append(ZNaive_i_final)
			# with np.errstate(divide='ignore', invalid='ignore'):
			# 	pbNaive = (1+(1/(ab *ZNaive_i)))**(-1)
			# # Interpolate f(t) for use in ODE
			# pbNaive_interp = interp1d(time_array, pbNaive, kind='linear', fill_value='extrapolate')
			# # Define the ODE: dN/dt = lambda * (1 - f(t)) * N
			# def dNdtNaive(t, N):
			#     return (lamA * (1 - pbNaive_interp(t)) - 1.) * N
			# # Initial condition
			# N0 = 1.0
			# # Solve the ODE over the time span of your data
			# solNaive = solve_ivp(dNdtNaive, t_span=(time_array[0], time_array[-1]), y0=[N0], t_eval=time_array)
			# # Result
			# N_ANaive_i = solNaive.y[0]  # solution N(t) evaluated at t_vals
			# N_A_peak_Naive_i = np.max(N_ANaive_i)
			# tNaives.append(tNaive)
			# N_A_peaks_Naive.append(N_A_peak_Naive_i)
			# with np.errstate(divide='ignore', invalid='ignore'):
			# 	ZNaive+=np.log(ZNaive_i)
			# if i%4 ==0:
			# 	axt.plot(time_array, ZNaive_i, color = my_blue, lw = 1, alpha = .2, ls = '-')
			# 	axt2.plot(time_array, N_ANaive_i, color = my_blue, lw = 1.5, alpha = 1., ls = '-')

			#----------MEMORY----------
			dataMemory_i = dataMemory[dataMemory['ens_id']==i].reset_index()
			dataMemory_i_t = ensemble_of_expansions_time(dataMemory_i, N_ens, p, time_array, lamB, C, dT)
			dataMemory_i_t = dataMemory_i_t[dataMemory_i_t['m']==1].reset_index()
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

			DNAs.append((N_A_peak_Memory_i))#/N_A_peaks_Naive[i]))
			with np.errstate(divide='ignore', invalid='ignore'):
				ZMemory+=np.log(ZMemory_i)
		
		lamZ = np.max([lamB, lamA/pmem*(beta_r-1)-lamB])
		exponent, intercept = np.polyfit(np.log(np.array(Zs)), np.log(np.array(DNAs)), deg=1)
		exponents.append(exponent)
		print(lamA/((lamA/pmem)*(beta_r*0.9-1)-lamB), -exponent)
		axP.scatter(np.array(Zs)/np.min(Zs), np.array(DNAs)/np.max(DNAs), color = my_pmem_colors[i_exp], label = r'$%.1f$'%pmem)
		axP.plot(np.logspace(0, 1.3), 0.5*(np.logspace(0, 1.3))**(exponent), color = my_pmem_colors[i_exp], ls= '--')#, label = r'$N_A^c\sim\bar{\cal Z}^{\lambda_A/\lambda_{\hat{\cal Z}}}$')

	axexp.plot(pmems, -np.array(exponents), marker = 'D', color = my_blue, ls = '')
	x = np.where(lamA/np.linspace(1, 4, 50)*(beta_r-1)-lamB > lamB, lamA/(lamA/np.linspace(1, 4, 50)*(beta_r-1)-lamB), lamA/lamB)
	axexp.plot(np.linspace(1, 4, 50), x, color = 'k', ls = '--', label = r'$\hat\lambda_{\cal Z} = \max{[\lambda_B, \frac{\lambda_A}{p_{\textrm{mem}}}(\beta^*-1) - \lambda_B]}$')


	output_plot = '../../../Figures/'+project+'/'+subproject
	os.makedirs(output_plot, exist_ok=True)

	my_plot_layout(ax = axP, xscale='log', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# axRP.set_xticks([])
	# axRP.set_yticks([])
	# axRP.set_ylim(bottom = 4e-2)
	# axRP.set_xlim(left = -0.3, right = 7.3)
	axP.legend(fontsize = 20, loc = 0, title = r'$p_\textrm{mem}$', title_fontsize = 22)
	figP.savefig(output_plot + '/NAcZbar.pdf')

	my_plot_layout(ax = axexp, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# axexp.set_xticks([])
	# axexp.set_yticks([])
	# axexp.set_ylim(bottom = 4e-2)
	# axexp.set_xlim(left = -0.3, right = 7.3)
	axexp.legend(fontsize = 20, loc = 0)
	figexp.savefig(output_plot + '/exponents.pdf')

	
if __name__ == "__main__":
    main()

