import sys
sys.path.append('../../my_lib/')
from funcs import*
from matplotlib.colors import LogNorm
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=80, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = -1)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--p', type=float, default=3.0, help="# steps.")
	parser.add_argument('--pmem', type=float, default=1.0, help="# steps for memory.")
	parser.add_argument('--DDE', type=float, default=0.0, help="DDE value.")
	parser.add_argument('--new', type=int, default=0, help="run Z values again.")
	parser.add_argument('--exp', type=int, default=1, help="experiment.")
	parser.add_argument('--one_WT', type=int, default = 1)
	args = parser.parse_args()

	# Parameters -----------------------------------------------------

	N_ant = args.N_ant
	N_ens = args.N_ens
	N_inf = args.N_inf
	N_evo = args.N_evo
	N_epi = args.N_epi
	L0 = args.L0
	l = args.l
	p = args.p
	pmem = args.pmem
	DDE = args.DDE
	new = args.new

	if N_evo == -1:
		N_evo = 'R'


	E_lim = -11.  # Threshold for the sum of entries
	t_lim = 8.  # Threshold for the sum of entries
	chunk_size = 1e6  # Size of each chunk
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
	time_array = np.linspace(0, T, int((T-0)/dT))#[::100]
	colors_inf = plt.cm.jet(np.linspace(0,1,N_inf))
	colors_mut = [my_blue, my_red]

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	project = 'memory_response'
	subproject = 'multi-epitope'
	subproject = 'Z_NA_dynamics'
	one_WT = args.one_WT


	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{int(k_step)}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv", converters={"antigen": literal_eval})
	if one_WT:
		WTs = antigens_data.iloc[[0]]
	else:
		WTs = antigens_data
	output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
	os.makedirs(output_plot, exist_ok=True)

	N_ant = len(WTs[:])

	for index, row in WTs.iterrows():
		antigen1 = row['antigen']
		a1 = index

		figZ, axZ = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
		figN_A, axN_A = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.95, 'bottom':.15, 'top': 0.94})
		# antigen_seq = from_aa_to_i(antigen1, energy_model, '../../')
		#--------------------------Energy Motif--------------------------
		motif = get_motif(antigen1, energy_model, '../../')*1.2

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

		dataNaive = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})
		dataMemory = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/DDE_%.1f_pmem_%.1f'%(DDE, pmem) + '/activated_repertoire.csv')#, converters={"N_t": literal_eval})
		E_ms = []

		ZNaive = np.zeros_like(time_array)
		NNaive = np.zeros_like(time_array)
		ZMemory = np.zeros_like(time_array)
		NMemory = np.zeros_like(time_array)
		N_ANaive = np.zeros_like(time_array)
		N_AMemory = np.zeros_like(time_array)

		# N_ens = 2
		for i in tqdm(range(N_ens)):
			#----------NAIVE----------
			dataNaive_i = dataNaive[dataNaive['ens_id']==i].reset_index()
			E_ms.append(np.min(dataNaive_i['E']))
			dataNaive_i_t = ensemble_of_expansions_time(dataNaive_i, N_ens, p, time_array, lamB, C, dT)
			N_t_array = np.vstack(dataNaive_i_t['N_t'].values)  # shape: (num_rows, list_length)
			E_array = np.array(dataNaive_i_t['E'].values)       # shape: (num_rows,)

			# ZNaive_i = dataNaive_i_t.apply(lambda row: np.array(row['N_t']) / np.exp(row['E']), axis=1).sum()
			# NNaive_i = dataNaive_i_t.apply(lambda row: np.array(row['N_t']), axis=1).sum()

			ZNaive_i = (N_t_array / np.exp(E_array[:, None])).sum(axis=0)
			NNaive_i = N_t_array.sum(axis=0)

			with np.errstate(divide='ignore', invalid='ignore'):
				pbNaive = (1+(1/(1e-9 *ZNaive_i)))**(-1)

			# Interpolate f(t) for use in ODE
			pbNaive_interp = interp1d(time_array, pbNaive, kind='linear', fill_value='extrapolate')

			# Define the ODE: dN/dt = lambda * (1 - f(t)) * N
			def dNdtNaive(t, N):
				return (lamA * (1 - pbNaive_interp(t)) - 1.5) * N
			
			# Initial condition
			N0 = 1.0
			
			# Solve the ODE over the time span of your data
			solNaive = solve_ivp(dNdtNaive, t_span=(time_array[0], time_array[-1]), y0=[N0], t_eval=time_array)

			# Result
			N_ANaive_i = solNaive.y[0]  # solution N(t) evaluated at t_vals

			if i%20==0:
				axZ.plot(ZNaive_i, N_ANaive_i, color = my_blue, lw = .2, alpha = .2, ls = '', marker = 'o', ms = 2)
				axN_A.plot(time_array, N_ANaive_i, color = my_blue, lw = 2, alpha = .5, ls = '-')

			#----------MEMORY----------
			dataMemory_i = dataMemory[dataMemory['ens_id']==i].reset_index()
			# E_ms.append(np.min(dataMemory_i['E']))
			dataMemory_i_t = ensemble_of_expansions_time(dataMemory_i, N_ens, p, time_array, lamB, C, dT)
			N_t_array = np.vstack(dataMemory_i_t['N_t'].values)  # shape: (num_rows, list_length)
			E_array = np.array(dataMemory_i_t['E'].values)       # shape: (num_rows,)

			ZMemory_i = (N_t_array / np.exp(E_array[:, None])).sum(axis=0)
			NMemory_i = N_t_array.sum(axis=0)

			with np.errstate(divide='ignore', invalid='ignore'):
				pbMemory = (1+(1/(1e-9 *ZMemory_i)))**(-1)

			# Interpolate f(t) for use in ODE
			pbMemory_interp = interp1d(time_array, pbMemory, kind='linear', fill_value='extrapolate')

			# Define the ODE: dN/dt = lambda * (1 - f(t)) * N
			def dNdtMemory(t, N):
				return (lamA * (1 - pbMemory_interp(t)) - 1.5) * N

			# Initial condition
			N0 = 1.0

			# Solve the ODE over the time span of your data
			solMemory = solve_ivp(dNdtMemory, t_span=(time_array[0], time_array[-1]), y0=[N0], t_eval=time_array)

			# Result
			N_AMemory_i = solMemory.y[0]  # solution N(t) evaluated at t_vals

			if i%20==0:
				axZ.plot(ZMemory_i, N_AMemory_i, color = my_purple, lw = .2, alpha = .2, ls = '', marker = 'o', ms = 2)
				axN_A.plot(time_array, N_AMemory_i, color = my_purple, lw = 2, alpha = .5, ls = '-')
							
			with np.errstate(divide='ignore', invalid='ignore'):
				NNaive+=np.log(NNaive_i)
				NMemory+=np.log(NMemory_i)
				ZNaive+=np.log(ZNaive_i)
				ZMemory+=np.log(ZMemory_i)
				N_ANaive+=np.log(N_ANaive_i)
				N_AMemory+=np.log(N_AMemory_i)

		NNaive = np.exp(NNaive/N_ens)
		NMemory = np.exp(NMemory/N_ens)
		ZNaive = np.exp(ZNaive/N_ens)
		ZMemory = np.exp(ZMemory/N_ens)
		N_ANaive = np.exp(N_ANaive/N_ens)
		N_AMemory = np.exp(N_AMemory/N_ens)

		axN_A.plot(time_array, N_ANaive, color = my_blue, lw = 4, alpha = .8, ls = '-', label = r'$\textrm{Naive}$')
		axN_A.plot(time_array, N_AMemory, color = my_purple, lw = 4, alpha = .8, ls = '-', label = r'$\textrm{Memory}$')

		axZ.plot(ZNaive, N_ANaive, color = my_blue, lw = 4, alpha = .8, ls = '-', marker = '', ms = 5, label = r'$\textrm{Naive}$')
		axZ.plot(ZMemory, N_AMemory, color = my_purple, lw = 4, alpha = .8, ls = '-', marker = '', ms = 5, label = r'$\textrm{Memory}$')

		lamZ = np.max([lamB, lamA/p*(beta_r-1)])
		logging.info(f'lamZ: {lamZ}')

		Z_mf = Kstep**(1/v)/(Kd_r)**(1+1/v)*np.exp(lamZ*(time_array - t0))
		# axZ.plot(time_array[(time_array<8) & (time_array>t0)], Z_mf[(time_array<8) & (time_array>t0)], alpha = .5, color = my_blue, lw = 2, ls = '--')#, label = r'$\textrm{Mean-field}\, \textrm{approx}$')
		Z_mf = Kstep**(1/v)/(Kd_r)**(1+1/v)*np.exp(lamB*(time_array - t0))
		# axZ.plot(time_array[(time_array<8) & (time_array>t0)], Z_mf[(time_array<8) & (time_array>t0)], alpha = .5, color = my_blue, lw = 2, ls = ':')#, label = r'$\textrm{Mean-field}\, \textrm{approx}$')

		# axZ.plot(np.logspace(0, 4.1, 100), (1/(Kd_r))*np.logspace(0, 4.1, 100), alpha = 1, color = 'k', lw = 2, ls = ':')#, label = r'$\textrm{Mean-field}\, \textrm{approx}$')
		
		if pmem == 4.0:
			Z_mf = 0.1*C*Kstep**(1/v)/(Kd_r)**(1+1/v)*np.exp(lamB*(time_array - t0))
			# axZ.plot(time_array[(time_array<6) & (time_array>t0)], Z_mf[(time_array<6) & (time_array>t0)], alpha = .5, color = my_purple, lw = 2, ls = '--')#, label = r'$\textrm{Mean-field}\, \textrm{approx}$')
		
		# axZ.hlines(1/Kd_r, 0, T, color = 'k', ls = '--')

		my_plot_layout(ax = axZ, xscale='log', yscale= 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30)
		# axZ.set_xticks([])
		# axZ.set_yticks([])
		axZ.set_ylim(bottom = 1e4, top = 2e13)
		# axZ.set_xlim(left = 1, right = 8)
		axZ.legend(fontsize = 30, loc = 2)
		figZ.savefig(output_plot + '/Z_t_p-%.1f_pmem_%.1f_DDE_%.1f.pdf'%(p, pmem, DDE))

		my_plot_layout(ax = axN_A, xscale='linear', yscale= 'log', ticks_labelsize= 40, x_fontsize=30, y_fontsize=30)
		# axN_A.set_xticks([])
		# axN_A.set_yticks([])
		axN_A.set_ylim(bottom = 1e4, top = 2e13)
		axN_A.set_xlim(left = 1, right = 10)
		axN_A.legend(fontsize = 30, loc = 2)
		figN_A.savefig(output_plot + '/N_A_t_p-%.1f_pmem_%.1f_DDE_%.1f.pdf'%(p, pmem, DDE))


if __name__ == "__main__":
	main()

