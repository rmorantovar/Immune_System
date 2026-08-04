import sys
sys.path.append('../../lib/')
from functions_memory import*
# from classes import*
#from functions_2 import*
from scipy.signal import savgol_filter

# CODE TO CALCULATE PHENOTYPE EXPONENT FOR 1ST AND 2ND INFECTION

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 0)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--L0', type=int, default=10**6, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--new', type=int, default=0, help="run Z values again.")
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

	E_lim = -8.  # Threshold for the sum of entries
	t_lim = 8.  # Threshold for the sum of entries
	chunk_size = 1e6  # Size of each chunk
	p = 3
	k_step = 720
	n_jobs = -1

	lamA = 6.0
	lamB = 3 * np.log(2) #(days)^-1
	lamB = 2.
	dT = 0.05
	C = 1e4
	k_on = 1e6*24*3600; #(M*days)^-1
	b0 = 1e5*10
	time_array = np.linspace(0, 15, int((10-0)/dT))
	colors_inf = plt.cm.jet(np.linspace(0,1,N_inf))
	colors_mut = [my_blue, my_red]

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	#energy_model = 'MJ2'
	# antigen = args.antigen
	# epitopes = antigen.split('-')
	# l=len(epitopes[0])
	
	project = 'memory'
	subproject = 'multi-epitope'
	subproject = 'PS'
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}"
	pars_dir_1 = f"/L0-{L0}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
	antigens = antigens_data['antigen']

	N_ant = len(antigens)
	
	beta_star_l_th_1 = []
	beta_star_l_sim_1 = []
	beta_star_c_th_1 = []
	beta_star_c_sim_1 = []

	beta_star_l_th_2 = []
	beta_star_l_sim_2 = []
	beta_star_c_th_2 = []
	beta_star_c_sim_2 = []

	beta_star_l_th_2_memory = []
	beta_star_l_sim_2_memory = []
	beta_star_c_th_2_memory = []
	beta_star_c_sim_2_memory = []

	for a, antigen in enumerate(tqdm(antigens[::1])):
		antigen_seq = from_aa_to_i(antigen, energy_model, '../../')
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
		# print('beta_r = %.1f'%beta_r, 'K_r = %.1e'%Kd_r, 'E_r = %.1f'%E_r)
		#--------------------------Proofreading properties--------------------------
		# beta_step, E_step, Kd_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
		# print('beta_step = %.2f'%beta_step)

		D0_lineages_p = {}
		DR_lineages_p = {}

		D0_cells_p = {}
		DR_cells_p = {}

		D_theory = {}

		D_approx = {}
		# print(a, antigen)
		#-----------------Loading data----------------------------
		data1, return_data_type1 = get_data(folder_path = '../../out/memory/'+str(L0)+'/'+str(N_ens), data_type = 'cells_antigen-%d_1'%a)

		if(return_data_type1):
			energies_lineages1 = data1[0]
			energies_cells1 = data1[1]
			L_act1 = data1[2]
			final_times1 = data1[3]
		else:
			try:
				data1 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a+1) + '/activated_repertoire.csv')
				# print('File 1 read')
				energies_lineages1 =[]
				energies_cells1 =[]
				L_act1 = 0
				final_times1 = []
				counter = 0
				for i_ens in np.arange(N_ens):
					data_active = data1.loc[data1['ens_id']==i_ens]
					t_act_data = np.min(data_active['t'])
					# data_active = data_active.loc[data_active['t']<(t_act_data+1.2+0.2*(p-1))] # it was 1.0 + 0.1*...
					activation_times = np.array(data_active['t'])
					energies  = np.array(data_active['E'])
					#---------------------------- B cell linages ----------------------
					clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lamB, C, dT)
					#--------------------------t_C filter-------------------------
					lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
					clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
					#-------Simulations-------
					if(len(energies_C)>0):
						L_act_i = len(energies_C)
						counter+=1
						energies_lineages1+=list(energies_C)
						final_times1.append(np.max(activation_times_C))
						for j in range(len(energies_C)):
							energies_cells1+=([energies_C[j]]*int(clone_sizes_C[j, -1]))

					L_act1+=L_act_i

				os.makedirs('../../out/memory/'+str(L0)+'/'+str(N_ens), exist_ok=True)
				f = open('../../out/memory/'+str(L0)+'/'+str(N_ens)+'/processed_data_cells_antigen-%d_1.pkl'%a, 'wb')
				pickle.dump([energies_lineages1, energies_cells1, L_act1/counter, final_times1], f, pickle.HIGHEST_PROTOCOL)

			except FileNotFoundError:
				print(f'Data file does not exist')
				continue
		
		# Second infection
		data2, return_data_type2 = get_data(folder_path = '../../out/memory/'+str(L0)+'/'+str(N_ens), data_type = 'cells_antigen-%d_2'%a)

		if(return_data_type2):
			energies_lineages2 = data2[0]
			energies_cells2 = data2[1]
			L_act2 = data2[2]
			final_times2 = data2[3]
			energies_lineages2_memory = data2[4]
			energies_cells2_memory = data2[5]
			L_act2_memory = data2[6]
			final_times2_memory = data2[7]
		else:
			try:
				data2 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a+1) + '/%d'%(a+1) + '/activated_repertoire.csv')
				# print('File 2 read')
				energies_lineages2 =[]
				energies_cells2 =[]
				L_act2 = 0
				final_times2 = []
				energies_lineages2_memory =[]
				energies_cells2_memory =[]
				L_act2_memory = 0
				final_times2_memory = []
				counter = 0
				for i_ens in np.arange(N_ens):
					data_active = data2.loc[data2['ens_id']==i_ens]
					data_memory = data_active.loc[data_active['m']==1]
					
					t_act_data = np.min(data_active['t'])
					# data_active = data_active.loc[data_active['t']<(t_act_data+1.2+0.2*(p-1))] # it was 1.0 + 0.1*...
					activation_times = np.array(data_active['t'])
					energies  = np.array(data_active['E'])
					#---------------------------- B cell linages ----------------------
					clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lamB, C, dT)
					#--------------------------t_C filter-------------------------
					lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
					clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
					#-------Simulations-------
					if(len(energies_C)>0):
						L_act_i = len(energies_C)
						counter+=1
						energies_lineages2+=list(energies_C)
						final_times2.append(np.max(activation_times_C))
						for j in range(len(energies_C)):
							energies_cells2+=([energies_C[j]]*int(clone_sizes_C[j, -1]))

					L_act2+=L_act_i

					t_act_data_memory = np.min(data_memory['t'])
					# data_memory = data_memory.loc[data_memory['t']<(t_act_data_memory+1.2+0.2*(p-1))] # it was 1.0 + 0.1*...
					activation_times_memory = np.array(data_memory['t'])
					energies_memory  = np.array(data_memory['E'])
					#---------------------------- B cell linages ----------------------
					clone_sizes_memory = get_clones_sizes_C(len(activation_times_memory), time_array, activation_times_memory, lamB, C, dT)
					#--------------------------t_C filter-------------------------
					lim_size_memory = np.max([int(np.max(clone_sizes_memory[:, -1])*0.01), 2])
					clone_sizes_C_memory, activation_times_C_memory, energies_C_memory, filter_C_memory, n_C_memory = apply_filter_C(clone_sizes_memory, activation_times_memory, energies_memory, lim_size_memory)
					#-------Simulations-------
					if(len(energies_C_memory)>0):
						L_act_i_memory = len(energies_C_memory)
						energies_lineages2_memory+=list(energies_C_memory)
						final_times2_memory.append(np.max(activation_times_C_memory))
						for j in range(len(energies_C_memory)):
							energies_cells2_memory+=([energies_C_memory[j]]*int(clone_sizes_C_memory[j, -1]))

					L_act2_memory+=L_act_i_memory

				os.makedirs('../../out/memory/'+str(L0)+'/'+str(N_ens), exist_ok=True)
				f = open('../../out/memory/'+str(L0)+'/'+str(N_ens)+'/processed_data_cells_antigen-%d_2.pkl'%a, 'wb')
				pickle.dump([energies_lineages2, energies_cells2, L_act2/counter, final_times2, energies_lineages2_memory, energies_cells2_memory, L_act2_memory/counter, final_times2_memory], f, pickle.HIGHEST_PROTOCOL)
			except FileNotFoundError:
				# print(f'Unprocessed data does not exist')
				continue

		n_coarse = 1
		#-----------
		bins = np.linspace(np.min(energies_lineages1), -8, 40)
		data_lineages1 = np.histogram(energies_lineages1, density = False, bins = bins)#, cumulative = True, alpha = 0)
		bins = np.linspace(np.min(energies_cells1), -8, 40)
		data_cells1 = np.histogram(energies_cells1, density = False, bins = bins)#, cumulative = True, alpha = 0)
		
		Es_l1 = data_lineages1[1][:-1]
		cum_func = np.interp(Es_l1, Es_l1, np.cumsum(data_lineages1[0]))
		cum_func = np.cumsum(data_lineages1[0])
		# cum_func = savgol_filter(np.cumsum(data_lineages1[0]), int(len(Es_l1)/2), 2)
		P_l1 = np.diff(cum_func)/np.diff(Es_l1)#/(len(energies_lineages1))
		P_l1 /= np.sum(P_l1*np.diff(Es_l1))
		# plt.plot(Es_l1[:-1], P_l1)
		P_l1 = (np.diff(np.log(cum_func))/np.diff(Es_l1))*cum_func[:-1]#/(len(energies_lineages1))
		P_l1 /= np.sum(P_l1*np.diff(Es_l1))

		# P_l1 = np.interp(Es[:-2], data_lineages1[1][:-2], P_l1)

		Es_c1 = data_cells1[1][:-1]
		P_c1 = np.diff(np.cumsum(data_cells1[0]))/np.diff(Es_c1)#/(len(energies_cells1))
		P_c1 /= np.sum(P_c1*np.diff(Es_c1))

		size_l1 = len(Es_l1[:-1][(Es_l1[:-1]<Es_l1[:-1][np.argmax(P_l1)]) & (P_l1!=0)][:])
		size_c1 = len(Es_c1[:-1][(Es_c1[:-1]<Es_c1[:-1][np.argmax(P_c1)]) & (P_c1!=0)][:])

		# P_c1 = np.interp(Es[:-2], data_cells1[1][:-2], P_c1)
		Q_01 = Q0[::n_coarse][:-1]#*L0#/len(energies_lineages1)

		u_on, p_a, R, QR12 = calculate_QR(Q0, k_on, k_step, np.exp(lamA*(np.mean(final_times1)*1.0))/N_A, Es, p, lamA, b0, dE)
		Q_R12 = QR12[::n_coarse][:-1]#*L0#/len(energies_lineages1)
		Q_R12 = Q_R12/np.sum(Q_R12*dE[::n_coarse][:-1])

		#-----------
		bins = np.linspace(np.min(energies_lineages2), -11, 40)
		data_lineages2 = np.histogram(energies_lineages2, density = False, bins = bins)#, cumulative = True, alpha = 0)
		bins = np.linspace(np.min(energies_cells2), -11, 40)
		data_cells2 = np.histogram(energies_cells2, density = False, bins = bins)#, cumulative = True, alpha = 0)
		
		Es_l2 = data_lineages2[1][:-1]
		P_l2 = np.diff(np.cumsum(data_lineages2[0]))/np.diff(Es_l2)#/(len(energies_lineages2))
		P_l2 /= np.sum(P_l2*np.diff(Es_l2))
		#P_l2 = np.interp(Es[:-2], data_lineages2[1][:-2], P_l2)

		Es_c2 = data_cells2[1][:-1]
		P_c2 = np.diff(np.cumsum(data_cells2[0]))/np.diff(Es_c2)#/(len(energies_cells2))
		P_c2 /= np.sum(P_c2*np.diff(Es_c2))

		size_l2 = len(Es_l2[:-1][(Es_l2[:-1]<Es_l2[:-1][np.argmax(P_l2)]) & (P_l2!=0)][:])
		size_c2 = len(Es_c2[:-1][(Es_c2[:-1]<Es_c2[:-1][np.argmax(P_c2)]) & (P_c2!=0)][:])

		#-----------
		bins = np.linspace(np.min(energies_lineages2_memory), -11, 30)
		data_lineages2_memory = np.histogram(energies_lineages2_memory, density = False, bins = bins)#, cumulative = True, alpha = 0)
		bins = np.linspace(np.min(energies_cells2_memory), -11, 30)
		data_cells2_memory = np.histogram(energies_cells2_memory, density = False, bins = bins)#, cumulative = True, alpha = 0)
		
		Es_l2_memory = data_lineages2_memory[1][:-1]
		P_l2_memory = np.diff(np.cumsum(data_lineages2_memory[0]))/np.diff(Es_l2_memory)#/(len(energies_lineages2))
		P_l2_memory /= np.sum(P_l2_memory*np.diff(Es_l2_memory))
		#P_l2 = np.interp(Es[:-2], data_lineages2[1][:-2], P_l2)

		Es_c2_memory = data_cells2_memory[1][:-1]
		P_c2_memory = np.diff(np.cumsum(data_cells2_memory[0]))/np.diff(Es_c2_memory)#/(len(energies_cells2))
		P_c2_memory /= np.sum(P_c2_memory*np.diff(Es_c2_memory))

		size_l2_memory = len(Es_l2_memory[:-1][(Es_l2_memory[:-1]<Es_l2_memory[:-1][np.argmax(P_l2_memory)]) & (P_l2_memory!=0)][:])
		size_c2_memory = len(Es_c2_memory[:-1][(Es_c2_memory[:-1]<Es_c2_memory[:-1][np.argmax(P_c2_memory)]) & (P_c2_memory!=0)][:])

		try:

			popt_l_1, pcov_l_1 = curve_fit(my_linear_func, Es_l1[:-1][(Es_l1[:-1]<Es_l1[:-1][np.argmax(P_l1)]) & (P_l1!=0)][:int(4*size_l1/6)], np.log(P_l1[(Es_l1[:-1]<Es_l1[:-1][np.argmax(P_l1)]) & (P_l1!=0)][:int(4*size_l1/6)]))
			popt_c_1, pcov_c_1 = curve_fit(my_linear_func, Es_c1[:-1][(Es_c1[:-1]<Es_c1[:-1][np.argmax(P_c1)]) & (P_c1!=0)][:int(4*size_c1/6)], np.log(P_c1[(Es_c1[:-1]<Es_c1[:-1][np.argmax(P_c1)]) & (P_c1!=0)][:int(4*size_c1/6)]))

			popt_l_2, pcov_l_2 = curve_fit(my_linear_func, Es_l2[:-1][(Es_l2[:-1]<Es_l2[:-1][np.argmax(P_l2)]) & (P_l2!=0)][:int(4*size_l2/6)], np.log(P_l2[(Es_l2[:-1]<Es_l2[:-1][np.argmax(P_l2)]) & (P_l2!=0)][:int(4*size_l2/6)]))
			popt_c_2, pcov_c_2 = curve_fit(my_linear_func, Es_c2[:-1][(Es_c2[:-1]<Es_c2[:-1][np.argmax(P_c2)]) & (P_c2!=0)][:int(4*size_c2/6)], np.log(P_c2[(Es_c2[:-1]<Es_c2[:-1][np.argmax(P_c2)]) & (P_c2!=0)][:int(4*size_c2/6)]))

			popt_l_2_memory, pcov_l_2_memory = curve_fit(my_linear_func, Es_l2_memory[:-1][(Es_l2_memory[:-1]<Es_l2_memory[:-1][np.argmax(P_l2_memory)]) & (P_l2_memory!=0)][:int(4*size_l2_memory/6)], np.log(P_l2_memory[(Es_l2_memory[:-1]<Es_l2_memory[:-1][np.argmax(P_l2_memory)]) & (P_l2_memory!=0)][:int(4*size_l2_memory/6)]))
			popt_c_2_memory, pcov_c_2_memory = curve_fit(my_linear_func, Es_c2_memory[:-1][(Es_c2_memory[:-1]<Es_c2_memory[:-1][np.argmax(P_c2_memory)]) & (P_c2_memory!=0)][:int(4*size_c2_memory/6)], np.log(P_c2_memory[(Es_c2_memory[:-1]<Es_c2_memory[:-1][np.argmax(P_c2_memory)]) & (P_c2_memory!=0)][:int(4*size_c2_memory/6)]))

			beta_star_l_th_1.append(beta_r)
			beta_star_l_sim_1.append(popt_l_1[1])
			beta_star_c_th_1.append(beta_r - (lamB*p)/(lamA*1))
			beta_star_c_sim_1.append(popt_c_1[1])

			beta_star_l_th_2.append(beta_r)
			beta_star_l_sim_2.append(popt_l_2[1])
			beta_star_c_th_2.append(beta_r - (lamB*p)/(lamA*1))
			beta_star_c_sim_2.append(popt_c_2[1])

			beta_star_l_th_2_memory.append(beta_r)
			beta_star_l_sim_2_memory.append(popt_l_2_memory[1])
			beta_star_c_th_2_memory.append(beta_r - (lamB*p)/(lamA*1))
			beta_star_c_sim_2_memory.append(popt_c_2_memory[1])

		except:
			print('Fail to fit exponent in population %d!'%a)

		if a%40==0:
			# print(beta_r, popt_l_1[1])
			fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
			ax.plot(Es[::n_coarse][:-2][::1], Q_01[::1], color = 'grey', ls = '--')#, label = r'$\Omega_{0}$', lw = 2)
			# ax.plot(Es[::n_coarse][:-2][::1], Q_R1[::1], color = 'black', ls = '--')#, label = r'$\Omega_{\textrm{act}}$')
			ax.plot(Es[::n_coarse][:-2][::1], Q_R12[::1], color = 'black', ls = '--')#, label = r'$\Omega_{\textrm{act}}$')

			ax.plot(Es_l1[:-1][(Es_l1[:-1]<Es_l1[:-1][np.argmax(P_l1)]) & (P_l1!=0)][:int(4*size_l1/6)], (P_l1[(Es_l1[:-1]<Es_l1[:-1][np.argmax(P_l1)]) & (P_l1!=0)])[:int(4*size_l1/6)], color = my_blue, alpha = .8, ls = '', marker = 'o', ms = 4)
			ax.plot(Es_c1[:-1][(Es_c1[:-1]<Es_c1[:-1][np.argmax(P_c1)]) & (P_c1!=0)][:int(4*size_c1/6)], (P_c1[(Es_c1[:-1]<Es_c1[:-1][np.argmax(P_c1)]) & (P_c1!=0)])[:int(4*size_c1/6)], color = my_blue, alpha = .8, ls = '', marker = 's', ms = 4)	
			ax.plot(Es_l1[:-1], P_l1, color = my_blue, alpha = .6, ls = '', marker = 'o', ms = 4)
			ax.plot(Es_c1[:-1], P_c1, color = my_blue, alpha = .6, ls = '', marker = 's', ms = 4)	

			ax.plot(Es_l2_memory[:-1][(Es_l2_memory[:-1]<Es_l2_memory[:-1][np.argmax(P_l2_memory)]) & (P_l2_memory!=0)][:int(4*size_l2_memory/6)], (P_l2_memory[(Es_l2_memory[:-1]<Es_l2_memory[:-1][np.argmax(P_l2_memory)]) & (P_l2_memory!=0)])[:int(4*size_l2_memory/6)], color = my_red, alpha = .8, ls = '', marker = 'o', ms = 4)
			ax.plot(Es_c2_memory[:-1][(Es_c2_memory[:-1]<Es_c2_memory[:-1][np.argmax(P_c2_memory)]) & (P_c2_memory!=0)][:int(4*size_c2_memory/6)], (P_c2_memory[(Es_c2_memory[:-1]<Es_c2_memory[:-1][np.argmax(P_c2_memory)]) & (P_c2_memory!=0)])[:int(4*size_c2_memory/6)], color = my_red, alpha = .8, ls = '', marker = 's', ms = 4)	
			ax.plot(Es_l2_memory[:-1], P_l2_memory, color = my_red, alpha = .6, ls = '', marker = 'o', ms = 4)
			ax.plot(Es_c2_memory[:-1], P_c2_memory, color = my_red, alpha = .6, ls = '', marker = 's', ms = 4)	

			# my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
			# # ax.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
			# ax.set_xlim(left = -22, right = 0)
			# ax.set_ylim(bottom = 1e-5)
			# #ax.set_yticks([1, 0.1, 0.01, 0.001])
			# #ax.set_yticklabels([1, 0.1, 0.01])
			# fig.savefig('../../../Figures/memory/statistics/densities_'+energy_model+'_%d.pdf'%a)

			# ax.plot(Es[:-1], np.exp(beta_r * Es[:-1]) * (Q_R12[Es[:-2]<(E_r)][-1]) / (np.exp(beta_r * (E_r))), color = my_blue, ls = '--', lw = '3', alpha = 1)
			# ax.plot(Es[:-1], np.exp((beta_r - (lamB*p)/(lamA*1)) * Es[:-1]) * (Q_R12[Es[:-2]<(E_r+2)][-1]) / (np.exp((beta_r - (lamB*p)/(lamA*1)) * (E_r+2))), color = my_blue, ls = ':', lw = '3', alpha = 1)

			ax.plot(Es_l1, np.exp(popt_l_1[1] * Es_l1) * (P_l1[(Es_l1[:-1]<Es_l1[:-1][np.argmax(P_l1)]) & (P_l1!=0)])[int(2*size_l1/5)] / (np.exp(popt_l_1[1] * (Es_l1[int(2*size_l1/5)]))), color = my_blue, ls = '--', lw = '3', alpha = .6)
			ax.plot(Es_c1, np.exp(popt_c_1[1] * Es_c1) * (P_c1[(Es_c1[:-1]<Es_c1[:-1][np.argmax(P_c1)]) & (P_c1!=0)])[int(2*size_c1/5)] / (np.exp(popt_c_1[1] * (Es_c1[int(2*size_c1/5)]))), color = my_blue, ls = ':', lw = '3', alpha = .6)

			ax.plot(Es_l2_memory, np.exp(popt_l_2_memory[1] * Es_l2_memory) * (P_l2_memory[(Es_l2_memory[:-1]<Es_l2_memory[:-1][np.argmax(P_l2_memory)]) & (P_l2_memory!=0)])[int(2*size_l2_memory/5)] / (np.exp(popt_l_2_memory[1] * (Es_l2_memory[int(2*size_l2_memory/5)]))), color = my_red, ls = '--', lw = '3', alpha = .6)
			ax.plot(Es_c2_memory, np.exp(popt_c_2_memory[1] * Es_c2_memory) * (P_c2_memory[(Es_c2_memory[:-1]<Es_c2_memory[:-1][np.argmax(P_c2_memory)]) & (P_c2_memory!=0)])[int(2*size_c2_memory/5)] / (np.exp(popt_c_2_memory[1] * (Es_c2_memory[int(2*size_c2_memory/5)]))), color = my_red, ls = ':', lw = '3', alpha = .6)

			my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 22, x_fontsize=22, y_fontsize=22 )
			# ax.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
			ax.set_xlim(left = -22, right = 0)
			ax.set_ylim(bottom = 1e-3, top = 2e0)
			#ax.set_yticks([1, 0.1, 0.01, 0.001])
			#ax.set_yticklabels([1, 0.1, 0.01])
			fig.savefig('../../../Figures/memory/statistics/'+str(L0)+'/'+str(N_ens)+'/densities_'+energy_model+'_log_%d.pdf'%a)
			plt.close()

		# print('\n')
	
	x_l = np.linspace(0.2, 3.2, 50)
	fig_l, ax_l = plt.subplots(figsize=(5, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.12, 'top': 0.9})
	ax_l.scatter(beta_star_l_th_1, beta_star_l_sim_1, color = my_purple)
	ax_l.plot(x_l, x_l, color = 'k')
	my_plot_layout(ax = ax_l, xscale='linear', yscale= 'linear', ticks_labelsize= 22, x_fontsize=22, y_fontsize=22 )
	# ax_l.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
	ax_l.set_xlim(left = 0.4, right = 3.2)
	ax_l.set_ylim(bottom = 0.4, top = 3.2)
	#ax_l.set_yticks([1, 0.1, 0.01, 0.001])
	#ax_l.set_yticklabels([1, 0.1, 0.01])
	fig_l.savefig('../../../Figures/memory/statistics/'+str(L0)+'/'+str(N_ens)+'/exp_lineages_'+energy_model+'.pdf')
	plt.close()

	x_c = np.linspace(.2, 3.2, 50)
	fig_c, ax_c = plt.subplots(figsize=(5, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.12, 'top': 0.9})
	ax_c.scatter(beta_star_c_th_1, beta_star_c_sim_1, color = my_purple)
	ax_c.plot(x_c, x_c, color = 'k')
	my_plot_layout(ax = ax_c, xscale='linear', yscale= 'linear', ticks_labelsize= 22, x_fontsize=22, y_fontsize=22 )
	# ax_c.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
	ax_c.set_xlim(left = 0.4, right = 3.2)
	ax_c.set_ylim(bottom = 0.4, top = 3.2)
	#ax_c.set_yticks([1, 0.1, 0.01, 0.001])
	#ax_c.set_yticklabels([1, 0.1, 0.01])
	fig_c.savefig('../../../Figures/memory/statistics/'+str(L0)+'/'+str(N_ens)+'/exp_cells_'+energy_model+'.pdf')
	plt.close()

	x_sim = np.linspace(0.2, 3.4, 50)
	fig_sim, ax_sim = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	ax_sim.scatter(beta_star_l_sim_1, beta_star_c_sim_1, color = my_purple)
	ax_sim.plot(x_sim, x_sim - (lamB*p)/(lamA*1), color = 'k', ls = '--')
	ax_sim.plot(x_sim, x_sim, color = 'k')
	my_plot_layout(ax = ax_sim, xscale='linear', yscale= 'linear', ticks_labelsize= 22, x_fontsize=22, y_fontsize=22 )
	# ax_sim.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
	ax_sim.set_xlim(left = 0, right = 3.5)
	ax_sim.set_ylim(bottom = 0, top = 3.5)
	#ax_sim.set_yticks([1, 0.1, 0.01, 0.001])
	#ax_sim.set_yticklabels([1, 0.1, 0.01])
	fig_sim.savefig('../../../Figures/memory/statistics/'+str(L0)+'/'+str(N_ens)+'/exp_sim_l1c1_'+energy_model+'.pdf')
	plt.close()

	x_sim = np.linspace(.2, 4.4, 50)
	fig_sim, ax_sim = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	ax_sim.scatter(beta_star_c_sim_1, beta_star_l_sim_2, color = my_purple)
	ax_sim.scatter(beta_star_c_sim_1, beta_star_l_sim_2_memory, marker = 's', color = my_green)
	ax_sim.plot(x_sim, x_sim*(x_sim + (lamB*p)/(lamA*1)) , color = 'k', ls = '--')
	ax_sim.plot(x_sim, x_sim , color = 'k')
	my_plot_layout(ax = ax_sim, xscale='linear', yscale= 'linear', ticks_labelsize= 22, x_fontsize=22, y_fontsize=22 )
	# ax_sim.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
	ax_sim.set_xlim(left = 0, right = 4.5)
	ax_sim.set_ylim(bottom = 0, top = 4.5)
	#ax_sim.set_yticks([1, 0.1, 0.01, 0.001])
	#ax_sim.set_yticklabels([1, 0.1, 0.01])
	fig_sim.savefig('../../../Figures/memory/statistics/'+str(L0)+'/'+str(N_ens)+'/exp_sim_c1l2_'+energy_model+'.pdf')
	plt.close()

	x_sim = np.linspace(.2, 4.4, 50)
	fig_sim, ax_sim = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	ax_sim.scatter(beta_star_l_sim_2, beta_star_c_sim_2, color = my_purple)
	ax_sim.scatter(beta_star_l_sim_2_memory, beta_star_c_sim_2_memory, color = my_green)
	ax_sim.plot(x_sim, x_sim - (lamB*p)/(lamA*1) , color = 'k', ls = '--')
	ax_sim.plot(x_sim, x_sim , color = 'k')
	my_plot_layout(ax = ax_sim, xscale='linear', yscale= 'linear', ticks_labelsize= 22, x_fontsize=22, y_fontsize=22 )
	# ax_sim.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
	ax_sim.set_xlim(left = 0, right = 4.5)
	ax_sim.set_ylim(bottom = 0, top = 4.5)
	#ax_sim.set_yticks([1, 0.1, 0.01, 0.001])
	#ax_sim.set_yticklabels([1, 0.1, 0.01])
	fig_sim.savefig('../../../Figures/memory/statistics/'+str(L0)+'/'+str(N_ens)+'/exp_sim_l2c2_'+energy_model+'.pdf')
	plt.close()

	# D_temp1 = []
	# ts = np.linspace(t_lim-3.0, t_lim, 30)
	# for i_t, t in enumerate(ts):
	# 	u_on, p_a, R, QR1 = calculate_QR(Q0, k_on, k_step, np.exp(lamA*(t))/N_A, Es, p, lamA, b0, dE)
	# 	Q_R1 = QR1[::n_coarse][:-1]#*L0#/len(energies_lineages1)
	# 	Q_R1 = Q_R1/np.sum(Q_R1*dE[::n_coarse][:-1])
	# 	#D = abs(np.sum((np.cumsum(P_l1[P_l1!=0]))*np.log((np.cumsum(P_l1[P_l1!=0]))/(np.cumsum(Q_R1[P_l1!=0]))) * dE[::n_coarse][:-1][P_l1!=0]))
	# 	D = abs(np.sum(((P_l1[P_l1!=0]))*np.log(((P_l1[P_l1!=0]))/((Q_R1[P_l1!=0]))) * dE[::n_coarse][:-1][P_l1!=0]))
	# 	D_temp1.append(D)

	# t_optimal1 = ts[np.argmin(D_temp1)]
	# print(t_optimal1, np.mean(final_times1))

	# u_on, p_a, R, QR1 = calculate_QR(Q0, k_on, k_step, np.exp(lamA*(t_optimal1))/N_A, Es, p, lamA, b0, dE)
	# Q_R1 = QR1[::n_coarse][:-1]#*L0#/len(energies_lineages1)
	# Q_R1 = Q_R1/np.sum(Q_R1*dE[::n_coarse][:-1])

	# D_temp2 = []
	# ts = np.linspace(t_lim-3.0, t_lim, 30)
	# for i_t, t in enumerate(ts):
	# 	u_on, p_a, R, QR2 = calculate_QR(Q0, k_on, k_step, np.exp(lamA*(t))/N_A, Es, p, lamA, b0, dE)
	# 	Q_R2 = QR2[::n_coarse][:-1]#*L0#/len(energies_lineages2)
	# 	Q_R2 = Q_R2/np.sum(Q_R2*dE[::n_coarse][:-1])
	# 	#D = abs(np.sum((np.cumsum(P_l2[P_l2!=0]))*np.log((np.cumsum(P_l2[P_l2!=0]))/(np.cumsum(Q_R2[P_l2!=0]))) * dE[::n_coarse][:-1][P_l2!=0]))
	# 	D = abs(np.sum(((P_l2[P_l2!=0]))*np.log(((P_l2[P_l2!=0]))/((Q_R2[P_l2!=0]))) * dE[::n_coarse][:-1][P_l2!=0]))
	# 	D_temp2.append(D)

	# t_optimal2 = ts[np.argmin(D_temp2)]
	# u_on, p_a, R, QR2 = calculate_QR(Q0, k_on, k_step, np.exp(lamA*(t_optimal2))/N_A, Es, p, lamA, b0, dE)
	# Q_R2 = QR2[::n_coarse][:-1]#*L0#/len(energies_lineages2)
	# Q_R2 = Q_R2/np.sum(Q_R2*dE[::n_coarse][:-1])

if __name__ == "__main__":
    main()
