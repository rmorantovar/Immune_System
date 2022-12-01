import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'
#For simulation in C++
# ./EF_response_v2.x -a 6 -b 0.5 -k 1 -t 0 -T 6.5 -E MJ -C 10000 -B 100000000 -s TACNSEYPNTTKCGRWYC -q 2.0
# ./EF_response_v3.x -a 6 -b 0.5 -k 1 -t 0 -T 7.5 -E TCRen -C 10000 -B 100000000 -s EYTACNSEYPNTTKCGRWYCGRYPN -q 2.0 --ensemble -N 10
#--------------- PARAMETERS ---------------------
N_ens = 1
N_r = 1e8
T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_pr = 1/(60*5) #s^-1
k_pr = k_pr*3600 # hour^-1
k_pr = k_pr*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*1000
#N_c = 1e5
#E_ms = -27.63
E_ms = -25
C = 1e4
AA = 1

kappas = [1, 3]
t_fs = [2.48, 4.15]

my_red = np.array((228,75,41))/256.
my_purple = np.array((125,64,119))/256.
my_purple2 = np.array((116,97,164))/256.
my_green = np.array((125,165,38))/256.
my_blue = np.array((76,109,166))/256.
my_gold = np.array((215,139,45))/256.
my_brown = np.array((182,90,36))/256.
my_blue2 = np.array((80,141,188))/256.
my_yellow = np.array((246,181,56))/256.
my_green2 = np.array((158,248,72))/256.
my_cyan = 'tab:cyan'

antigen_color = my_yellow/256.

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_red, my_blue2, my_green, my_gold, my_brown])

#colors_kappa = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_kappa = np.flip(['tab:blue','tab:green','tab:red'])
colors_kappa = []
for i in range(len(color_list)):
        colors_kappa.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(kappas)):
    colors_R.append(['tab:grey', colors_kappa[i], colors_kappa[i], colors_kappa[i]])

time = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

# antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
# antigen = 'FMLFMAVFVMTSWYC'
# antigen = 'FTSENAYCGR'
# antigen = 'TACNSEYPNTTK'
#antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN'
#antigen = 'TACNSEYPNTTKCGRWYC'
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20

L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, Text_files_path)
print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

#--------------------------Repertoire properties--------------------------
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, N_r)
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')

min_E = -17.8
min_E = -19.5
max_E = -8

fig, ax = plt.subplots(figsize = (12, 1), linewidth = 6, gridspec_kw={'left':0.1, 'right':.9, 'bottom':.01, 'top': .3})
col_map = 'cividis'
#mpl.colorbar.ColorbarBase(ax, cmap=col_map, orientation = 'vertical')
# and create another colorbar with:
colorbar = mpl.colorbar.ColorbarBase(ax, cmap=plt.get_cmap(col_map + '_r'), orientation = 'horizontal')
colorbar.set_ticks(np.linspace(0, 1, 5))
#ax.set_xticks(np.linspace(0, 1, 5))
ax.set_xticklabels([r'$%.0f \cdot 10^{%d}$'%(10**(np.log10(np.exp(min_E + i*(max_E - min_E)/4))%1), int(np.log10(np.exp(min_E + i*(max_E - min_E)/4)))) for i in np.arange(0, 5, 1)], fontsize = 38)
ax.xaxis.tick_top()
fig.savefig("../../Figures/1_Dynamics/Trajectories/Muller/colorbar.pdf")

print('Loops...')
#--------------------------Loops--------------------------
for i_kappa, kappa in enumerate(kappas):
	print('--------')
	print('kappa = %.2f...'%kappa)
	beta_kappa, E_kappa, Kd_kappa = get_kappa_properties(betas, Q0, Es, dE, kappa)
	beta_act = np.min([beta_r, beta_kappa])
	m_bar_theory = np.array([np.sum(N_r*calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t))/N_A, Es, kappa, lambda_A, N_c, dE)[3]*dE) for t in time])
	t_act_theory = time[m_bar_theory>1][0]
	for rep in [0, 1, 2, 3]:
		#fig_muller, ax_muller = plt.subplots(figsize=(9/1.5,4/1.5), linewidth = 0, gridspec_kw={'left':0.005, 'right':.995, 'bottom':.02, 'top': 0.98}, dpi = 700, edgecolor = 'black')
		fig_muller, ax_muller = plt.subplots(figsize=(4.5,3*0.8), linewidth = 0, gridspec_kw={'left':0.005, 'right':.995, 'bottom':.02, 'top': 0.98}, dpi = 700, edgecolor = 'black')
		ax_muller.spines["top"].set_linewidth(3)
		ax_muller.spines["left"].set_linewidth(3)
		ax_muller.spines["right"].set_linewidth(3)
		ax_muller.spines["bottom"].set_linewidth(3)
		#-----------------Loading data----------------------------
		parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_pr/24, kappa, N_c, linear, N_ens)+energy_model
		#data = pd.read_csv(Text_files_path + 'Dynamics/Trajectories/'+parameters_path+'/energies%d.txt'%rep, sep = '\t', header=None)
		data = get_data(folder_path = Text_files_path + 'Dynamics/Trajectories/L%d/'%L+parameters_path, rep = rep)
		#-----------------Filtering data----------------------------
		min_e_data = np.min(data[0])
		max_e_data = np.max(data[0])

		data_active = data.loc[data[1]==1]
		print('Activated clones before filter:%d'%len(data_active[0]))
		t_act_data = np.min(data_active[3])
		print('t_act_data: %.2f'%t_act_data)
		data_active = data_active.loc[data_active[3]<(t_act_data+1.0+0.1*(kappa-1))]
		activation_times = np.array(data_active[3])
		print('Activated clones after filter:%d'%len(activation_times))
		energies  = np.array(data_active[0])
		ar1, ar2 = np.histogram(activation_times, bins = time)
		m_data = np.cumsum(ar1)
		#---------------------------- B cell linages ----------------------
		clone_sizes = get_clones_sizes_C(int(m_data[-1]), time, activation_times, lambda_B, C, dT)
		#-----------------------------Activation time------------------------
		t_act = get_t_act(time, N_r, Q0, k_on, k_pr, lambda_A, Es, dE, kappa, N_c)
		#-----------------------------m(t)-----------------------------------
		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf_sim)/N_A, Es, kappa, lambda_A, N_c, dE)
		m_f_expected = np.sum(N_r*QR*dE)
		print('Activated clones expected at Tf_sim:%.d'%m_f_expected)

		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*(t_act_data+1.3))/N_A, Es, kappa, lambda_A, N_c, dE)
		m_f_expected = np.sum(N_r*QR*dE)
		print('Activated clones expected at t_act_data + 1.3:%.d'%m_f_expected)

		#---------------------- Total Pop size ----------------------
		total_pop = np.sum(clone_sizes, axis = 0)
		total_pop_active = total_pop - total_pop[0] + 1
		t_C = time[total_pop_active<(C-1)][-1] # Calculate time for reaching carrying capacity
		#--------------------------t_C filter-------------------------
		print('Applying filter...')
		lim_size = 2
		clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
		ax_muller.vlines([t_act_theory, t_fs[i_kappa]], 0, 1, color = 'black', linewidth = 2, alpha = 1, linestyle = '--')
		print('min time : %.2f'%np.min(activation_times_C))
		#ax_muller.vlines(np.max(activation_times_C), 0, 1, color = 'black', linewidth = 2, alpha = .8, linestyle = ':')
		print('max time : %.2f'%np.max(activation_times_C))
		print('Applying filter...')
		lim_size = 20
		clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
		print('Activated clones:', np.shape(clone_sizes_C))
		sort_inds = activation_times_C.argsort()
		total_pop_active = np.sum(clone_sizes_C, axis = 0) #C
		bcell_freqs = clone_sizes_C/total_pop_active
		bcell_freqs = clone_sizes_C/(np.sum(clone_sizes_C[:,-1]))
		entropy = -np.sum(bcell_freqs*np.log(bcell_freqs), axis = 0)

		#------------------------- Stackplots -------------------------
		greys = plt.cm.get_cmap(col_map, 50)
		min_bell_freq = np.min(bcell_freqs[:,-1])
		
		delta_E = max_E - min_E
		counter_final = 0
		for c in np.flip(range(len(clone_sizes_C[:,0]))):
			color_c = greys(int(50*(1-abs((energies_C[c]-min_E)/delta_E))))
			ax_muller.stackplot(time, [(bcell_freqs[c, -1] - bcell_freqs[c, :])/2 + np.ones_like(bcell_freqs[0, :])*np.sum(bcell_freqs[:c, -1]), bcell_freqs[c, :], (bcell_freqs[c, -1] - bcell_freqs[c, :])/2], colors = ['white', color_c, 'white']);
			if activation_times_C[c] in activation_times_C[sort_inds[:3]]:
				ax_muller.scatter(activation_times_C[c], (bcell_freqs[c, -1] - bcell_freqs[c, 0])/2 + np.sum(bcell_freqs[:c, -1]), marker = 'D', edgecolor='black', linewidth=1, facecolor = 'white', s = 60, zorder = 20)
			# if (activation_times_C[c] in activation_times_C[sort_inds[-1:]]):
			# 	if(np.random.randint(0, 2)>0 and counter_final==0):
			# 		ax_muller.scatter(activation_times_C[c], (bcell_freqs[c, -1] - bcell_freqs[c, 0])/2 + np.sum(bcell_freqs[:c, -1]), marker = 'o', edgecolor='black', linewidth=1, facecolor = 'white', s = 60, zorder = 20)
			# 		counter_final+=1

		cumsum_freqs = np.cumsum(bcell_freqs, axis = 0)

		my_plot_layout(ax = ax_muller, ticks_labelsize=38, yscale = 'linear')
		ax_muller.set_yticks([])
		#ax_muller.set_xticks(np.arange(Tf))
		ax_muller.set_xticks([])
		ax_muller.set_xlim(T0, Tf-1)
		ax_muller.set_ylim(0, 1)
		fig_muller.savefig('../../Figures/1_Dynamics/Trajectories/Muller/B_cell_clones_kappa-%.2f_%d_'%(kappa, rep)+energy_model+'.pdf', edgecolor=fig_muller.get_edgecolor())
		plt.close(fig_muller)



