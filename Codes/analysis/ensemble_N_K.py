import sys
sys.path.append('../library/')
from Immuno_models import*
plt.rcParams['text.usetex'] = True


N_A = 6.02214076e23
colors = ['tab:blue','tab:red']
colors_fit = ['darkblue', 'darkred']
growth_models = ['exponential']#, 'linear']

Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/"

M1 = np.loadtxt(Text_files_path+'MJ.txt', skiprows= 1, usecols=range(1,21)).tolist()
M2 = (np.loadtxt(Text_files_path+'MJ2.txt', skiprows= 1, usecols=range(1,21))).tolist()
M3 = np.loadtxt(Text_files_path+'BLOSUM62.txt', skiprows= 1, max_rows = 23, usecols=range(1,24)).tolist()
Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
Alphabet2 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w']
Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)

Matrix = 'MJ2'
#Matrix = 'MM'

if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()

antigen = 'CMFILVWYAGTSQNEDHRKPFMRTP'
antigen = 'FMLFMAVFVMTSWYC'
antigen = 'TACNSEYPNTTK'
antigen = 'TACNSEYPNTTKCGRWYC'
#antigen = 'FTSENAYCGR'
#antigen = 'MRTAYRNG'
#antigen = 'MRTAY'

L=len(antigen)

N_ens = 100
N_r = 5e4
N_r = 1e6
#N_r = 1e6
T0 = 0
Tf = 6
#Tf = 8
dT = .1
days = np.arange(0, Tf, 1)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_pr = .1 # hour^-1
k_pr = k_pr*24 #days^-1
thetas = [1, 1.5, 2]
colors_theta = ['darkred', 'olive', 'navy']
colors_theta = ['darkred', 'olive', 'darkblue']
lambda_B = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e4
E_ms = -28
C = 2e4

print('k_pr/k_on = %.1e'%(k_on/k_pr)**(-1))

#----------------------------------------------------------------

antigen_list = [i for i in antigen]
antigen_seq = np.array([], dtype = int)
for i, aa in enumerate(antigen_list):
    index = Alphabet_list.index(aa)
    antigen_seq = np.append(antigen_seq, int(index))
PWM_data = M2[:,antigen_seq]

#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

Es, dE, Q0, betas = calculate_Q0(0.01, 50, 200000, PWM_data, E_ms, L)
Kds = np.exp(Es[:-1])

beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
Kd_r = np.exp(E_r)

print('beta_r = %.2f'%beta_r)


E_pr = Es[:-1][Ks<(k_pr/k_on)][-1]
Kd_pr = np.exp(E_pr)
beta_pr = betas[Ks<Kd_pr][-1]
print('beta_pr = %.2f'%beta_pr)
#----------------------------------------------------------------

lambda_Bs = np.array([np.flip([.5])*lambda_A, np.flip([.5])*lambda_A], dtype=object)
n_bins = 15
d=20
energy_model = 'MJ'
colors_gm = np.array([plt.cm.Oranges(np.linspace(0,1,len(lambda_Bs[0])+2)),plt.cm.Reds(np.linspace(0,1,len(lambda_Bs[1])+2)) ], dtype=object)
FIG, AX = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
FIG2, AX2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
FIG3, AX3 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

delta_Nb = lambda t, tb, Nb, N, lambda_B, C: lambda_B*Nb*(1-(N/C))*np.heaviside(t-tb, 1)

for i_theta, theta in enumerate(thetas):
	print('theta = %.2f...'%theta)
	beta_q = betas[betas>theta][-1]
	E_q = Es[betas>theta][-1]
	Kd_q = np.exp(E_q)
	Kd_act = np.max([Kd_q, Kd_r])
	for j, gm in enumerate(growth_models):
		fig0, ax0 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig, ax = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		fig2, ax2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
		for n_lambda_B, lambda_B in enumerate(lambda_Bs[j]):

			#--------------------------m(t)---------------------------
			u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time[0])/N_A, Es, theta, lambda_A, N_c, dE)
			M_r = N_r*N_c*np.sum(Q0*p_a*dE)
			m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time])

			t_act = time[m_bar>1][0]
			t_C = t_act+1.2
			#--------------------------------------------------------

			parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, N_r)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_pr-%.6f_theta-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 0.5, k_pr/24, theta, j, N_ens)+energy_model
			data = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/energies_ensemble.txt', sep = '\t', header=None)
			data2 = pd.read_csv(Text_files_path + 'Dynamics/Ensemble/'+parameters_path+'/summary_ensemble.txt', sep = '\t', header=None)

			min_e_data = np.min(data[0])
			max_e_data = np.max(data[0])

			energies0 = data[0]
			data_Es0 = ax0.hist(energies0, bins = np.linspace(min_e_data, max_e_data, n_bins), density = False, color = colors_gm[j][n_lambda_B+2],histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
			counts0 = data_Es0[0]
			Es_array_data0 = data_Es0[1][:-1][counts0!=0]
			counts0 = counts0[counts0!=0]

			ax.plot(np.exp(Es_array_data0), counts0/N_ens, color = colors_gm[j][n_lambda_B+2], alpha = .8, marker = 'o', ms = 8, linestyle = '')
			#popt, pcov = curve_fit(f = my_linear_func , xdata = np.log(Es_array_data0[0:4]), ydata= np.log(counts0)[0:4] )
			#beta_act2 = popt[1]
			beta_act = np.min([theta, beta_r])

			print('beta_q = %.2f'%beta_q)
			print('beta_act = %.2f'%(beta_act))

			u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*Tf)/N_A, Es, theta, lambda_A, N_c, dE)

			ax.plot(np.exp(Es[:-1]), Q0*N_r, color = colors_gm[j][n_lambda_B+2])
			ax.plot(np.exp(Es[:-1]), QR*N_r, color = colors_gm[j][n_lambda_B+2])
			ax.plot(np.exp(Es_array_data0[0:10]), (counts0[0]/N_ens)*(np.exp(Es_array_data0[0:10])/np.exp(Es_array_data0[0]))**(beta_act), color = colors_gm[j][n_lambda_B+2], linewidth = 5, linestyle = ':', marker = '', ms = 15, alpha = .8)
			#ax.plot(np.exp(Es_array_data0[0:10]), (counts0[0]/N_ens)*(np.exp(Es_array_data0[0:10])/np.exp(Es_array_data0[0]))**(beta_act2), color = colors_gm[j][n_lambda_B+2], linewidth = 5, linestyle = ':', marker = '', ms = 15, alpha = .4)
			
			clone_sizes_binned = np.zeros(n_bins-1)
			var_clone_sizes_binned = np.zeros(n_bins-1)
			max_clone_sizes_binned = np.zeros(n_bins-1)

			for i_ens in tqdm(np.arange(N_ens)):

				data_i = data.loc[data[4]==i_ens]
				data_active = data_i.loc[data_i[1]==1]

				data_active = data_active.loc[data_active[3]<t_C]
				energies = data_active[0]
				data_plasma = data_active.loc[data_active[2]==1]
				data_GC = data_active.loc[data_active[2]==0]
				activation_times = np.array(data_active[3]) #Can be changed for data_plasma
				
				#clone_sizes = np.exp(lambda_B*(Tf - activation_times))
				clone_sizes = np.ones((len(activation_times), len(time)))

				# #-----t_C filter-------
				filter_C = activation_times<t_C
				clone_sizes = clone_sizes[filter_C, :]
				activation_times = activation_times[filter_C]
				energies = energies[filter_C]

				for i_t, t in enumerate(time[:-1]):
					for i in np.arange(0, len(activation_times)):
						tb = activation_times[i]
						Nb = clone_sizes[i, i_t]# * np.heaviside(tb-t)
						N = np.sum(clone_sizes[:, i_t]) - np.sum(clone_sizes[:, 0])
						clone_sizes[i, i_t+1] = Nb + delta_Nb(t, tb, Nb, N, lambda_B, C)*dT

				clone_sizes_f = clone_sizes[:,-1]

				Kds = np.exp(data_active[0]) #Can be changed for data_plasma

				data_Es = ax0.hist(energies, bins = np.linspace(min_e_data, max_e_data, n_bins), density = False, color = colors_gm[j][n_lambda_B+2],histtype = 'step', zorder=10, align = 'mid', linewidth = 2, alpha = 0)
				counts = data_Es[0]
				Es_array_data = data_Es[1][:-1]#[counts!=0]
				#counts = counts[counts!=0]

				clone_sizes_binned_i = np.zeros([len(counts)])
				act_times_binned = np.zeros([len(counts)])
				clone_sizes_binned_2= np.zeros([len(counts)])
				var_clone_sizes_binned_i = np.zeros([len(counts)])
				max_clone_sizes_binned_i = np.zeros([len(counts)])

				#print(int(len(counts)), len(Es_array_data))
				for i in np.arange(int(len(counts))-1):
						clone_sizes_binned_i[i] = np.mean( np.concatenate((clone_sizes_f[(energies>=Es_array_data[i]) & (energies<Es_array_data[i+1])], np.array([0]) )) )#/N_ens
						var_clone_sizes_binned_i[i] = np.var( np.concatenate((clone_sizes_f[(energies>=Es_array_data[i]) & (energies<Es_array_data[i+1])], np.array([0]) )) )#/N_ens
						#act_times_binned[i] = np.mean( np.concatenate((activation_times[(energies>=Es_array_data[i]) & (energies<Es_array_data[i+1])], np.array([0]) )) )#/N_ens
						#clone_sizes_binned_2[i] += np.mean( np.concatenate(((clone_sizes[(energies>=Es_array_data[i]) & (energies<Es_array_data[i+1])])**2, np.array([0]) )) )#/N_ens
						max_clone_sizes_binned_i[i] = np.max(clone_sizes_f[(energies>=Es_array_data[i]) & (energies<Es_array_data[i+1]) ], initial=1)#/N_ens

				clone_sizes_binned += clone_sizes_binned_i/N_ens
				var_clone_sizes_binned += var_clone_sizes_binned_i/N_ens
				max_clone_sizes_binned += max_clone_sizes_binned_i/N_ens

				max_clone_size = 1# np.max(clone_sizes_binned_i)
				
				ax.plot(np.exp(Es_array_data), counts, color = colors_gm[j][n_lambda_B+2], alpha = .8, marker = '^', ms = 5, linestyle = '')
			
				#Es_array_data = Es_array_data[clone_sizes_binned_i!=0]
				#act_times_binned = act_times_binned[clone_sizes_binned_i!=0]
				#max_clone_sizes_binned_i = max_clone_sizes_binned_i[clone_sizes_binned_i!=0]
				#clone_sizes_binned_i = clone_sizes_binned_i[clone_sizes_binned_i!=0]

				#-------Simulations-------
				ax2.plot(np.exp(Es_array_data[:]), clone_sizes_binned_i[:]/max_clone_size, color = 'orange', linewidth =5, linestyle = '', marker = 's', ms = 5, alpha = .05)
				#ax2.plot(np.exp(Es_array_data[:]), max_clone_sizes_binned_i[:]/max_clone_size, linewidth =5, linestyle = '', marker = '*', ms = 8,  color = colors_theta[i_theta])
				AX.plot(np.exp(Es_array_data[:]), clone_sizes_binned_i[:]/max_clone_size, linewidth = 5, linestyle = '', marker = 's', ms = 5, color = colors_theta[i_theta], alpha = .05)
				#AX.plot(np.exp(Es_array_data[:]), max_clone_sizes_binned_i[:]/max_clone_size, linewidth =5, linestyle = '', marker = '*', ms = 8,  color = colors_theta[i_theta])
				AX.vlines([Kd_pr, Kd_q, Kd_r], 1, 1e2, linestyles = ['-',':', '--'], color = ['grey', colors_theta[i_theta], 'gray'])

				AX2.scatter(np.exp(energies), activation_times, color = colors_theta[i_theta])
				AX2.vlines([Kd_pr, Kd_q, Kd_r], 4, 6, linestyles = ['-',':', '--'], color = ['grey', colors_theta[i_theta], 'gray'])

				AX3.scatter(np.exp(energies), clone_sizes_f, color = colors_theta[i_theta])
				AX3.vlines([Kd_pr, Kd_q, Kd_r], 4, 6, linestyles = ['-',':', '--'], color = ['grey', colors_theta[i_theta], 'gray'])

			ax2.plot(np.exp(Es_array_data[:]), clone_sizes_binned[:]/max_clone_size, color = 'orange', linewidth =3, linestyle = '-', marker = 's', ms = 4, alpha = 1, label = r'$%.2f$'%theta)
			AX.plot(np.exp(Es_array_data[:]), clone_sizes_binned[:]/max_clone_size, linewidth =3, linestyle = '-', marker = 's', ms = 4, color = colors_theta[i_theta], alpha = 1, label = r'$%.2f$'%theta)

			#-------Theory-------
			cross_over = 0# np.where(clone_sizes_binned_i==max_clone_size)[0][0]
			#ax2.plot(np.exp(Es_array_data[0:cross_over+1]), (clone_sizes_binned_i[cross_over]/max_clone_size)*(np.exp(Es_array_data[0:cross_over+1])/np.exp(Es_array_data[cross_over]))**((lambda_B/lambda_A)), color = 'orange', linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)
			#ax2.plot(np.exp(Es_array_data[:]), (clone_sizes_binned_i[-1]/max_clone_size)*(np.exp(Es_array_data[:])/np.exp(Es_array_data[-1]))**((lambda_B/lambda_A)*(-theta)), color = 'orange', linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)
			ax2.vlines([Kd_pr, Kd_q, Kd_r], 4e-3, 1.5, linestyles = ['-',':', '--'], color = 'grey')

			#AX.plot(Kds_array_data[0:cross_over+1], (clone_sizes_binned_i[cross_over]/max_clone_size)*(Kds_array_data[0:cross_over+1]/Kds_array_data[cross_over])**((lambda_B/lambda_A)), linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8)
			#AX.plot(Kds_array_data[cross_over:], (clone_sizes_binned_i[cross_over]/max_clone_size)*(Kds_array_data[cross_over:]/Kds_array_data[cross_over])**((lambda_B/lambda_A)*(-theta)), linewidth =3, linestyle = '--', marker = '', ms = 15, alpha = .8, color = colors_theta[i_theta])
			#AX.vlines(Kd_q, AX.get_ylim()[0], AX.get_ylim()[1])

	#AX.vlines(Kd_r, AX.get_ylim()[0], AX.get_ylim()[1])
	my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	#ax.legend(title=r'$\lambda_A/\lambda_B$', fontsize = 30, title_fontsize = 35)
	ax.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
	ax.set_ylim(bottom = 1e-7)
	fig.savefig('../../Figures/1_Dynamics/Ensemble/Q_K_theta-%.1f.pdf'%theta)

	my_plot_layout(ax = ax2, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	#ax2.legend(fontsize = 30, title_fontsize = 35)
	ax2.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
	#ax2.set_ylim(bottom = 5e-5, top = 2e0)
	ax2.set_yticks([1, 0.1, 0.01, 0.001])
	#ax2.set_yticklabels([1, 0.1, 0.01, 0.001])
	fig2.savefig('../../Figures/1_Dynamics/Ensemble/N_vs_K_theta-%.1f.pdf'%theta)


my_plot_layout(ax = AX, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
AX.legend(fontsize = 30, title_fontsize = 35, title = r'$\theta$')
AX.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#AX.set_ylim(bottom = 5e-5, top = 2e0)
#AX.set_yticks([1, 0.1, 0.01, 0.001])
#AX.set_yticklabels([1, 0.1, 0.01])
FIG.savefig('../../Figures/1_Dynamics/Ensemble/N_vs_K.pdf')

my_plot_layout(ax = AX2, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#AX2.legend(fontsize = 30, title_fontsize = 35, title = r'$\theta$')
#AX2.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#AX2.set_ylim(bottom = 5e-5, top = 2e0)
#AX2.set_yticks([1, 0.1, 0.01, 0.001])
#AX2.set_yticklabels([1, 0.1, 0.01])
FIG2.savefig('../../Figures/1_Dynamics/Ensemble/act_T_vs_K.pdf')

my_plot_layout(ax = AX3, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
#AX3.legend(fontsize = 30, title_fontsize = 35, title = r'$\theta$')
#AX3.set_xlim(left = np.exp(E_ms+2), right = np.exp(E_ms+29))
#AX3.set_ylim(bottom = 5e-5, top = 2e0)
#AX3.set_yticks([1, 0.1, 0.01, 0.001])
#AX3.set_yticklabels([1, 0.1, 0.01])
FIG3.savefig('../../Figures/1_Dynamics/Ensemble/clone_sizes_vs_K.pdf')






