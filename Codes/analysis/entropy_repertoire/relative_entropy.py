import sys
sys.path.append('../library/')
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Immune_system/repertoire_entropy/H/'

#--------------- PARAMETERS ---------------------
L0s = np.logspace(4, 19, 16)
T0 = 0
Tf = 20
Tf_sim = 7
#Tf = 10
dT = 0.02
lambda_A = 6
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*10
E_m = -24
C = 1e4
AA = 1

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['MJ']
energy_model = 'MJ'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0

E_lims = [-5.0, -5.5, -6.0, -6.5, -7, -7.5, -8, -8.5, -9.0]
t_lims = [4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8.0]
ps = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5.0]
ps = np.linspace(1, 10, 80)
color_list = np.array([my_blue2, my_green, my_purple2, my_red, my_gold])


antigen = 'EYTACNSEYPNTTKCGRWYCGRYPN' #L=25
# antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
# antigen = 'RMTWPVAAERTMYYTPRGVA' #L=20
# antigen = 'MNTYRREAETPWWSYTPMGE' #L=20
#antigen = 'TACNSEYPNTTKCGRWYC' #L=18

L=len(antigen)
# print('--------')
# print('L=%d'%(L))
#----------------------------------------------------------------
energy_model = 'TCRen'
#energy_model = 'MJ2'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, energy_model, '../../in/')
print(Alphabet)

# print('min_e_PWM=%.2f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
# print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 5000, PWM_data, E_m, L)
Kds = np.exp(Es[:-1])

fig_D, ax_D = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_p, ax_p = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_t, ax_t = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_step, ax_step = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
fig_K_star, ax_K_star = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

for L0 in tqdm(L0s):
	# print('--------')
	# print('L0 = %.0e...'%L0)
		#--------------------------Repertoire properties--------------------------
	beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L0)
	# print('L0 = %.e'%L0, 'beta_r = %.1f'%beta_r, 'K_r = %.1e'%Kd_r, 'E_r = %.1f'%E_r)
	ps = [beta_r]
	k_step = Kd_r*k_on/2

	#--------------------------Proofreading properties--------------------------
	beta_step, E_step, Kd_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
	# print('beta_step = %.2f'%beta_step)


	t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
	# print('Loops...')
	#--------------------------Loops--------------------------
	#fig_K, ax_K = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})


	D_theory = {}
	# fig_loglike, ax_loglike = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

	for i_p, p in enumerate(ps):
		# fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

		m_bar_theory = np.array([np.sum(L0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)[3]*dE) for t in time_array])
		t_act_theory = time_array[m_bar_theory>1][0] 
		
		# print('p = %.2f...'%p)
		beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
		# print('K_p = %.1e...'%Kd_p, 'E_p = %.1e...'%E_p)
		# print('t_act_theoty:', t_act_theory)
		if p==1:
			t_act_1 = t_act_theory
			# print(t_act_1)
		

		n_coarse = 1
		bins = Es[::n_coarse]
		Q_0 = Q0[::n_coarse][:-1]#*L0#/len(energies_lineages)

		# ts = np.linspace(t_lim-3.0, t_lim, 30)
		# for i_t, t in enumerate(ts):
		# 	u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, N_c, dE)
		# 	Q_R = QR[::n_coarse][:-1]#*L0#/len(energies_lineages)
		# 	Q_R = Q_R/np.sum(Q_R*dE[::n_coarse][:-1])
		# 	#D = abs(np.sum((np.cumsum(P_l[P_l!=0]))*np.log((np.cumsum(P_l[P_l!=0]))/(np.cumsum(Q_R[P_l!=0]))) * dE[::n_coarse][:-1][P_l!=0]))


		u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t_act_theory))/N_A, Es, p, lambda_A, N_c, dE)
		Q_R = QR[::n_coarse][:-1]#*L0#/len(energies_lineages)
		Q_R = Q_R/np.sum(Q_R*dE[::n_coarse][:-1])
		
		# ax.plot(Es[::n_coarse][:-2][::1], Q_0[::1], color = 'grey', ls = '--', label = r'$\Omega_{0}$', lw = 2)
		# ax.plot(Es[::n_coarse][:-2][::1], Q_R[::1], color = 'black', ls = '--', label = r'$\Omega_{\textrm{act}}$')
		
		# ax.plot(Es[:-1], np.exp(beta_r * Es[:-1]) * (Q0[Es[:-1]<E_r][-1]*1e8/len(energies_lineages)) / (np.exp(beta_r * E_r)), color = my_blue, ls = '--')
		# ax.plot(Es[:-1], np.exp((beta_r - (lambda_B*p)/(lambda_A)) * Es[:-1]) * (Q0[Es[:-1]<(E_r+1)][-1]*1e8/len(energies_lineages)) / (np.exp((beta_r - (lambda_B*p)/(lambda_A)) * (E_r+1))), color = my_red, ls = ':')
		# ax.plot(Es[:-1], np.exp((-4+beta_r) * Es[:-1]) * (Q0[Es[:-1]<(E_r+1.4)][-1]*1e8/len(energies_lineages)) / (np.exp((-4+beta_r) * (E_r+1.4))), color = my_blue, ls = '--')

		# ax_loglike.plot(Es[::n_coarse][:-2], np.log10((Q_R)/(Q_0)), ls = '--', alpha = 1, label = r'$%.1f$'%p, lw = 2)

		D_theory[p] = np.sum((Q_R[Q_R!=0])*np.log10((Q_R[Q_R!=0])/(Q_0[Q_R!=0])) * dE[::n_coarse][:-1][Q_R!=0])# - np.sum((P_c[P_c!=0])*np.log((Q_R[P_c!=0])/(1)) * dE[::n_coarse][:-1][P_c!=0])
		
		# my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		# ax.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
		# ax.set_xlim(left = -22, right = -5+1+8)
		# ax.set_ylim(bottom = 2e-7, top = 8e-1)
		# #ax.set_yticks([1, 0.1, 0.01, 0.001])
		# #ax.set_yticklabels([1, 0.1, 0.01])
		# fig.savefig('../../../Figures/repertoire_entropy/theory/Pi_vs_Omega_act_'+energy_model+'_p-%.1f_log.pdf'%(p))

		# my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
		# ax.legend(fontsize = 14, title_fontsize = 24, title = r'$\textrm{level}$', loc = 0)
		# ax.set_xlim(left = -22, right = -5+1+8)
		# ax.set_ylim(bottom = 2e-7, top = np.max(Q_R)+.2)
		# #ax.set_yticks([1, 0.1, 0.01, 0.001])
		# #ax.set_yticklabels([1, 0.1, 0.01])
		# fig.savefig('../../../Figures/repertoire_entropy/theory/Pi_vs_Omega_act_'+energy_model+'_p-%.1f.pdf'%(p))

	# my_plot_layout(ax = ax_loglike, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	# ax_loglike.legend(fontsize = 16, title_fontsize = 20, loc = 0, title = r'$p$')
	# #ax_loglike.set_title(r'$\hat{\Pi}_{\textrm{act}}/\Omega_{\textrm{act}}$', fontsize = 22)
	# ax_loglike.set_xlim(left = -21, right = -5+1)
	# ax_loglike.set_ylim(bottom = 2e-5)
	# #ax_loglike.set_yticks([1, 0.1, 0.01, 0.001])
	# #ax_loglike.set_yticklabels([1, 0.1, 0.01])
	# fig_loglike.savefig('../../../Figures/repertoire_entropy/theory/Pi_vs_Omega_act_distance_'+energy_model+'.pdf')

	#ax_D.plot(ps, D_theory.values(), label = r'$%d$'%np.log10(L0), ls = '-')
	ax_D.scatter(L0, D_theory.values(), marker = 'D', label = r'$%.1f$ ; $%.1e$ ; $%.1f$'%(beta_r, k_step/(24*60), t_act_theory))
	ax_p.scatter(L0, beta_r, marker = 'o', color = my_blue)
	ax_t.scatter(L0, t_act_theory, marker = 'o', color = my_red)
	ax_step.scatter(L0, k_step/(24*60), marker = 'o', color = my_green)
	ax_K_star.scatter(L0, Kd_r, marker = 'o', color = my_purple)


my_plot_layout(ax = ax_D, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=24, y_fontsize=24)
ax_D.legend(fontsize = 14, title_fontsize = 16, loc = 0, title = r'$p$ ; $k_{\mathrm{step}}[\mathrm{ min}^{-1}]$ ; $t_{\mathrm{{\scriptsize act}}}[\mathrm{ days}]$')
ax_D.set_ylabel(r'$D({\Omega_{{\textrm{act}}}|\Omega_{0}})$')
#ax_D.set_xlim(left = -21, right = -7)
#ax_D.set_ylim(bottom = 2e-7, top = 3e-2)
ax_D.set_yticks(range(3, 16, 2))
#ax_D.set_yticklabels([1, 0.1, 0.01])
fig_D.savefig('../../../Figures/repertoire_entropy/theory/D_'+energy_model+'.pdf')

my_plot_layout(ax = ax_p, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=24, y_fontsize=24)
# ax_p.legend(fontsize = 14, title_fontsize = 16, loc = 0, title = r'$p$ ; $k_{\mathrm{step}}[\mathrm{ min}^{-1}]$ ; $t_{\mathrm{{\scriptsize act}}}[\mathrm{ days}]$')
#ax_p.set_ylabel(r'$D({\Omega_{{\textrm{act}}}|\Omega_{0}})$')
#ax_p.set_xlim(left = -21, right = -7)
#ax_p.set_ylim(bottom = 2e-7, top = 3e-2)
#ax_p.set_yticks(range(3, 16, 2))
#ax_p.set_yticklabels([1, 0.1, 0.01])
fig_p.savefig('../../../Figures/repertoire_entropy/theory/p_'+energy_model+'.pdf')

my_plot_layout(ax = ax_t, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=24, y_fontsize=24)
# ax_t.legend(fontsize = 14, title_fontsize = 16, loc = 0, title = r'$p$ ; $k_{\mathrm{step}}[\mathrm{ min}^{-1}]$ ; $t_{\mathrm{{\scriptsize act}}}[\mathrm{ days}]$')
#ax_t.set_ylabel(r'$D({\Omega_{{\textrm{act}}}|\Omega_{0}})$')
#ax_t.set_xlim(left = -21, right = -7)
#ax_t.set_ylim(bottom = 2e-7, top = 3e-2)
#ax_t.set_yticks(range(3, 16, 2))
#ax_t.set_yticklabels([1, 0.1, 0.01])
fig_t.savefig('../../../Figures/repertoire_entropy/theory/t_'+energy_model+'.pdf')

my_plot_layout(ax = ax_step, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=24, y_fontsize=24)
# ax_step.legend(fontsize = 14, title_fontsize = 16, loc = 0, title = r'$p$ ; $k_{\mathrm{step}}[\mathrm{ min}^{-1}]$ ; $t_{\mathrm{{\scriptsize act}}}[\mathrm{ days}]$')
#ax_step.set_ylabel(r'$D({\Omega_{{\textrm{act}}}|\Omega_{0}})$')
#ax_step.set_xlim(left = -21, right = -7)
#ax_step.set_ylim(bottom = 2e-7, top = 3e-2)
#ax_step.set_yticks(range(3, 16, 2))
#ax_step.set_yticklabels([1, 0.1, 0.01])
fig_step.savefig('../../../Figures/repertoire_entropy/theory/step_'+energy_model+'.pdf')

my_plot_layout(ax = ax_K_star, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=24, y_fontsize=24)
# ax_K_star.legend(fontsize = 14, title_fontsize = 16, loc = 0, title = r'$p$ ; $k_{\mathrm{step}}[\mathrm{ min}^{-1}]$ ; $t_{\mathrm{{\scriptsize act}}}[\mathrm{ days}]$')
#ax_K_star.set_ylabel(r'$D({\Omega_{{\textrm{act}}}|\Omega_{0}})$')
#ax_K_star.set_xlim(left = -21, right = -7)
#ax_K_star.set_ylim(bottom = 2e-7, top = 3e-2)
#ax_K_star.set_yticks(range(3, 16, 2))
#ax_K_star.set_yticklabels([1, 0.1, 0.01])
fig_K_star.savefig('../../../Figures/repertoire_entropy/theory/K_star_'+energy_model+'.pdf')



