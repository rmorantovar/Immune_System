import sys
sys.path.append('../library/')
from Immuno_models import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

RT = 0.593 #kcal/mol
N_A = 6.02214076e23
#----------------------------------------------------------------------
#Matrix = 'BLOSUM62'
Matrix = 'MJ2'
#Matrix = 'MM'

if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
if(Matrix == 'MM'):
    M2 = (np.loadtxt(Text_files_path + Matrix + '.txt')+1)*e0
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
if(Matrix == 'BLOSUM62'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,25))
    Alphabet = np.array(['A', 'R'  ,'N' , 'D'  ,'C' , 'Q'  ,'E'  ,'G'  ,'H' , 'I'  ,'L'  ,'K'  ,'M' , 'F' , 'P' , 'S'  ,'T' , 'W' , 'Y' , 'V' , 'B' , 'Z'  ,'X',  '*'])
L_alphabet = len(Alphabet)
#----------------------------------------------------------------------
fig_p_a, ax_p_a = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.18, 'right':.95, 'bottom':.15})

file_name = 'PWM_Adams_etal_2016_1.pkl'

N_r = 2e5
N_r = 1e6
T0 = 2
Tf = 5.8
Tf = 6
dT = 0.1
days = np.arange(0, Tf, 1)
time = np.linspace(T0, Tf, int((Tf-T0)/dT))
lambda_A = 6 #days^-1
k_pr = .1 # hour^-1
k_pr = k_pr*24 #days^-1
qs = [1, 2]
beta = 1*lambda_A
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e3
E_ms = -27

es = np.linspace(E_ms, E_ms+30, 100)
de = es[1]-es[0]

antigen_aa  = 'TACNSEYPNTTK'
#antigen_aa = 'TANSEYPNTK'
#antigen_aa = 'MRTAYRNG'
#antigen_aa = 'MRTAY'

L=len(antigen_aa)

antigen_list = [i for i in antigen_aa]
antigen = np.array([], dtype = int)
for i, aa in enumerate(antigen_list):
    index = Alphabet_list.index(aa)
    antigen = np.append(antigen, int(index))
PWM_data = M2[:,antigen]

#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

min_E_data = np.sum([np.min(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + E_ms
avg_E_data = np.sum([np.mean(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + E_ms
var_E_data = np.sum([np.var(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))])
max_E_data = np.sum([np.max(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + E_ms

Es, dE, Q0, lambdas = calculate_Q0(0.01, 50, PWM_data, E_ms, L)

Kd = np.exp(Es[:-1])
beta_r = lambdas[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]
E_r = Es[:-1][np.cumsum(Q0*dE)<(1/N_r)][-1]

beta_act = lambdas[np.where(k_on*Kd>k_pr)][0]

#--------- Plot p_a --------------
for q in np.arange(1, 4):
    ax_p_a.plot(Kd, 1/(1+(Kd*k_on/k_pr)**q), label = '%d'%q)
#ax_p_a.vlines([Kd[np.where(lambdas[:-1]<1)][0], Kd[np.where(lambdas[:-1]<beta_act)][0]], 1e-6, 2e0, color = ['grey', 'black'], linestyle = ':', linewidth = 2)
ax_p_a.hlines([1], np.min(Kd), np.max(Kd), color= 'grey', linestyle = '--')
my_plot_layout(ax=ax_p_a, xlabel = r'$K_{d} [M]$', ylabel = r'$p_a$', xscale = 'log', yscale = 'log', x_fontsize = 34, y_fontsize = 34)
ax_p_a.set_title('$k_pr=%.1f$ $\mathrm{hours}^{-1}$'%(k_pr/24), fontsize = 25)
leg = ax_p_a.legend(title = '$q$', title_fontsize = 25, fontsize = 20, loc = 3)
leg._legend_box.align = "center"
ax_p_a.set_xlim(left = np.min(Kd), right = np.max(Kd))
#ax_p_a.set_ylim(1e-5, 1.2)
fig_p_a.savefig('../../Figures/7_Recognition/activation_single_cell.pdf')
#---------------------------------

#--------- Plot QR --------------
fig_Q0_QR, ax_Q0_QR = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_Q0_QR_n, ax_Q0_QR_n = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
fig_Q0_QR_n2, ax_Q0_QR_n2 = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
for q in [1, 2]:
	#ax_densities_R = ax_Q0_QR.inset_axes([0, 0.7, 0.35, 0.35])
	colors_data = np.array(['indigo', 'indianred'])
	labels_data = np.array(['H1', 'H3', 'MHC'])
	#p_a = (1/(1 + (k_on*np.exp(Es[:-1])/k_pr)**q) )
	#ax_densities_R.plot(Kd, Q0*N_r, linewidth = 2, alpha = .4, linestyle = ':', marker = '', ms = 6, color = colors_data[0], label = labels_data[0])
	ax_Q0_QR.plot(Kd, Q0*N_r, linewidth = 5, alpha = .4, linestyle = '-', marker = '', ms = 6, color = 'grey', label = labels_data[0])
	ax_Q0_QR_n.plot(Kd, Q0, linewidth = 5, alpha = .4, linestyle = '-', marker = '', ms = 6, color = 'grey', label = labels_data[0])
	#ax_Q0_QR.plot(Kd, Q0*N_r, linewidth = 2, alpha = .4, linestyle = ':', marker = '', ms = 6, color = colors_data[0], label = labels_data[0])
	if q == 1:
		for r, rho_A in enumerate((np.logspace(np.log10(np.exp(lambda_A*T0*.9)), np.log10(np.exp(lambda_A*Tf)), 3)/N_A)):
			u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, rho_A, Es, q, lambda_A, N_c, dE)
			#ax_Q0_QR.plot(Kd, R, linewidth = 2, alpha = .2*r, linestyle = '--', marker = '', ms = 6, color = colors_data[0])
			ax_Q0_QR.plot(Kd, N_r*QR, linewidth = 4, alpha = .3*(r+1), linestyle = '-', marker = '', ms = 6, color = 'indianred')
			if(r==2):
				ax_Q0_QR_n.plot(Kd, QR/np.sum(QR*dE), linewidth = 5, alpha = .3*(r+1), linestyle = '-', marker = '', ms = 6, color = 'indianred')
				D_KL = np.sum(dE*(QR/np.sum(QR*dE))*np.log2((QR/np.sum(QR*dE))/Q0))
				ax_Q0_QR_n2.plot(Kd, (QR/np.sum(QR*dE))/Q0, linewidth = 3, alpha = .3*(r+1), linestyle = '-', marker = '', ms = 6, color = 'indianred', label = '%.2f'%(D_KL))
			#ax_Q0_QR.vlines(Kd[np.where((QR)==np.max(QR))[0]], ax_Q0_QR.get_ylim()[0], 1, color = colors_data[0], linestyle = '-', alpha = .8)
	if q == 2:
		for r, rho_A in enumerate((np.logspace(np.log10(np.exp(lambda_A*T0*1.37)), np.log10(np.exp(lambda_A*Tf)), 3)/N_A)):
			u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_pr, rho_A, Es, q, lambda_A, N_c, dE)
			#ax_Q0_QR.plot(Kd, R, linewidth = 2, alpha = .2*r, linestyle = '--', marker = '', ms = 6, color = colors_data[0])
			ax_Q0_QR.plot(Kd, N_r*QR, linewidth = 5, alpha = .3*(r+1), linestyle = '-', marker = '', ms = 6, color = 'dodgerblue')
			if(r==2):
				ax_Q0_QR_n.plot(Kd, QR/np.sum(QR*dE), linewidth = 5, alpha = .3*(r+1), linestyle = '-', marker = '', ms = 6, color = 'dodgerblue')
				D_KL = np.sum(dE*(QR/np.sum(QR*dE))*np.log2((QR/np.sum(QR*dE))/Q0))
				ax_Q0_QR_n2.plot(Kd, (QR/np.sum(QR*dE))/Q0, linewidth = 3, alpha = .3*(r+1), linestyle = '-', marker = '', ms = 6, color = 'dodgerblue', label = '%.2f'%(D_KL))
			#ax_Q0_QR.vlines(Kd[np.where((QR)==np.max(QR))[0]], ax_Q0_QR.get_ylim()[0], 1, color = colors_data[0], linestyle = '-', alpha = .8)
	#ax_Q0_QR_n2.vlines()		
	my_plot_layout(ax=ax_Q0_QR, xscale = 'log', yscale = 'log', ticks_labelsize = 30)
	my_plot_layout(ax=ax_Q0_QR_n, xscale = 'log', yscale = 'linear', ticks_labelsize = 30)
	my_plot_layout(ax=ax_Q0_QR_n2, xscale = 'log', yscale = 'log', ticks_labelsize = 30)
	ax_Q0_QR.set_ylim(bottom = N_r**(-1/2), top =N_r)
	ax_Q0_QR_n.set_ylim(top =.5)
	#ax_Q0_QR_n2.set_ylim()
	ax_Q0_QR_n2.legend(fontsize = 30, title = r'$D_{KL}$', title_fontsize = 35, loc = 3)

fig_Q0_QR.savefig('../../Figures/7_Recognition/Q0_QR.pdf')
fig_Q0_QR_n.savefig('../../Figures/7_Recognition/Q0_QR_n.pdf')
fig_Q0_QR_n2.savefig('../../Figures/7_Recognition/Q0_QR_n2.pdf')
#---------------------------------

#--------- Plot M_r --------------
#K_d_array2 = Kd[[int(i) for i in np.logspace(0, np.log10(len(Kd)), 20)]]
K_d_array2 = Kd[::50]
for q in np.arange(1, 4):
	fig_M_r, ax_M_r = plt.subplots(figsize=(10,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})
	rho_A_array = np.logspace(10, 15, 40)
	rho_A_array = np.logspace(np.log10(np.exp(lambda_A*T0)), np.log10(np.exp(lambda_A*Tf)), 40)
	x, y = np.meshgrid(K_d_array2, rho_A_array)
	#z = np.log10(Q0*(y*k_on/lambda_A)*((k_pr/(k_pr + (k_on*x)**3 ) ))*N_r)
	p_e2 = (y/N_A)*k_on/lambda_A
	p_a2 = 1/(1 + (k_on*x/k_pr)**q)
	#z = np.log10(N_r*Q0[[int(i) for i in np.logspace(0, np.log10(len(Kd)), 20)]]*(1-np.exp(-p_e2*p_a2*N_c)))
	z = np.log10(N_r*Q0[::50]*(1-np.exp(-p_e2*p_a2*N_c)))

	cs = ax_M_r.contourf(x, y/1e3, z, cmap = plt.cm.RdBu_r, levels = np.linspace(-8.5,5,400), vmax =8.5 , vmin = -8.5)
	cs2 = ax_M_r.contour(cs, levels=[0], colors='k', linestyles = 'dashed', linewidths = 4)
	my_plot_layout(ax=ax_M_r, yscale='log', xscale='log', x_fontsize = 35, y_fontsize = 35, ticks_labelsize= 30)

	ax_M_r.vlines([Kd[np.cumsum(N_r*Q0*dE)<1][-1], Kd[lambdas[:-1]<q][0], Kd[np.where(lambdas[:-1]<beta_act)][0]], ax_M_r.get_ylim()[0], ax_M_r.get_ylim()[1], colors = ['black','grey','darkred'], linestyles = ['--', ':', ':'], linewidths = [2, 2, 2])

	cbar = fig_M_r.colorbar(cs, ticks=np.linspace(-8,4,7))
	#cbar.ax.set_ylabel(r'$\log_{10}{Q}$', fontsize = 25)
	cbar.set_ticklabels(np.linspace(-8,4,7))
	cbar.ax.tick_params(labelsize = 30)
	cbar.add_lines(cs2)
	cbar.lines[-1].set_linestyles(cs2.linestyles)
	cbar.lines[-1].set_linewidths(2.5)
	#ax_M_r.set_xlim(np.exp(E_ms+1), 1e-5)
	ax_M_r.set_ylim(np.min(rho_A_array/1e3)*100, np.max(rho_A_array/1e3))
	ax_M_r.set_xlim(K_d_array2[5], 2e-5)

	fig_M_r.savefig('../../Figures/7_Recognition/recognition_q-%d.pdf'%(q))
#---------------------------------


