import sys
sys.path.append('../../lib/')
from functions_1 import*
from classes import*
# from functions_2 import*


Epitopes = [['TACNSYPNTAKCRWYR', 'TACNSYPNTAKCRWYR'],
			['TACNSYPNTAKCRWYR', 'TACNIYNNTAKCRWYR'], 
			['TACNSYPNTAKCRWYR', 'TACNEYPNTAKCHDVR'],
			['TACNSYPNTAKCRWYR', 'QRCNSYPNTAICMTYR'],
			['TACNSYPNTAKCRWYR', 'TYDNSNPTTFKERWYM']] #6
antigens = [0, 1, 2, 3, 4]
# antigens = ['TACNSEYPNTTRAKCGRWYR', 'TVCNSRYPNTTRLKFGRWYR', 'TVCPSRYENTIRLKFGRWYR'] #4
# antigens = ['TACNSEYPNTTRAKCGRWYR', 'TACRSEYPNTTRAKCGRKYR', 'TACRSEYPEFTRAKCGRKYR'] #2
# antigens = ['TACNSEYPNTTRAKCGRWYR', 'TACNSEYPNTTRAKCGRWYR', 'TACNSEYPNTTRAKCGRWYR'] #0
energy_model = 'TCRen'
L0 = 1e6
N_ens = 1
T0 = 0
Tf = 10
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
b0 = 1e5*10
C = 1e4
AA = 1
p = 3

colors_inf = ['gray', my_green, my_red]

for antigen in antigens:
	fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	epitopes = Epitopes[antigen]
	Q0_total = np.zeros(400000-1)
	for e, epitope in enumerate(epitopes):
		print(epitope)
		l = len(epitope)
		E_m = -3
		epitope_seq = from_aa_to_i(epitope, energy_model, '../../')
		# Calculate motif
		motif = get_motif(epitope_seq, energy_model, '../../')
		# Normalize motif
		for i in range(l):
			E_m+=np.min(motif[:, i], axis=0)
			motif[:, i] -= np.min(motif[:, i], axis=0)

		# Calculate Q0, Es, dE, betas
		Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
		Q0_total += Q0
		
		ax.plot(np.exp(Es), Q0*L0, color = colors_inf[e+1], ls = ':')

		for t in [6.0, 6.5]:
			u_on, p_a, R, QR = calculate_QR(Q0[:-1], k_on, k_step, np.exp(lambda_A*t)/N_A, Es, p, lambda_A, b0, dE)
			ax.plot(np.exp(Es[:-1]), QR*L0, color = colors_inf[e+1])

		u_on, p_a, R, QR = calculate_QR(Q0[:-1], k_on, k_step, np.exp(lambda_A*7)/N_A, Es, p, lambda_A, b0, dE)
		ax.plot(np.exp(Es[:-1]), QR*L0, color = colors_inf[e+1], label = r'$%d$'%(e+1))
		
		

		ax.hlines([1, 1/(20**16)*L0], np.exp(-21), np.exp(-4), color = 'k', ls = '--', alpha = .6)

	my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
	ax.legend(fontsize = 18, title_fontsize = 22, loc = 0, title = r'$\mathrm{Epitope}$')
	# ax.set_xlim(left = np.exp(-21), right = np.exp(-13+8.5))
	ax.set_ylim(bottom = 1e-10)#, top = 5e5)
	fig.savefig('../../../Figures/memory/epitopes/Q0_'+energy_model+'_'+str(antigen)+'.pdf')


