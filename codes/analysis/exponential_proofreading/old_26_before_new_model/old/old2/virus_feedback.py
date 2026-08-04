import sys
sys.path.append('../../lib/')
from functions_memory import*
import random

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/primary_response/'

lamA = 6
lamZ = 3
alpha = 2e-10
delta = 1
L0 = 1e8
time = np.linspace(0, 7, 1000)
dt = time[1] - time[0]
colors = [my_blue, my_gold, my_green]

antigen = 'TACNSEYPNTTRAKCGRWYC' #l=20
#antigen = 'CMFILVWYAGTSQNEDHRKPFMRTPTRMCW' #l=30
l=len(antigen)
print('--------')
#----------------------------------------------------------------
model = 'TCRen'
#model = 'MJ2'
#--------------------------Energy Motif--------------------------
antigen_seq = from_aa_to_i(antigen, model, '../../')
motif = get_motif(antigen_seq, model, Text_files_path)
#Change values by the minimum
E_m = -3
for i in np.arange(l):
	E_m+=np.min(motif[:,i], axis=0)
	motif[:,i]-=np.min(motif[:,i], axis=0)
print('min_e_PWM=%.4f'%(np.sum([np.min(motif[:,i]) for i in range(len(motif[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(motif[:,i]) for i in range(len(motif[0,:]))])))
#Change values by the minimum
for i in np.arange(l):
    motif[:,i]-=np.min(motif[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_m, l)
Kds = np.exp(Es[:-1])

fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
ax.plot(Kds, P_min_e_Q0(L0, Q0, dE))
ax.plot(Kds, np.cumsum(P_min_e_Q0(L0, Q0, dE)*dE))

my_plot_layout(ax=ax, yscale = 'log', xscale = 'log', ticks_labelsize = 30)
#ax_Q_.set_xticks([])
ax.set_xlim(right = 1.2) #use 1e-3 for other plots
ax.set_ylim(bottom = 1e-15, top = 1.5)
fig.savefig('../../../Figures/memory_response/Virus_feedback/Q0_L0-%.0e_'%(L0)+model+'.pdf')
plt.close(fig)

fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

v = np.ones_like(time)
v_maxs = []
Z0s = []
for i_z in range(50):
	r1 = random.uniform(0, 1)
	r2 = Es[:-1][np.cumsum(P_min_e_Q0(L0, Q0, dE)*dE)<r1][-1]
	Z0_i = np.exp(r2)**-1
	Z0s.append(Z0_i)
	Z = Z0_i + 1e4*Z0_i/(1+Z0_i*np.exp(-lamZ*(time)))
	Z = Z0_i*np.exp(lamZ*(time - 2))
	lamA_mod = (lamA*(1-1/(1+(alpha*Z)**-1)) - delta*(1/(1+(alpha*Z)**-1)))
	for t in range(1, len(time)):
		lamA_mod_t = (lamA*(1-1/(1+(alpha*Z[t])**-1)) - delta)#*(1/(1+(alpha*Z[t])**-1)))
		v[t] = v[t-1] + v[t-1]*lamA_mod_t*dt
	v_max = np.max(v)
	v_maxs.append(v_max)
	# v = np.exp((lamA*(1-1/(1+(alpha*Z)**-1)) - delta)*time)
	# ax.plot(time, Z, color = my_blue, ls = '--', lw = 1, alpha = .2)
	ax.plot(time, v, color = my_red, ls = '-', lw = 1, alpha = .2)
	# ax.plot(time, lamA_mod, lw = 1, alpha = 0.2, color = my_purple)

Z0 = np.exp(np.sum(Es[:-1]*P_min_e_Q0(L0, Q0, dE)*dE))**-1
Z = Z0 + 1e4*Z0/(1+Z0*np.exp(-lamZ*(time))/1)
Z = Z0*np.exp(lamZ*(time - 2))
lamA_mod = (lamA*(1-1/(1+(alpha*Z)**-1)) - delta*(1/(1+(alpha*Z)**-1)))
for t in range(1, len(time)):
	lamA_mod_t = (lamA*(1-1/(1+(alpha*Z[t])**-1)) - delta)#*(1/(1+(alpha*Z[t])**-1)))
	v[t] = v[t-1] + v[t-1]*lamA_mod_t*dt
# v = np.exp((lamA*(1-1/(1+(alpha*Z)**-1)) - delta)*time)
# ax.plot(time, Z, label = r'$\cal Z$', color = my_blue, ls = '--', lw = 3)
ax.plot(time, v, label = r'$10^{%d} $ '%(int(np.log10(1/Z0))), color = my_red, ls = '-', lw = 3)
# ax.plot(time, lamA_mod, lw = 3, color = my_purple)

my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
ax.legend(fontsize = 22, title_fontsize = 24, loc = 0, title = r'$K^*\, [{\rm M}]$')
ax.set_xlim(left = 0, right = 7)
ax.set_ylim(bottom = 100)#, top = 1e15)
# ax.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
# ax.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
fig.savefig('../../../Figures/memory_response/Virus_feedback/virus_feedback.pdf')


fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

# lamAs = [4, 6, 8]
# for i_A, lamA in enumerate(lamAs):
# 	Z0s = np.logspace(6, 8, 20)
# 	v_maxs = []
# 	for i_z, Z0 in enumerate(Z0s):
# 		Z = Z0*np.exp(lamZ*(time-2))
# 		v = np.exp(lamA*(1-alpha*Z)*time)
# 		v_max = np.max(v)
# 		v_maxs.append(v_max)

ax.plot(Z0s, v_maxs, label = r'$%d$ '%(lamA), color = my_gold, ls = '', lw = 3, marker = '^', ms = 10)
# ax.plot(Z0s, Z0s**(- lamA/lamZ), label = r'$%d$ '%(lamA), color = my_gold, ls = '', lw = 3, marker = '^', ms = 10)

my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
ax.legend(fontsize = 22, title_fontsize = 24, loc = 0, title = r'$\lambda_A \,[{\rm days}^{-1}]$')
ax.set_xlim(left = 1e6, right = 1e8)
ax.set_ylim(bottom = np.min(v_maxs)*0.5, top = np.max(v_maxs)*2)
# ax.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
# ax.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
fig.savefig('../../../Figures/memory_response/Virus_feedback/max_viral_load_Z0.pdf')


