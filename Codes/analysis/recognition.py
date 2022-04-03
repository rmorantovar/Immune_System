import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
from Immuno_models import*
import scipy.special as sc
import pickle
from tqdm import tqdm

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
fig_p_a, ax_p_a = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15, 'bottom': .15, 'right':.8})
fig_Q_P0_R, ax_Q_P0_R = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15, 'bottom': .15, 'right':.8})
ax_densities_R = ax_Q_P0_R.inset_axes([0.7, 0.7, 0.3, 0.3])
fig_M_r, ax_M_r = plt.subplots(figsize=(12,8), gridspec_kw={'left':0.15, 'bottom': .15, 'right':.95})


file_name = 'PWM_Adams_etal_2016_1.pkl'
L = 12
gamma = 1/(3600*1) #s^-1
k_on = 1e6 #(M*s)^-1
alpha = 1 # days^-1
alpha = alpha/(60*60*24) # s^-1
e_MS = -27
N_r = 1e5
N_c = 1e0

es = np.linspace(e_MS, e_MS+50, 100)
de = es[1]-es[0]


#file_PWM = open(Text_files_path + 'Data/'+ file_name, 'rb')
#log10Kd_WT = -8.91169118
#PWM_data = (pickle.load(file_PWM) - log10Kd_WT)
#PWM_data = (pickle.load(file_PWM) - log10Kd_WT)*np.log(10)
#file_PWM.close()

antigen_aa = 'TACNSEYPNTTK'
antigen_list = [i for i in antigen_aa]
antigen = np.array([], dtype = int)
for i, aa in enumerate(antigen_list):
    index = Alphabet_list.index(aa)
    antigen = np.append(antigen, int(index))
PWM_data = M2[:,antigen]

#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

min_E_data = np.sum([np.min(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + e_MS
avg_E_data = np.sum([np.mean(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + e_MS
var_E_data = np.sum([np.var(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))])
max_E_data = np.sum([np.max(PWM_data[:,i]) for i in np.arange(len(PWM_data[0,:]))]) + e_MS

#---------------------------------------------------------------------- #Calculate Free energy from Z(PWM)
Tmin = .01
Tmax = 50
Ts = np.linspace(Tmin, Tmax, 20000)
lambdas = 1/Ts[:-1]
F_PWM = -Ts*np.log(Z_PWM(PWM_data, Ts))
Us = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts)) + e_MS
print(Us[0])
dU = np.diff(Us)
#---------------------------------------------------------------------- #Calculate entropy functions and distributions
Omega = 20**L
Z = Z_PWM(PWM_data, 1)
S = np.cumsum(lambdas[:-1]*dU)
S_gaussian = -(Us[:-1]-avg_E_data)**2/(2*var_E_data) + L*np.log(20) - 0.5*np.log(2*np.pi*var_E_data)

P0 = 1/Omega
Q0 = np.exp(S)*P0
Q0_gaussian = np.exp(S_gaussian)*P0

P_eq = np.exp(-(Us[:-1]-e_MS))/Z
Q_eq = np.exp(S)*P_eq
Q_guassian = np.exp(S_gaussian)*P_eq

KL_D = -np.sum(Q_eq*np.log(Q_eq)*dU)


K_d_array = np.exp(Us[:-1])

lambda_r = lambdas[np.where(k_on*K_d_array>gamma)][0]

#--------- Plot p_a --------------
for q in np.arange(1, 4):
    ax_p_a.plot(K_d_array, 1/(1+(K_d_array*k_on/gamma)**q), label = '%d'%q)
ax_p_a.vlines([K_d_array[np.where(lambdas[:-1]<1)][0], K_d_array[np.where(lambdas[:-1]<lambda_r)][0]], 1e-6, 2e0, color = ['grey', 'black'], linestyle = ':', linewidth = 2)
ax_p_a.hlines([1, .5], 1e-12, 1e-5, color= 'grey', linestyle = '--')
my_plot_layout(ax=ax_p_a, xlabel = r'$K_{d} [M]$', ylabel = r'$p_a$', xscale = 'log', yscale = 'log', x_fontsize = 34, y_fontsize = 34)
ax_p_a.set_title('$\gamma=%.1e$ $\mathrm{min}^{-1}$'%(gamma*60), fontsize = 25)
leg = ax_p_a.legend(title = '$q$', title_fontsize = 25, fontsize = 20, loc = 3)
leg._legend_box.align = "center"
ax_p_a.set_xlim(np.exp(-25), 1e-5)
ax_p_a.set_ylim(1e-5, 1.2)
fig_p_a.savefig('../../Figures/7_Recognition/activation_single_cell.pdf')
#---------------------------------

#--------- Plot Q_R --------------
colors_data = np.array(['indigo', 'indianred'])
labels_data = np.array(['H1', 'H3', 'MHC'])
p_a = (1/(1 + (k_on*np.exp(Us[:-1])/gamma)**2 ) )
ax_densities_R.plot(K_d_array, Q0*N_r, linewidth = 2, alpha = .4, linestyle = ':', marker = '', ms = 6, color = colors_data[0], label = labels_data[0])
ax_Q_P0_R.plot(K_d_array, Q0*N_r, linewidth = 2, alpha = .4, linestyle = ':', marker = '', ms = 6, color = colors_data[0], label = labels_data[0])
for r, rho_A in enumerate((1e3*np.logspace(8, np.log10(np.exp(35)), 5)/N_A)):
	p_e = rho_A*k_on/alpha #probability of engagement up to time t
	#Q_R = Q0*p_e*p_a*N_r
	Q_R = N_r*Q0*(1-np.exp(-p_e*p_a*N_c))
	#ax_Q_P0_R.plot(K_d_array, p_e*p_a, linewidth = 2, alpha = .1*r, linestyle = '--', marker = '', ms = 6, color = colors_data[0])
	ax_Q_P0_R.plot(K_d_array, 1-np.exp(-p_e*p_a*N_c), linewidth = 2, alpha = .2*r, linestyle = '--', marker = '', ms = 6, color = colors_data[0])
	ax_Q_P0_R.plot(K_d_array, Q_R, linewidth = 2, alpha = .2*r, linestyle = '-', marker = '', ms = 6, color = colors_data[0])
	ax_Q_P0_R.vlines(K_d_array[np.where((Q_R)==np.max(Q_R))[0]], ax_Q_P0_R.get_ylim()[0], np.max(Q_R), color = colors_data[0], linestyle = '-', alpha = .8)
ax_Q_P0_R.vlines([K_d_array[np.where(lambdas[:-1]<1)][0], K_d_array[np.where(lambdas[:-1]<lambda_r)][0]], ax_Q_P0_R.get_ylim()[0], np.max(Q_R), colors = ['grey', 'black'], linestyle = ':', linewidth = 2)

#ax_Q_P0.plot(Us[:-1][::100], np.log(R)[::100], linewidth = 2, alpha = .6, linestyle = '-', marker = '', ms = 6, color = 'darkred')
ax_Q_P0_R.hlines([1, N_r], np.exp(min_E_data), 1e-2, color = 'grey', linestyle = '--')
my_plot_layout(ax=ax_Q_P0_R, xscale = 'log', yscale = 'log', ylabel = r'${Q(K_d)}$', xlabel = r'${K_d}$ $[M]$')
#ax_Q_P0_R.legend(loc = 4, fontsize = 20, title = r'Sequence', title_fontsize = 16)
ax_Q_P0_R.set_ylim(bottom = N_r**(-1/2), top = 2*np.max(Q_R))
ax_Q_P0_R.set_xlim(left = np.exp(min_E_data-1), right = 1e-3)
#ax_Q_P0_R.set_xlim(right = 0, left = -10.5);
ax_densities_R.hlines([1, N_r], np.exp(min_E_data), np.exp(avg_E_data), color = 'grey', linestyle = '--' )
my_plot_layout(ax=ax_densities_R, xscale = 'log', yscale = 'log', ylabel = r'$N_r\cdot Q_0$', xlabel =r'$\log_{10}{K_d} [M]$', ticks_labelsize = 10, x_fontsize=14, y_fontsize=16)
#ax_densities_R.set_xlim(right=0)
#ax_densities_R.xaxis.tick_top()
#ax_densities_R.xaxis.set_label_position('top') 
fig_Q_P0_R.savefig('../../Figures/7_Recognition/Q_P0_R.pdf')
#---------------------------------

#--------- Plot M_r --------------
rho_A_array = np.logspace(9, 15, 40)
x, y = np.meshgrid(K_d_array, rho_A_array)
#z = np.log10(Q0*(y*k_on/alpha)*((gamma/(gamma + (k_on*x)**3 ) ))*N_r)
p_e2 = (1e3*y/N_A)*k_on/alpha
p_a2 = 1/(1 + (k_on*x/gamma)**2)
z = np.log10(N_r*Q0*(1-np.exp(-p_e2*p_a2*N_c)))


cs = ax_M_r.contourf(x, y, z, cmap = plt.cm.RdBu_r, levels = np.linspace(-6,4,500), vmax =6 , vmin = -6)
cs2 = ax_M_r.contour(cs, levels=[0], colors='k', linestyles = 'dashed', linewidths = 4)
my_plot_layout(ax=ax_M_r, xlabel=r'$K_{d}$ $[M]$', ylabel=r'$\rho_A$ $[\mathrm{copies}/mL]$', yscale='log', xscale='log', x_fontsize = 34, y_fontsize = 34)


print(lambda_r)
ax_M_r.plot(K_d_array, 5e-3*K_d_array**(2-.78), linewidth = 2, color = 'dimgray', linestyle = 'dashed')
ax_M_r.plot(K_d_array, 1.5e-24*K_d_array**-lambda_r, linewidth = 2, color = 'dimgray', linestyle = 'dashed')
ax_M_r.text(2e-8, 1e-13, r'$\sim K_d^{q-\lambda(\epsilon_\gamma)}$', fontsize = 26 )
#ax_M_r.text(2e-11, 3e-16, r'$\sim K_d^{-\lambda(\epsilon_r)}$', fontsize = 26 )

cbar = fig_M_r.colorbar(cs, ticks=np.linspace(-6,4,5))
cbar.ax.set_ylabel(r'$\log_{10}{M_r}$', fontsize = 25)
cbar.set_ticklabels(np.linspace(-6,4,5))
cbar.ax.tick_params(labelsize = 25)
cbar.add_lines(cs2)
cbar.lines[-1].set_linestyles(cs2.linestyles)
cbar.lines[-1].set_linewidths(2.5)
ax_M_r.set_xlim(np.exp(e_MS+1), 1e-5)
ax_M_r.set_ylim(np.min(rho_A_array), np.max(rho_A_array))

fig_M_r.savefig('../../Figures/7_Recognition/recognition.pdf')
#---------------------------------


