import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

E_ms = -25
antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
L=len(antigen)
#----------------------------------------------------------------
model = 'TCRen'
#--------------------------Energy Motif--------------------------
PWM_data, M, Alphabet = get_motif(antigen, model, Text_files_path)
print('min_e_PWM=%.4f'%(np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
print('mean_e_PWM=%.4f'%(np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])))
#Change values by the minimum
for i in np.arange(L):
    PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)

#--------------------------Entropy function--------------------------
Es, dE, rho_01, rho_02, betas = calculate_Q0_2(0.05, 50, 400000, PWM_data, E_ms, L)

fig, ax = plt.subplots(figsize=(4.5*2,3.0*2), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.1, 'top': 0.96})

ax.plot(Es[:-1], rho_01)
ax.plot(Es, rho_02, '--')

my_plot_layout(ax = ax, yscale = 'log', xscale = 'linear')

fig.savefig('../../Figures/0_Shape_Space/entropy.pdf')