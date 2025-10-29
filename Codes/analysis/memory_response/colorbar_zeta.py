import sys
sys.path.append('../../my_lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Immune_System/primary_response/in/'
#For simulation in C++
# ./EF_response_v2.x -a 6 -b 0.5 -k 1 -t 0 -T 6.5 -E MJ -C 10000 -B 100000000 -s TACNSEYPNTTKCGRWYC -q 2.0
# ./EF_response_v3.x -a 6 -b 0.5 -k 1 -t 0 -T 7.5 -E TCRen -C 10000 -B 100000000 -s EYTACNSEYPNTTKCGRWYCGRYPN -q 2.0 --ensemble -N 10
#--------------- PARAMETERS ---------------------
N_ens = 1
L_0 = 1e9
T0 = 0
Tf = 9
Tf_sim = 7
#Tf = 10
dT = 0.05
lambda_A = 6
k_step = 1/(60*2) #s^-1
k_step = k_step*3600 # hour^-1
k_step = k_step*24 #days^-1
lambda_B = 3 * np.log(2) #(days)^-1
k_on = 1e6*24*3600; #(M*days)^-1
N_c = 1e5*10
#N_c = 1e5
#E_ms = -27.63
E_ms = -24
C = 1e4
AA = 1

ps = [1, 3]
t_ten = [2.76, 4.37]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
#color_list = np.array([(228,75,41), (125,165,38), (76,109,166), (215,139,45)])
color_list = np.array([my_red, my_blue2, my_green, my_gold, my_brown])

#colors_p = np.flip(['tab:blue', 'tab:red', 'tab:blue'])
#colors_p = np.flip(['tab:blue','tab:green','tab:red'])
colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

#colors_R = [['tab:grey', 'tab:grey', 'tab:blue', 'tab:blue'], ['tab:grey', 'tab:grey', 'tab:green', 'tab:green'], ['tab:grey', 'tab:grey', 'tab:red', 'tab:red'], ['tab:red', 'tab:red', 'tab:red', 'tab:red']]
colors_R = []
for i in range(len(ps)):
    colors_R.append(['tab:grey', colors_p[i], colors_p[i], colors_p[i]])

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
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
beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
print('beta_r = %.1f'%beta_r)

#--------------------------Proofreading properties--------------------------
beta_pr, E_pr, Kd_pr = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_pr = %.2f'%beta_pr)

t_prime = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*N_c))
print('--------')

min_E = -17.8
min_E = -19.5
max_E = -10

min_E = np.log(1e-9)
max_E = np.log(1e-4)

# color_vals = np.linspace(0, 1, 9)
# cmap = plt.get_cmap('Reds')
# my_colors = [cmap(val) for val in color_vals] 

ps = np.linspace(0.3, 1.4, 5)

project = 'memory_response'
subproject = 'data'

output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
os.makedirs(output_plot, exist_ok=True)

fig, ax = plt.subplots(figsize = (5*1.62, 1), linewidth = 6, gridspec_kw={'left':0.1, 'right':.9, 'bottom':.1, 'top': .3})
col_map = 'rainbow'
#mpl.colorbar.ColorbarBase(ax, cmap=col_map, orientation = 'vertical')
# and create another colorbar with:
colorbar = mpl.colorbar.ColorbarBase(ax, cmap=plt.get_cmap(col_map + '_r'), orientation = 'horizontal')
# colorbar.set_ticks([np.linspace(0, 1, 9)])
#ax.set_xticks(np.linspace(0, 1, 5))
#ax.set_xticklabels([r'$%.0f \cdot 10^{%d}$'%(10**(np.log10(np.exp(min_E + i*(max_E - min_E)/4))%1), int(np.log10(np.exp(min_E + i*(max_E - min_E)/4)))) for i in np.arange(0, 5, 1)], fontsize = 30)
# ax.set_xticklabels([r'$10^{%d}$'%(int(np.log10(np.exp(min_E + i*(max_E - min_E)/5)))) for i in np.arange(0, 6, 1)], fontsize = 30)
ax.set_xticks([i for i in np.linspace(0, 1, 5)], [r'$%.1f$'%i for i in ps], fontsize = 26)
ax.xaxis.tick_top()
fig.savefig(output_plot + "/colorbar_zeta.pdf")


