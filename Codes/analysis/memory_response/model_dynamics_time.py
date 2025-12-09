import sys
sys.path.append('../../my_lib/')
# from functions_2 import*
from funcs import*
plt.rcParams['text.usetex'] = True

Text_files_path = '/Users/robertomorantovar/Library/CloudStorage/Dropbox/Research/Immune_system/'

#--------------- PARAMETERS ---------------------
N_ens = 1
L_0s = [1e9]
T0 = 0
Tf = 8
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
#b0 = 1e5
#E_ms = -27.63
E_ms = -24
C = 1e4

ps = [1, 2, 2.5, 3, 4]
ps = [1, 3]

transparency_n = [1]

color_list = np.array([my_blue, my_gold, my_green, my_red, my_purple2, my_brown, my_blue2, my_yellow, my_purple, my_green2])#
color_list = np.array([my_blue2, my_red, my_green])
color_list = np.array([my_red, my_blue2, my_brown, my_red, my_gold])

colors_p = []
for i in range(len(color_list)):
        colors_p.append(np.array(color_list[i]))

colors_R = []
for i in range(len(ps)):
    colors_R.append([colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i], colors_p[i]])

time_array = np.linspace(T0, Tf, int((Tf-T0)/dT))
energy_models = ['TRCen']
energy_model = 'TCRen'
models_name = ['exponential']#, 'linear',]
growth_models = [0]
linear = 0
alpha_p = [.6, 1]

project = 'memory_response'
subproject = 'multi-epitope'
subproject = 'model_dynamics_time'

output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
os.makedirs(output_plot, exist_ok=True)

antigen = 'TACNSEYPNTTRAKCGRWYC' #L=20
L=len(antigen)
print('--------')
print('L=%d'%(L))
#----------------------------------------------------------------
model = 'TCRen'
#--------------------------Energy Motif--------------------------
motif = get_motif(from_aa_to_i(antigen, energy_model, Text_files_path), model, Text_files_path)
#Change values by the minimum
for i in np.arange(L):
    motif[:,i]-=np.min(motif[:,i], axis=0)

avg_E = np.sum([np.mean(motif[:,i]) for i in range(len(motif[0,:]))]) + E_ms
var_E = np.sum([np.var(motif[:,i]) for i in range(len(motif[0,:]))])

#--------------------------Entropy function--------------------------
Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif, E_ms, L)
Kds = np.exp(Es[:-1])
beta_1, E_1, Kd_1 = get_p_properties(betas, Q0, Es, dE, 1)
#--------------------------Proofreading properties--------------------------
beta_step, E_step, K_step = get_proofreading_properties(betas, Q0, Es, dE, k_step, k_on)
print('beta_step = %.2f'%beta_step, 'K_step = %.2e'%K_step)

t0 = 1/lambda_A*np.log((lambda_A*N_A)/(k_on*b0))
print('--------')
#--------------------------Loops--------------------------

fig_antigen, ax_antigen = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
ax_antigen.plot(time_array, np.exp(lambda_A*time_array)/(1e0), linewidth = 5, color = antigen_color, ls = '--')

# Define the ODE: dN/dt = lambda * (1 - f(t)) * N
def dNdtNaive(t, N):
    pb = (1+(1e-9/(1e5*60*60*24*np.exp(lambda_B*(t))/N_A)))**(-1)
    return (lambda_A * (1 - pb) - 2*pb) * N
# Initial condition
N0 = 1.0
# Solve the ODE over the time span of your data
solNaive = solve_ivp(dNdtNaive, t_span=(time_array[0], time_array[-1]), y0=[N0], t_eval=time_array)
# Result
N_A_real = solNaive.y[0]  # solution N(t) evaluated at t_vals

ax_antigen.plot(time_array, N_A_real, linewidth = 5, color = antigen_color)

my_plot_layout(ax=ax_antigen, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
ax_antigen.set_xlim(right = Tf, left = T0)
ax_antigen.set_xticks([])
#ax_antigen.set_xlim(right = 1e-2, left = 1e-11) #use 1e-3 for other plots
ax_antigen.set_ylim(bottom = 1e2, top = 2e13)
#ax_antigen.legend(title = r'$\p$', title_fontsize = 34, fontsize = 32)
fig_antigen.savefig(output_plot + '/antigen.pdf')
plt.close(fig_antigen)

fig_antigen, ax_antigen = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
ax_antigen.vlines(4, 1, np.exp(lambda_A*4)/(1e0), linewidth = 5, color = antigen_color)
ax_antigen.plot(time_array[time_array>4], np.ones_like(time_array[time_array>4])*np.exp(lambda_A*4)/(1e0), linewidth = 5, color = antigen_color)

my_plot_layout(ax=ax_antigen, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
ax_antigen.set_xlim(right = Tf, left = T0)
ax_antigen.set_xticks([])
#ax_antigen.set_xlim(right = 1e-2, left = 1e-11) #use 1e-3 for other plots
ax_antigen.set_ylim(bottom = 1e2, top = 2e10)
#ax_antigen.set_ylim(bottom = 1, top = 1e7)
#ax_antigen.legend(title = r'$\p$', title_fontsize = 34, fontsize = 32)
fig_antigen.savefig(output_plot + '/antigen2.pdf')
plt.close(fig_antigen)

print('Loops...')
for L_0 in L_0s:
    print('________')
    print('L_0 = %.0e'%L_0)
    #--------------------------Repertoire properties--------------------------
    beta_r, E_r, Kd_r = get_repertoire_properties(betas, Q0, Es, dE, L_0)
    print('beta_r = %.1f'%beta_r, 'K_d_r = %.2e'%Kd_r)

    fig_R, ax_R =plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.2, 'top': 0.96})
    fig_K, ax_K = plt.subplots(figsize=(8*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_L, ax_L = plt.subplots(figsize=(8*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_N_b, ax_N_b = plt.subplots(figsize=(8.0*1.62,8), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

    U_threshold = 20
    N_threshold = U_threshold**(0.5*lambda_B/lambda_A)

    for i_p, p in enumerate(ps):
        fig_U_p, ax_U_p =plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.2, 'top': 0.96})
        fig_N_b_p, ax_N_b_p = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.2, 'top': 0.96})
        fig_N_b_p2, ax_N_b_p2 = plt.subplots(figsize=(8.0*1.62,8*0.6), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.2, 'top': 0.96})
        print('--------')
        print('p = %.2f...'%p)
        beta_p, E_p, Kd_p = get_p_properties(betas, Q0, Es, dE, p)
        print('beta_p = %.2f...'%p, 'K_p = %.2e...'%Kd_p)
        
        #-----------------Loading data----------------------------
        parameters_path = 'L-%d_Nbc-%d_Antigen-'%(L, L_0)+antigen+'_lambda_A-%.6f_lambda_B-%.6f_k_step-%.6f_theta-%.6f_Nc-%.6f_linear-%d_N_ens-%d_'%(lambda_A, 3.0, k_step/24, p, b0, linear, N_ens)+model

        #--------------------------m_bar(t)---------------------------
        u_on, p_a, R, QR = calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*time_array[0])/N_A, Es, p, lambda_A, b0, dE)
        M_r = L_0*b0*np.sum(Q0*p_a*dE)
        
        m_bar = np.array([np.sum(L_0*calculate_QR(Q0, k_on, k_step, np.exp(lambda_A*(t))/N_A, Es, p, lambda_A, b0, dE)[3]*dE) for t in time_array]) 
        m_bar_approx = ((k_on*M_r)/(N_A*lambda_A))*(np.exp(lambda_A*time_array))
        t_act = time_array[m_bar<1][-1]
        print('t_act:', t_act)

        ax_L.plot(time_array, m_bar, linewidth = 5, linestyle = '-', color = colors_p[i_p])

        #---------------------------- B cell linages ----------------------
        t_f = Tf
        n_t = int(m_bar[time_array<t_f][-1])
        activation_times = np.array([time_array[m_bar<=n][-1] for n in range(1, 5000, 5)])
        clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times + np.log(U_threshold)/(lambda_A), lambda_B, C, dT)
        #--------------------------t_C filter-------------------------
        lim_size = np.max([int(np.max(clone_sizes[:, -1])*0.01), 2])
        clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times + np.log(U_threshold)/(lambda_A), np.ones_like(activation_times), lim_size)
        bar_N_1 = clone_sizes_C[0, -1]

        ax_N_b.hlines(C, T0, Tf, color = 'grey', linestyle = ':')
        my_colors = [my_green_a, my_green_b, my_green_c]
        
        
        total_population = np.zeros_like(time_array)
        for i_k, k in enumerate([0, 1, 10]):
            t_1 = activation_times_C[k*10] - np.log(U_threshold)/(lambda_A) + 0.05
            t_Th = activation_times_C[k*10]
            print(t_1, t_Th)
            K_act =  K_step*np.exp((lambda_A/p)*(activation_times_C[k*10] - np.log(U_threshold)/(lambda_A) - t0))
            ax_U_p.plot(time_array, np.cumsum(N_A_real*dT)*b0*k_on/N_A*(K_step/(K_step+K_act))**p, linewidth = 3, color = my_colors[i_k], linestyle= '-', alpha = .8)
           
            # ax_N_b.plot(time_array, clone_sizes_C[k*10, :]-np.heaviside(activation_times_C[k*10] - time_array , 1), linewidth = 1.5, color = colors_p[i_p], linestyle= '-', alpha = .8)
            
            ax_N_b_p.plot(time_array[np.exp(lambda_B*0.5*(time_array - t_1))<=(N_threshold)], np.exp(lambda_B*0.5*(time_array[np.exp(lambda_B*0.5*(time_array - t_1))<=(N_threshold)] - t_1)), linewidth = 3, color = my_colors[i_k], linestyle= '-', alpha = .8)
            total_population[np.exp(lambda_B*0.5*(time_array - t_1))<=(N_threshold)] += np.exp(lambda_B*0.5*(time_array[np.exp(lambda_B*0.5*(time_array - t_1))<=(N_threshold)] - t_1)) * np.heaviside(time_array[np.exp(lambda_B*0.5*(time_array - t_1))<=(N_threshold)] - t_1 , 1)
            ax_N_b_p.plot(time_array[time_array>1.01*t_Th], 0.96*(N_threshold)*(clone_sizes_C[k*10, time_array>1.01*t_Th][:len(time_array[time_array>1.01*t_Th])]), linewidth = 3, color = my_colors[i_k], linestyle= '-', alpha = .8)
            total_population[time_array>1.01*t_Th] += 0.96*(N_threshold)*(clone_sizes_C[k*10, time_array>1.01*t_Th][:len(time_array[time_array>1.01*t_Th])]) * np.heaviside(time_array[time_array>1.01*t_Th] - t_Th , 1)

            ax_N_b_p2.plot(time_array, np.exp(lambda_B*3/activation_times_C[9*k*10]*(time_array - activation_times_C[1*10])), linewidth = 3, color = my_colors[i_k], linestyle= '-', alpha = .8)
        
        # ax_N_b_p.plot(time_array, total_population, linewidth = 4, color = 'darkgreen', linestyle= '-', alpha = .5)

        #--------------------------K_front(t)---------------------------
        K_front = (k_step/k_on)*np.exp(lambda_A*(time_array[time_array>t0] - t0)/p)
        t_p = time_array[time_array>t0][K_front>Kd_p][0]

        ax_K.plot(time_array[time_array>t0], K_front, alpha = .5, color = colors_p[i_p], linewidth = 5, linestyle = '--')
        ax_K.plot(time_array[time_array>1.15][0], Kd_r, markerfacecolor = 'grey', marker = '*', ms = 22, markeredgecolor='black', alpha = 0.8, zorder = 20)
        ax_K.plot(time_array[time_array>1.15][0], Kd_p, markerfacecolor = colors_p[i_p], marker = 'o', ms = 18, markeredgecolor='black', alpha = alpha_p[i_p], zorder = 20)

        if(p>beta_r):
            #ax_K.hlines(Kd_r, T0, t_act, color = colors_p[i_p], linewidth = 5, linestyle = '-', alpha = 1)
            ax_K.plot(time_array[time_array>=(t_act+0.2)], (Kd_r)*np.exp(lambda_A*(time_array[time_array>=(t_act+0.2)] - (t_act+0.2))/p), alpha = 1, color = colors_p[i_p], linewidth = 5, linestyle = '-')
        else:
            ax_K.hlines(Kd_p, t_act+0.15, t_p, color = colors_p[i_p], linewidth = 5, linestyle = '-', alpha = 1)
            ax_K.plot(time_array[time_array>=t_p], (Kd_p)*np.exp(lambda_A*(time_array[time_array>=t_p] - t_p)/p), alpha = 1, color = colors_p[i_p], linewidth = 5, linestyle = '-')

        ax_U_p.hlines([1, U_threshold], T0, Tf, color = 'grey', linestyle = [':', '--'], linewidth = 2)
        ax_N_b_p.hlines([(1, N_threshold)], T0, Tf, color = 'grey', linestyle = [':', '--'], linewidth = 2)


        my_plot_layout(ax=ax_N_b_p, yscale = 'log', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
        if p==1:
            ax_N_b_p.set_xticks([])
        ax_N_b_p.set_xlim(right = Tf, left = T0+1)
        ax_N_b_p.set_ylim(bottom = 1e0, top = 1e3)
        #ax_N_b.set_ylim(bottom = 1e0, top = C*1.1)
        fig_N_b_p.savefig(output_plot + '/Bcell_p-%.1f_L0-%.0e_'%(p, L_0)+energy_model+'.pdf')
        plt.close(fig_N_b_p)

        my_plot_layout(ax=ax_N_b_p2, yscale = 'log', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
        # ax_N_b_p2.set_xticks([])
        ax_N_b_p2.set_xlim(right = Tf, left = T0+1)
        ax_N_b_p2.set_ylim(bottom = 1e0, top = 5e2)
        #ax_N_b.set_ylim(bottom = 1e0, top = C*1.1)
        fig_N_b_p2.savefig(output_plot + '/Bcell2_p-%.1f_L0-%.0e_'%(p, L_0)+energy_model+'.pdf')
        plt.close(fig_N_b_p2)

        my_plot_layout(ax=ax_U_p, yscale = 'log', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
        ax_U_p.set_xticks([])
        ax_U_p.set_xlim(right = Tf, left = T0+1)
        ax_U_p.set_ylim(bottom = 1e-3, top = 2e5)
        fig_U_p.savefig(output_plot + '/U_p-%.1f_L0-%.0e_'%(p, L_0)+energy_model+'.pdf')
        plt.close(fig_U_p)

    # my_plot_layout(ax=ax_N_b, yscale = 'log', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    # ax_N_b.set_xlim(right = Tf-1, left = T0+1)
    # ax_N_b.set_ylim(bottom = 1e0, top = 8e2)
    # #ax_N_b.set_ylim(bottom = 1e0, top = C*1.1)
    # fig_N_b.savefig(output_plot + '/time/L%d/Bcell_Nr-%.0e_'%(L, L_0)+energy_model+'.pdf')
    # plt.close(fig_N_b)

    # my_plot_layout(ax=ax_K, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    # # ax_K.set_xticks([])
    # ax_K.set_xlim(right = Tf-1, left = T0+1)
    # ax_K.set_ylim(bottom = 9e-10, top = 2e-4)
    # fig_K.savefig(output_plot + '/time/L%d/K_Nr-%.0e_'%(L, L_0)+energy_model+'.pdf')
    # plt.close(fig_K)

    # my_plot_layout(ax=ax_L, yscale = 'log', xscale = 'linear', ticks_labelsize = 40, x_fontsize=30, y_fontsize=30 )
    # ax_L.set_xticks([])
    # ax_L.set_xlim(right = Tf-1, left = T0+1)
    # ax_L.set_ylim(bottom = 9e-1, top = 9e3)
    # ax_L.set_yticks([1e0, 1e2])
    # fig_L.savefig(output_plot + '/time/L%d/L_Nr-%.0e_'%(L, L_0)+energy_model+'.pdf')
    # plt.close(fig_L)





