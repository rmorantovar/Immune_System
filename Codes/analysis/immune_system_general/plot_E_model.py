import sys
sys.path.append('../library/')
from functions import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

#--------------- PARAMETERS ---------------------

my_E_m = -25
N_r = 1e9

print('min_E_PWM=%.2f'%(-24))
print('min_K_PWM=%.2e'%np.exp(-24))

print('--------')
#----------------------------------------------------------------
data_ready = 0
energy_models = ['TCRen', 'Gaussian', 'MJ2']
colors_models = [my_red, my_blue, my_green]
cmaps_model = ['Reds', 'Blues', 'Greens']

#fig_em, ax_em = plt.subplots(1, 2, figsize = (20, 8))
for i_model, energy_model in enumerate(energy_models):
    print('--------')
    print('Model = ', energy_model)
    print('--------')

    fig_Omega, ax_Omega = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    fig_betas, ax_betas = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    fig_l_sigma, ax_l_sigma = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    fig_K_r, ax_K_r = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    fig_K_m, ax_K_m = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    # fig4, ax4 = plt.subplots(figsize = (10, 8))
    # fig5, ax5 = plt.subplots(figsize = (10, 8))
    # fig6, ax6 = plt.subplots(figsize = (10, 8))
    # fig7, ax7 = plt.subplots(figsize = (10, 8))

    Motif, M, Alphabet = get_motif_em(antigen = 'R', energy_model = energy_model, M = [], read_matrix = True, Text_files_path = Text_files_path)

    df = pd.read_csv(Text_files_path + 'E_model/All_data_' + energy_model + '.csv')
    f = open(Text_files_path + 'E_model/trajectories_'+energy_model+'.pkl', 'rb')
    trajectories = pickle.load(f)
    l_trajectories = 10000-2

    E_trajectories = trajectories[0][:-(l_trajectories)].astype(float)
    rho_trajectories = trajectories[1][:-(l_trajectories)].astype(float)

    Es_average = trajectories[0][-(l_trajectories):].astype(float)
    rho_average = trajectories[1][-(l_trajectories):].astype(float)

    n_trajectories = int(len(trajectories[0][:-(l_trajectories)])/(l_trajectories))
    print(n_trajectories)
    for i in range(n_trajectories):
        Es_i = E_trajectories[i*(l_trajectories):(i+1)*(l_trajectories)]
        rho_i = np.array(rho_trajectories[i*(l_trajectories):(i+1)*(l_trajectories)])
        df_i = pd.DataFrame({'x': np.exp(Es_i[::1]), 'y':rho_i[::1]})
        sns.lineplot(data = df_i, x = 'x', y = 'y', lw = 1., alpha = .05, color = colors_models[i_model], ax = ax_Omega)
    
    df_avg = pd.DataFrame({'x': np.exp(Es_average[::1]), 'y':np.array(rho_average[::1])})
    sns.lineplot(data = df_avg, x = 'x', y = 'y', lw = 1.5, alpha = 1, color = colors_models[i_model], label = energy_model, ax = ax_Omega)
    
    sns.histplot(data = df, x = 'beta_r', ax = ax_betas, stat = 'probability', bins = np.linspace(1.6, 3.8, 30), color = colors_models[i_model])
    sns.histplot(data = df, x = 'K_r', ax = ax_K_r, stat = 'count', bins = np.logspace(-9, -6, 20), color = colors_models[i_model])
    sns.histplot(data = df, x = 'K_m', ax = ax_K_m, stat = 'count', bins = np.logspace(-13, -9, 20), color = colors_models[i_model])
    #ax_betas.vlines([np.mean(df['beta_r']), abs(np.log(5e-8))/(np.mean(df['l'])*np.mean(df['var'])**2)], 0, ax_betas.get_ylim()[1], color = colors_models[i_model], lw = 2, ls = ['-', '--'])
    #data2d = np.histogram2d(df['l'], df['var'], bins = [np.arange(15, 30), np.linspace(0, 2, 80)])
    #ax_l_sigma.contourf(data2d[2][:-1], data2d[1][:-1], data2d[0], cmap = cmaps_model[i_model], levels = 30)
    sns.histplot(data = df, x = 'var', y = 'l', discrete=(False, True),  bins = [np.linspace(.4, 1.4, 30), np.arange(14, 31)], cmap = cmaps_model[i_model], ax = ax_l_sigma)

   
    # #--------------------------Loops--------------------------
    # ax_Omega.scatter(np.exp(np.mean(np.log(K_ms))), (1/(20**l)), marker = 'o', color = colors_models[i_model], s = 28)
    # #ax_Omega.scatter(np.exp(np.mean(np.log(K_rs))), 10**rho_average[Es_average<(np.mean(np.log(K_rs)))][-1], marker = '*', color = colors_models[i_model], s = 28)

    my_plot_layout(ax=ax_Omega, xscale = 'log', yscale = 'log')
    ax_Omega.legend(fontsize = 20)
    ax_Omega.set_ylim(top = 1)
    ax_Omega.set_xlim(right = 1)
    fig_Omega.savefig('../../Figures/0_Shape_Space/energy_model_'+energy_model+'.pdf')

    # my_plot_layout(ax_Omega=ax_Omega, xscale = 'log', yscale = 'linear')
    # fig_Omega.savefig('../../Figures/0_Shape_Space/energy_model_2_'+energy_model+'.pdf')

    my_plot_layout(ax=ax_betas, xscale = 'linear', yscale = 'linear')
    ax_betas.set_xlim(left = 1.6, right = 3.8)
    #ax_betas.legend(fontsize = 20)
    #ax_Omega.set_xlim(left = np.exp(my_E_m-1))
    fig_betas.savefig('../../Figures/0_Shape_Space/betas_'+energy_model+'.pdf')

    my_plot_layout(ax=ax_K_r, xscale = 'log', yscale = 'linear')
    #ax_K_r.set_xlim(left = 1.5, right = 3.5)
    #ax_K_r.legend(fontsize = 20)
    #ax_K_r.set_xlim(left = np.exp(my_E_m-1))
    fig_K_r.savefig('../../Figures/0_Shape_Space/K_r_'+energy_model+'.pdf')

    my_plot_layout(ax=ax_K_m, xscale = 'log', yscale = 'linear')
    #ax_K_m.set_xlim(left = 1.5, right = 3.5)
    #ax_K_m.legend(fontsize = 20)
    #ax_K_m.set_xlim(left = np.exp(my_E_m-1))
    fig_K_m.savefig('../../Figures/0_Shape_Space/K_m_'+energy_model+'.pdf')


    my_plot_layout(ax=ax_l_sigma, xscale = 'linear', yscale = 'linear')
    #ax_l_sigma.set_xlim(left = 1.5, right = 3.5)
    #ax_l_sigma.legend(fontsize = 20)
    ax_l_sigma.set_yticks(range(17, 27, 3))
    ax_l_sigma.set_ylim(bottom = 17, top = 27)
    fig_l_sigma.savefig('../../Figures/0_Shape_Space/l_vs_std_'+energy_model+'.pdf')

    # my_plot_layout(ax_Omega=ax4, xscale = 'linear', yscale = 'linear')
    # ax4.set_xlim(left = 0.4, right = 1.6)
    # #ax4.legend(fontsize = 20)
    # #ax4.set_xlim(left = np.exp(my_E_m-1))
    # fig4.savefig('../../Figures/0_Shape_Space/alpha_'+energy_model+'.pdf')

    # my_plot_layout(ax_Omega=ax5, xscale = 'linear', yscale = 'linear')
    # ax5.set_xlim(left = 14, right = 30)
    # #ax5.legend(fontsize = 20)
    # #ax5.set_xlim(left = np.exp(my_E_m-1))
    # fig5.savefig('../../Figures/0_Shape_Space/ls_'+energy_model+'.pdf')

    # my_plot_layout(ax_Omega=ax6, xscale = 'linear', yscale = 'linear')
    # ax6.tick_params(labelsize = 24)
    # #ax6.set_xticks(np.arange(0, len(vars_data), 4) + 0.5)
    # #ax6.set_yticks(np.arange(0, len(ls_data)) + 0.5)
    # #ax6.set_xticklabels([r'$%.2f$'%i for i in np.linspace(0, 2, 20)])
    # #ax6.set_yticklabels([r'$%d$'%i for i in ls_data]);
    # #ax6.legend(fontsize = 20)
    # #ax6.set_xlim(left = np.exp(my_E_m-1))
    # fig6.savefig('../../Figures/0_Shape_Space/ls_vars_'+energy_model+'.pdf')

    # my_plot_layout(ax_Omega=ax7, xscale = 'linear', yscale = 'linear')
    # ax7.tick_params(labelsize = 24)
    # #ax7.set_xticks(np.arange(0, len(alphas_data), 2) + 0.5)
    # #ax7.set_yticks(np.arange(0, len(ls_data)) + 0.5)
    # #ax7.set_xticklabels([r'$%.2f$'%i for i in alphas[::2]])
    # #ax7.set_yticklabels([r'$%d$'%i for i in ls_data]);
    # #ax7.set_xlim(left = 14, right = 26)
    # #ax7.set_ylim(bottom = 0.4, top = 1.6)
    # #ax7.legend(fontsize = 20)
    # #ax7.set_xlim(left = np.exp(my_E_m-1))
    # fig7.savefig('../../Figures/0_Shape_Space/ls_alphas_'+energy_model+'.pdf')

    # #fig_em.savefig('../../Figures/0_Shape_Space/motifs.pdf')

