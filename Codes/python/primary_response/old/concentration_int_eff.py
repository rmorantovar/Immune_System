import sys
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/immune_system/primary_immune_response/'

# Define constants and parameters
D0 = 3*1e-8 * (60*60*24) # Diffusion coefficient [cm^2/days]
k_on = 1e6*24*3600; #(M*days)^-1
b0 = 10**5
u = .000  # Velocity
L = 10 # Length of the domain [cm]
T = 40  # Total simulation time [days]
Nx = 1000  # Number of spatial grid points
Nt = 10000  # Number of time steps
lamA = 6
r0 = 0.3 # cm

# Discretize the spatial and temporal dimensions
dx = L / Nx
dt = T / Nt

# Create spatial and temporal grids
xGrid = np.linspace(0, L, Nx + 1)
tGrid = np.linspace(0, T, Nt + 1)

# Initialize concentration matrix
C = np.zeros((Nt + 1, Nx + 1))

immunizations = ['infec']#, 'vacc']
colorsD = [my_purple, my_blue, my_green, my_red]
D0 = 1e-3
Ds = np.logspace(np.log10(0.1*D0), np.log10(100*D0), 20)


for immunization in immunizations:
    print(immunization)
    if immunization == 'vacc':
        C0s = np.logspace(12., 14., 3)
    else:
        C0s = [1]
    results = defaultdict(list)
    for C0 in C0s:
        for i_D, D in enumerate(tqdm(Ds)):
            r1 = np.sqrt((2*D/lamA)*np.log((lamA*N_A)/(0.8*k_on*b0)))
            rho1 = 1e10/10#(lamA*N_A/(k_on*b0))*10
            #print('D=%.0e'%D, 'r1 = %.2f'%r1, 'Na1 = %.1e'%rho1)
            r_array = np.linspace(0, L, Nx+1)
            if immunization == 'vacc':
                #C[0,0] = C0
                integral_in = np.zeros_like(tGrid)
                for i_t in range(1, Nt+1):
                    C[i_t,:] = (C0/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
                    #C[C<1]=0
                    integral_in[i_t] = np.sum(C[i_t,xGrid>r0]*r_array[xGrid>r0]**2*dx)*(4*np.pi)
            else:
                growth_rate = lamA
                #C[0,0] = C0
                integral_in = np.zeros_like(tGrid)
                for i_t in range(1, Nt+1):
                    C[i_t,:] = np.exp(lamA*i_t*dt)*(1/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
                    #C[C<1]=0
                    integral_in[i_t] = np.sum(C[i_t,xGrid>r0]*r_array[xGrid>r0]**2*dx)*(4*np.pi)
                # print(tGrid[integral_in<rho1][-10:])
                # print((integral_in)[integral_in<rho1][-10:])
            popt, pcov = curve_fit(f = my_linear_func, xdata = tGrid[integral_in>rho1][:15], ydata = np.log(integral_in)[integral_in>rho1][:15])
            tstar = tGrid[integral_in>(10*rho1)][0]
            rstar = np.sqrt(2*D*tstar)
            results['C0'].append(C0)
            results['D'].append(D)
            results['lamA'].append(popt[1])
            results['t'].append(tstar)
            results['r'].append(rstar/r0)
            #print(popt[1])

    results_db = pd.DataFrame(results)
    # print(results_db)

    if immunization == 'vacc':
        fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
        for D in Ds:
            ax.plot(results_db.groupby('C0')['r'].mean()/D , results_db.loc[results_db['D']==D]['lamA']/lamA, color = my_green, lw = 5)
        #ax.plot(np.logspace(np.log10(3e0), np.log10(3e1)) , 2*np.logspace(np.log10(3e0), np.log10(3e1))**(-2)*lamA, color = my_green, lw = 4, ls = '--')
        my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
        #ax.set_xlim(left = 0, right = 7)
        #ax.set_ylim(bottom = 1e0, top = 2e11)
        #ax.legend(fontsize = 20, title_fontsize = 24, loc = 0)#, labels = [r'$%.1e$'%i for i in C0s])
        fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_A_eff_' + immunization + '.pdf')

        fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
        for D in Ds:
            ax.plot(results_db.groupby('C0')['r'].mean()/D , results_db.loc[results_db['D']==D]['t'], color = my_green, lw = 5)
        #ax.plot(np.logspace(np.log10(3e0), np.log10(3e1)) , 2*np.logspace(np.log10(3e0), np.log10(3e1))**(-2)*lamA, color = my_green, lw = 4, ls = '--')
        my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
        #ax.set_xlim(left = 0, right = 7)
        #ax.set_ylim(bottom = 1e0, top = 2e11)
        #ax.legend(fontsize = 20, title_fontsize = 24, loc = 0)#, labels = [r'$%.1e$'%i for i in C0s])
        fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_' + immunization + '.pdf')

        plt.close()
    else:
        fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
        ax.plot(results_db['r'], results_db['lamA']-lamA, color = my_purple, lw = 5, ls = '', marker = 'o', ms = 8)
        ax.plot(results_db['r'], .1*results_db['r']**(-3), color = my_purple2)
        ax.plot(results_db['r'], .1*results_db['r']**(-2), color = my_purple2)
        my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
        #ax.set_xlim(left = 0, right = 7)
        #ax.set_ylim(bottom = 1e0, top = 2e11)
        #ax.legend(fontsize = 20, title_fontsize = 24, loc = 0)#, labels = [r'$%.1e$'%i for i in C0s])
        fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_A_eff_' + immunization + '.pdf')
        plt.close()

        fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
        ax.plot(results_db['r'], results_db['t'] - (1/lamA)*np.log(rho1*10), color = my_purple, lw = 5, ls = '', marker = 'o', ms = 8)
        ax.plot(results_db['r'], 0.05*results_db['r']**(-3), color = my_purple2)
        ax.plot(results_db['r'], 0.05*results_db['r']**(-2), color = my_purple2)
        my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
        #ax.set_xlim(left = 0, right = 7)
        #ax.set_ylim(bottom = 1e0, top = 2e11)
        #ax.legend(fontsize = 20, title_fontsize = 24, loc = 0)#, labels = [r'$%.1e$'%i for i in C0s])
        fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_' + immunization + '.pdf')
        plt.close()

        #ax.plot(Ds, lam_eff, lw = 4, ls = '-', label = r'$%.1e$'%C0)#, color = antigen_color)
        # ax2.plot(tGrid[:], integral_out, ls = '--', lw = 5, color = antigen_color, label = 'outer')

    # X, Y = np.meshgrid(Ds, C0s)
    # results_db = pd.DataFrame(results)
    # df_pivot_lamA = results_db.pivot_table(index='C0', 
    #                       columns='D', 
    #                       values='lamA', 
    #                       aggfunc='mean')
    # df_pivot_t = results_db.pivot_table(index='C0', 
    #                       columns='D', 
    #                       values='t', 
    #                       aggfunc='mean')
    # df_pivot_r = results_db.pivot_table(index='C0', 
    #                       columns='D', 
    #                       values='r', 
    #                       aggfunc='mean')

    # fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    # cs = ax.pcolor(X, Y, df_pivot_lamA, shading = 'nearest', edgecolors = 'face', cmap = 'plasma', vmin=0, vmax=20)
    # cbar = fig.colorbar(cs)
    # cbar.ax.tick_params(labelsize=30) 
    # if immunization == 'vacc':
    #     cbar.ax.set_yticks(range(2, 20, 2))
    # else:
    #     cbar.ax.set_yticks(range(2, 20, 2))
    # # ax.hlines(lamA, np.min(Ds), np.max(Ds))
    # my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    # #ax.set_xlim(left = 0, right = 7)
    # #ax.set_ylim(bottom = 1e0, top = 2e11)
    # #ax.legend(fontsize = 20, title_fontsize = 24, loc = 0)
    # fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_A_eff_' + immunization + '.pdf')
    # plt.close()

    # fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    # cs = ax.pcolor(X, Y, df_pivot_t, shading = 'nearest', edgecolors = 'face', cmap = 'plasma', vmin=0, vmax=41)
    # cbar = fig.colorbar(cs)
    # cbar.ax.tick_params(labelsize=30) 
    # if immunization == 'vacc':
    #     cbar.ax.set_yticks(range(0, 41, 5))
    # else:
    #     cbar.ax.set_yticks(range(0, 41, 5))
    # # ax.hlines(lamA, np.min(Ds), np.max(Ds))
    # my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    # #ax.set_xlim(left = 0, right = 7)
    # #ax.set_ylim(bottom = 1e0, top = 2e11)
    # #ax.legend(fontsize = 20, title_fontsize = 24, loc = 0)
    # fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_' + immunization + '.pdf')
    # plt.close()
    
    # fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    # cs = ax.pcolor(X, Y, df_pivot_r, shading = 'nearest', edgecolors = 'face', cmap = 'plasma', vmin=0, vmax=.9)
    # cbar = fig.colorbar(cs)
    # cbar.ax.tick_params(labelsize=30) 
    # if immunization == 'vacc':
    #     cbar.ax.set_yticks(np.linspace(.1, .8, 8))
    # else:
    #     cbar.ax.set_yticks(np.linspace(.1, .8, 8))
    # # ax.hlines(lamA, np.min(Ds), np.max(Ds))
    # my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    # #ax.set_xlim(left = 0, right = 7)
    # #ax.set_ylim(bottom = 1e0, top = 2e11)
    # #ax.legend(fontsize = 20, title_fontsize = 24, loc = 0)
    # fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/r_' + immunization + '.pdf')






