import sys
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/immune_system/primary_immune_response/'

# Define constants and parameters
D0 = 3*1e-8 * (60*60*24) # Diffusion coefficient [cm^2/days]
k_on = 1e6; #(M*days)^-1
Kstep = (0.5/60)/k_on;
Kstar = 10**-8;
b0 = 10**5
p = 3
u = .000  # Velocity
L = 10 # Length of the domain [cm]
T = 100  # Total simulation time [days]
Nx = 1000  # Number of spatial grid points
Nt = 10000  # Number of time steps
lamA = 6
lamB = 3 * np.log(2) #(days)^-1


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
# D0 = 1e-2
Ds = np.logspace(np.log10(0.01*D0), np.log10(100*D0), 12)
Nstars = [1e10/100, 1e10/1]
r0s = np.linspace(0.1, 1, 2)


for immunization in immunizations:
    print(immunization)
    fig_lam, ax_lam = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_t, ax_t = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

    results = defaultdict(list)
    for r0 in r0s:
        for i_D, D in enumerate(tqdm(Ds)):
            Nstar = ((1/N_A)*k_on*b0*(Kstar/Kstep)**(-p)/(lamA/(60*60*24)))**-1
            r_array = np.linspace(0, L, Nx+1)
            growth_rate = lamA
            NAeff = np.zeros_like(tGrid)
            for i_t in range(1, Nt+1):
                C[i_t,:] = np.exp(lamA*i_t*dt)*(1/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
                #C[C<1]=0
                NAeff[i_t] = np.sum(C[i_t,xGrid>r0]*r_array[xGrid>r0]**2*dx)*(4*np.pi)
            lamAeff_array = np.diff(np.log(NAeff))/np.diff(tGrid)
            R = 1-np.exp(-np.cumsum((1/N_A)*k_on*(60*60*24)*b0*(Kstar/Kstep)**(-p)*NAeff*dt))
            # popt, pcov = curve_fit(f = my_linear_func, xdata = tGrid[NAeff>Nstar/2][:5], ydata = np.log(NAeff)[NAeff>Nstar/2][:5])
            # tstar = tGrid[NAeff>(Nstar)][0]
            # lamAeff = lamAeff_array[NAeff[:-1]>Nstar][0]
            tstar = tGrid[R<(0.5)][-1]
            lamAeff = lamAeff_array[R[:-1]<0.5][-1]
            rstar = np.sqrt(D*tstar2)
            
            results['N'].append(Nstar)
            results['D'].append(D)
            results['r0'].append(r0)
            results['t'].append(tstar)
            results['lamA'].append(lamAeff)
            results['alpha'].append(r0/rstar)
            #print(popt[1])

    results_db = pd.DataFrame(results)
    # print(results_db)
        
    alpha_array = np.logspace(np.log10(0.1), np.log10(40), 100)
    
    markers = ['o', 's', '^', 'D', '*']
    for i_r0, r0 in enumerate(r0s):
        data_r0 = results_db.loc[results_db['r0']==r0]
        alpha_D_r_array = np.logspace(np.min(r0**(1/2)*np.log10(data_r0['alpha']/data_r0['D']**(1/4))), np.log10(3), 100)
        
        print(r0, np.array(data_r0['alpha'][data_r0['D']<3e-3])[-1], np.array(data_r0['t'][data_r0['D']<3e-3])[-1], np.array(data_r0['lamA'][data_r0['D']<3e-3])[-1]) 

        # ax_t.plot(data_r0['alpha'], data_r0['t']/((1/lamA)*np.log(Nstar)), lw = 5, ls = ''
            # , marker = markers[i_r0], ms = 8, color = my_purple, label = r'$%.1f$'%(r0))
        ax_t.plot(data_r0['alpha'], data_r0['t']/((1/lamA)*np.log(Nstar)), lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'black', label = r'$%.1f$'%(r0))

        # ax_lam.plot(data_r0['alpha'], (data_r0['lamA']/lamA), lw = 5, ls = ''
            # , marker = markers[i_r0], ms = 8, color = my_purple, label = r'$%.1f$'%(r0))
        ax_lam.plot(data_r0['alpha'], (data_r0['lamA']/lamA), lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'black', label = r'$%.1f$'%(r0))        

        ax.plot(data_r0['alpha'], lamB*p/(data_r0['lamA']*2.2), lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'black', label = r'$%.1f$'%(r0))
    
    
    ax_t.fill_between(alpha_array[(alpha_array>1.16) & (alpha_array<8.54)], 0.5*np.ones_like(alpha_array[(alpha_array>1.16) & (alpha_array<8.54)])
        , 7*np.ones_like(alpha_array[(alpha_array>1.16) & (alpha_array<8.54)]), color = my_purple2, alpha = .2)
    ax_t.vlines(10.78, 0.5, 7, color = 'k', ls = '--')
    ax_t.plot(alpha_array, 1 + 1/(4*np.log(Nstar))*alpha_array**(2) - 2*np.log(alpha_array)/np.log(Nstar) + np.log(4*3.1415)/(2*np.log(Nstar))
        , color = my_purple2, lw = 5, alpha = .6)
    #ax_t.plot(alpha_array, 0.01*alpha_array**2
    #   , color = my_purple2, lw = 5, alpha = .6). # To check scale behaviour a^2 for a>tilde a
    my_plot_layout(ax = ax_t, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_t.set_xlim(left = 0, right = 7)
    ax_t.set_ylim(bottom = 0.5, top = 7)
    # ax_t.set_yticks(range(1, 12, 3))
    ax_t.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in C0s])
    fig_t.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_' + immunization + '.pdf')
    plt.close()

    ax_lam.fill_between(alpha_array[(alpha_array>1.16) & (alpha_array<8.54)], 0.8*np.ones_like(alpha_array[(alpha_array>1.16) & (alpha_array<8.54)])
        , 2.05*np.ones_like(alpha_array[(alpha_array>1.16) & (alpha_array<8.54)]), color = my_purple2, alpha = .2)
    ax_lam.vlines(10.78, 0.8, 2.05, color = 'k', ls = '--')
    ax_lam.plot(alpha_array, 1+alpha_array**2/((4*np.log(Nstar) + 4*np.log(alpha_array**2) + 2*np.log(4*3.14)) + alpha_array**2)
        , color = my_purple2, lw = 5, alpha = .6)
    # ax_lam.plot(alpha_array, 1+alpha_array**2/((4*np.log(Nstar) + 4*np.log(alpha_array**2) + 2*np.log(4*3.14)) + alpha_array**2)
        # , color = my_purple2, lw = 5, alpha = .6)
    my_plot_layout(ax = ax_lam, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_lam.set_xlim(left = 0, right = 7)
    ax_lam.set_ylim(bottom = 0.8, top = 2.05)
    ax_lam.set_yticks(range(1, 3))
    ax_lam.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in C0s])
    fig_lam.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_' + immunization + '.pdf')

    ax.plot(alpha_array, lamB*p/(2.2*lamA*(1+alpha_array**2/((4*np.log(Nstar) + 4*np.log(alpha_array**2) + 2*np.log(4*3.14)) + alpha_array**2)))
        , color = my_purple2, lw = 5, alpha = .6)
    my_plot_layout(ax = ax, xscale='log', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax.set_xlim(left = 0, right = 7)
    #ax.set_ylim(bottom = 0.8)
    ax.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in C0s])
    fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/zeta_' + immunization + '.pdf')


  



  