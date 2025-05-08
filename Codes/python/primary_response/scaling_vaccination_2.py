import sys
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")
from scipy.optimize import fsolve

def equation(x, y):
    return  (x**-2) + np.log(x**-8) - y

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/immune_system/primary_immune_response/'

# Define constants and parameters
D0 = 3e-8*(60*60*24) # Diffusion coefficient [cm^2/day]
k_on = 1e6*(60*60*24); #(M*day)^-1
Kstep = (0.5*(60*24))/(k_on);
Kstar = 3*10**-8;
b0 = 10**5
p = 3
u = .000  # Velocity
L = 10 # Length of the domain [cm]
T = 400  # Total simulation time [days]
Nx = 1000  # Number of spatial grid points
Nt = 40000  # Number of time steps
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

immunizations = ['vacc']
colorsD = [my_purple, my_blue, my_red, my_green, my_purple, my_purple2]
# D0 = 1e-2
Ds = [D0]
# Ds = np.logspace(np.log10(0.01*D0), np.log10(100*D0), 5)
r0s = [0.1, 1.0]
N0s = np.logspace(8, 19, 11)
#N0s = [1e13]

uact = ((1/N_A)*k_on*b0*(Kstar/Kstep)**(-p))
print('%.1e, %.1e'%(Kstar, Kstep))

for immunization in immunizations:
    print(immunization)
    fig0, ax0 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

    fig_lam, ax_lam = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_t, ax_t = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

    fig_lam_N0, ax_lam_N0 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_t_N0, ax_t_N0 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_N0, ax_N0 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})

    results = defaultdict(list)
    for i_r0, r0 in enumerate(r0s):
        for i_D, D in enumerate(Ds):
            for N0 in tqdm(N0s):
                r_array = np.linspace(0, L, Nx+1)
                NAeff = np.zeros_like(tGrid)
                for i_t in range(1, Nt+1):
                    C[i_t,:] = (N0/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
                    #C[C<1]=0
                    NAeff[i_t] = np.sum(C[i_t,xGrid>r0]*r_array[xGrid>r0]**2*dx)*(4*np.pi)
                lamAeff_array = np.diff(np.log(NAeff))/np.diff(tGrid)
                R = 1-np.exp(-np.cumsum(uact*NAeff*dt))
                # Rapprox = 1-np.exp(-Nstar**-1*NAeff)
                tstar = tGrid[R<(0.5)][-1]
                lamAeff = lamAeff_array[R[:-1]<0.5][-1]
                # tstar_approx = tGrid[Rapprox<(0.5)][-1]
                # lamAeff_approx = lamAeff_array[Rapprox[:-1]<0.5][-1]
                rstar = np.sqrt(D*tstar)
                # rstar_approx = np.sqrt(D*tstar_approx)
    
                results['N0'].append(N0)
                results['D'].append(D)
                results['r0'].append(r0)
                results['t'].append(tstar)
                results['lamA'].append(lamAeff)
                results['alpha'].append(r0/rstar)

                # print(r0/rstar, 1/fsolve(equation, x0=0.00001, args=(np.log((N0*uact*4*r0**2)**4/(D*((4*3.1415)**1/2))**4),)))


    results_db = pd.DataFrame(results)
    N0_array = np.logspace(5, 19., 10000)
    
    
    markers = ['o', 's', '^', 'D', '*']

    for i_r0, r0 in enumerate(r0s):
        data_r0 = results_db.loc[(results_db['r0']==r0) & (results_db['D']==D0)]
        print(r0, np.array(data_r0['alpha'][data_r0['N0']<1e13])[-1], np.array(data_r0['t'][data_r0['N0']<1e13])[-1], np.array(data_r0['lamA'][data_r0['N0']<1e13])[-1]) 
        
        alpha_array = np.logspace(-1, np.log10(40), 100)

        factor = D0/r0**2

        ax_t.plot(data_r0['alpha'], data_r0['t']*factor, lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'k', label = r'$%.1f$'%(r0))
        # ax_t.plot(data_r0['alpha2'], data_r0['t2']*factor, lw = 5, ls = ''
            # , marker = markers[i_r0], ms = 8, color = 'grey')
        ax_t.plot(alpha_array, (alpha_array)**-2
        , color = my_blue2, lw = 5, alpha = .6)

        ax_lam.plot(data_r0['alpha'], data_r0['lamA']/factor, lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'k', label = r'$%.1f$'%(r0))
        # ax_lam.plot(data_r0['alpha2'], data_r0['lamA2']/factor, lw = 5, ls = ''
            # , marker = markers[i_r0], ms = 8, color = 'grey')
        ax_lam.plot(alpha_array, (alpha_array)**(4)/4
        , color = my_blue2, lw = 5, alpha = .6)

        ax.plot(data_r0['alpha'], lamB*p/(data_r0['lamA']*2.2), lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'k', label = r'$%.1f$'%(r0))
        # ax.plot(data_r0['alpha2'], lamB*p/(data_r0['lamA2']*2.2), lw = 5, ls = ''
            # , marker = markers[i_r0], ms = 8, color = 'grey')
        ax.plot(alpha_array, lamB*p/(2.2*(D*(alpha_array)**(4)/(4*r0**2)))
        , color = my_blue2, lw = 5, alpha = .6)

        # ---- For dosage plot ----

        alpha_array_N00 = np.sqrt(((np.log((N0_array*uact*4*r0**2/D0)**4/(4*3.1415)**2))))
        alpha_array_N0 = np.array([1/fsolve(equation, x0=0.001, args=(np.log((N0i*uact*4*r0**2)**4/(D*((4*3.1415)**1/2))**4),)) for N0i in N0_array])

        # alpha_array = np.linspace(1, 8, 100)
        N01 = np.exp(1**2/4)*np.sqrt(2*np.pi)*(1/uact)*D0/(4*r0**2)*1**2
        #print(N01, N0_array[alpha_array_N0>1][0])

        ax_t_N0.plot(data_r0['N0'], data_r0['t']*factor, lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'k', label = r'$%.1f$'%(r0))
        # ax_t_N0.plot(data_r0['N0'], data_r0['t2']*factor, lw = 5, ls = ''
            # , marker = markers[i_r0], ms = 8, color = 'grey')
        ax_t_N0.plot(N0_array, (alpha_array_N0)**(-2)
        , color = my_blue2, lw = 5, alpha = .6)
        ax_t_N0.vlines(N01, 5e-3, alpha_array_N0[N0_array<N01][-1]**(-2), color = 'k', ls = '--')
        # ax_t_N0.plot(N0_array, (alpha_array_N00)**(-2)
        # , color = my_blue2, lw = 5, alpha = .6, ls = '--')

        ax_lam_N0.plot(data_r0['N0'], (data_r0['lamA']/factor), lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'k', label = r'$%.1f$'%(r0))
        # ax_lam_N0.plot(data_r0['N0'], (data_r0['lamA2']/factor), lw = 5, ls = ''
            # , marker = markers[i_r0], ms = 8, color = 'grey')
        ax_lam_N0.plot(N0_array, (alpha_array_N0)**(4)/4
        , color = my_blue2, lw = 5, alpha = .6)
        ax_lam_N0.vlines(N01, 2e-4, alpha_array_N0[N0_array<N01][-1]**(4)/4, color = 'k', ls = '--')
        # ax_lam_N0.plot(N0_array, (alpha_array_N00)**(4)/4
        # , color = my_blue2, lw = 5, alpha = .6, ls = '--')

        ax_N0.plot(data_r0['N0'], lamB*p/(data_r0['lamA']*2.2), lw = 5, ls = ''
            , marker = markers[i_r0], ms = 8, color = 'k', label = r'$%.1f$'%(r0))
        # ax_N0.plot(data_r0['N0'], lamB*p/(data_r0['lamA2']*2.2), lw = 5, ls = ''
            # , marker = markers[i_r0], ms = 8, color = 'grey')
        ax_N0.plot(N0_array, lamB*p/(2.2*(D*(alpha_array_N0)**(4)/(4*r0**2)))
        , color = my_blue2, lw = 5, alpha = .6)
        # ax_N0.plot(N0_array, lamB*p/(2.2*(D*(alpha_array_N00)**(4)/(4*r0**2)))
        # , color = my_blue2, lw = 5, alpha = .6, ls = '--')
    

    ax0.plot(N0_array, alpha_array_N0, ls = '-',
     color = my_blue2, lw = 5, alpha = .6)
    ax0.plot(N0_array, alpha_array_N00, ls = '--', ms = 8, marker = '',
     color = my_blue, lw = 5, alpha = .6)

    # ax0.vlines(0.273548, 0, 1)
    #ax0.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    my_plot_layout(ax = ax0, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax0.set_xlim(left = 0.8)
    #ax0.set_ylim(bottom = -0.2, top = 3)
    #ax0.legend(fontsize = 16, title_fontsize = 20, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig0.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/alpha_' + immunization + '.pdf')

    ax_t.fill_between(alpha_array[(alpha_array>3.85) & (alpha_array<5.5)], 5e-4*np.ones_like(alpha_array[(alpha_array>3.85) & (alpha_array<5.5)])
        , 2e2*np.ones_like(alpha_array[(alpha_array>3.85) & (alpha_array<5.5)]), color = my_purple2, alpha = .2)
    ax_t.vlines(1, 5e-4, 2e2, color = 'k', ls = '--')
    my_plot_layout(ax = ax_t, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    # ax_t.set_xlim(left = 0.8)
    ax_t.set_ylim(bottom = 5e-4, top = 2e2)
    #ax_t.set_yticks(range(1, 12, 3))
    ax_t.legend(fontsize = 16, title_fontsize = 20, loc = 2, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig_t.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_' + immunization + '.pdf')
    plt.close()

    ax_lam.fill_between(alpha_array[(alpha_array>3.85) & (alpha_array<5.5)], 5e-6*np.ones_like(alpha_array[(alpha_array>3.85) & (alpha_array<5.5)])
        , 5e5*np.ones_like(alpha_array[(alpha_array>3.85) & (alpha_array<5.5)]), color = my_purple2, alpha = .2)
    ax_lam.vlines(1, 5e-6, 5e5, color = 'k', ls = '--')
    my_plot_layout(ax = ax_lam, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    # ax_lam.set_xlim(left = 0.8)
    ax_lam.set_ylim(bottom = 5e-6, top = 5e5)
    ax_lam.legend(fontsize = 16, title_fontsize = 20, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig_lam.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_' + immunization + '.pdf')
    
    my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax.set_xlim(left = 0, right = 7)
    #ax.set_ylim(bottom = 0.8)
    ax.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in C0s])
    fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/zeta_' + immunization + '.pdf')

    # ax_t_N0.fill_between(N0_array[(N0_array>1e12) & (N0_array<1e14)], 5e-3*np.ones_like(N0_array[(N0_array>1e12) & (N0_array<1e14)])
        # , 1e1*np.ones_like(N0_array[(N0_array>1e12) & (N0_array<1e14)]), color = my_purple2, alpha = .2)
    
    my_plot_layout(ax = ax_t_N0, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_t_N0.set_xlim(left = 8e3)
    ax_t_N0.set_ylim(bottom = 5e-3, top = 1.5e1)
    #ax_t.set_yticks(range(1, 12, 3))
    ax_t_N0.legend(fontsize = 16, title_fontsize = 20, loc = 2, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig_t_N0.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_N0_' + immunization + '.pdf')
    plt.close()

    # ax_lam_N0.fill_between(N0_array[(N0_array>1e12) & (N0_array<1e14)], 5e-3*np.ones_like(N0_array[(N0_array>1e12) & (N0_array<1e14)])
        # , 2e3*np.ones_like(N0_array[(N0_array>1e12) & (N0_array<1e14)]), color = my_purple2, alpha = .2)
    my_plot_layout(ax = ax_lam_N0, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax_lam_N0.set_xlim(left = 8e3)
    ax_lam_N0.set_ylim(bottom = 2e-4, top = 1e4)
    ax_lam_N0.legend(fontsize = 16, title_fontsize = 20, loc = 2, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig_lam_N0.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_N0_' + immunization + '.pdf')
    
    my_plot_layout(ax = ax_N0, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax.set_xlim(left = 0, right = 7)
    #ax.set_ylim(bottom = 0.8)
    ax_N0.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in C0s])
    fig_N0.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/zeta_N0_' + immunization + '.pdf')






