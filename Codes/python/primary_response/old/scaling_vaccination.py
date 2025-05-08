import sys
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")
from scipy.optimize import fsolve

def equation(x, y):
    return y + 1/(2*x**2) + 2*np.log(x)

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/immune_system/primary_immune_response/'

# Define constants and parameters
D0 = 3*1e-8 * (60*60*24) # Diffusion coefficient [cm^2/days]
k_on = 1e6*24*3600; #(M*days)^-1
b0 = 10**5
u = .000  # Velocity
L = 10 # Length of the domain [cm]
T = 200  # Total simulation time [days]
Nx = 100  # Number of spatial grid points
Nt = 20000  # Number of time steps
lamA = 6
r0 = 0.3

# Discretize the spatial and temporal dimensions
dx = L / Nx
dt = T / Nt

# Create spatial and temporal grids
xGrid = np.linspace(0, L, Nx + 1)
tGrid = np.linspace(0, T, Nt + 1)

# Initialize concentration matrix
C = np.zeros((Nt + 1, Nx + 1))

immunizations = ['vacc']
colorsD = [my_purple, my_blue, my_red, my_green]
D0 = 1e-2
Ds = np.logspace(np.log10(0.5*D0), np.log10(2*D0), 3)
Ds = [D0]
r0s = np.linspace(0.1, 1.0, 3)
r0s = [0.2, 0.5, 1]
N0s = np.logspace(10.5, 17., 10)
N0 = 1e13

for immunization in immunizations:
    print(immunization)
    fig_lam, ax_lam = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_t, ax_t = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    for i_D, D in enumerate(Ds):
        results = defaultdict(list)
        for i_r0, r0 in enumerate(r0s):
            for N0 in tqdm(N0s):
                Nstar = 1e10
                r_array = np.linspace(0, L, Nx+1)
                NAeff = np.zeros_like(tGrid)
                for i_t in range(1, Nt+1):
                    C[i_t,:] = (N0/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
                    #C[C<1]=0
                    NAeff[i_t] = np.sum(C[i_t,xGrid>r0]*r_array[xGrid>r0]**2*dx)*(4*np.pi)
                
                popt, pcov = curve_fit(f = my_linear_func, xdata = tGrid[NAeff>Nstar*0.95][:10], ydata = np.log(NAeff)[NAeff>Nstar*0.95][:10])
                tstar = tGrid[NAeff>(Nstar)][0]
                rstar = np.sqrt(2*D*tstar)
                results['N0'].append(N0)
                results['D'].append(D)
                results['lamA'].append(popt[1])
                results['r0'].append(r0)
                results['t'].append(tstar)
                results['alpha'].append(rstar/r0)
                #print(popt[1])

        results_db = pd.DataFrame(results)
        N0_array = np.logspace(np.log10(1.05*Nstar*np.sqrt(3.1415/2)), 18, 10000)
        alpha_array = np.array([fsolve(equation, x0=0.0001, args=(np.log((Nstar/y)*(3.1415/2)**(1/2)),)) for y in N0_array])
        alpha_array = 1/alpha_array
        alpha_array = np.logspace(np.log10(0.065), np.log10(10), 100)
        # alpha_array = np.sqrt(1/((2)*(np.log((N0_array/Nstar)*(2/3.1415)**(1/2)))))
        # alpha_array2 = np.sqrt(((((N0_array/Nstar)*(2/3.1415)**(1/2)))))
        #alpha_array2 = alpha_array**2
        markers = ['o', 's', '^', 'D', '*']
        # for i_D, D in enumerate(Ds):
        for i_r0, r0 in enumerate(r0s):
            data_r0 = results_db.loc[results_db['r0']==r0]

            factor = 2*D/r0**2

            ax_lam.plot(1/data_r0['alpha'], (data_r0['lamA']/factor), lw = 5, ls = ''
                , marker = markers[i_r0], ms = 8, color = my_blue, label = r'$%.1f$'%(r0))
            ax_lam.plot(1/alpha_array, (alpha_array)**(-4)/2
            , color = my_blue2, lw = 5, alpha = .6)
            
            ax_t.plot(1/data_r0['alpha'], data_r0['t']*factor, lw = 5, ls = ''
                , marker = markers[i_r0], ms = 8, color = my_blue, label = r'$%.1f$'%(r0))
            ax_t.plot(1/alpha_array, (alpha_array)**2
            , color = my_blue2, lw = 5, alpha = .6)
            

            # For dosage plot

            # ax_lam.plot(data_r0['N0'], (data_r0['lamA']/factor), lw = 5, ls = ''
            #     , marker = markers[i_r0], ms = 8, color = colorsD[i_D], label = r'$%.1f$'%(r0))
            # ax_lam.plot(N0_array, (alpha_array)**(4)/2
            # , color = my_blue2, lw = 5, alpha = .6)
            
            # ax_t.plot(data_r0['N0'], data_r0['t']*factor, lw = 5, ls = ''
            #     , marker = markers[i_r0], ms = 8, color = colorsD[i_D], label = r'$%.1f$'%(r0))
            # ax_t.plot(N0_array, (alpha_array)**(-2)
            # , color = my_blue2, lw = 5, alpha = .6)
    

    ax.plot(N0_array, (np.sqrt(((2)*(np.log((N0_array/Nstar)*(2/3.1415)**(1/2)))))),
     color = my_blue2, lw = 5, alpha = .6)

    # ax.vlines(0.273548, 0, 1)
    my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax.set_xlim(left = 0, right = 7)
    #ax.set_ylim(bottom = -0.2, top = 3)
    #ax.legend(fontsize = 16, title_fontsize = 20, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/alpha_' + immunization + '.pdf')


    # ax_lam.vlines(0.273548, 0, 1)
    my_plot_layout(ax = ax_lam, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_lam.set_xlim(left = 0, right = 7)
    #ax_lam.set_ylim(bottom = -0.2, top = 3)
    ax_lam.legend(fontsize = 16, title_fontsize = 20, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig_lam.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_' + immunization + '.pdf')

    # ax_t.vlines(0.273548, 0, 1)
    my_plot_layout(ax = ax_t, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_t.set_xlim(left = 0, right = 7)
    #ax_t.set_ylim(bottom = -0.5, top = 6)
    #ax_t.set_yticks(range(1, 12, 3))
    ax_t.legend(fontsize = 16, title_fontsize = 20, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig_t.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_' + immunization + '.pdf')
    plt.close()
    





