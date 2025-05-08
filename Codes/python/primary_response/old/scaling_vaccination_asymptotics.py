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
Nt = 4000  # Number of time steps
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
colorsD = [my_purple, my_blue, my_green, my_red]
D0 = 1e-3
Ds = np.logspace(np.log10(0.1*D0), np.log10(10*D0), 3)
r0s = np.linspace(0.01, .2, 6)
N0s = np.logspace(13., 16., 6)
N0 = 1e13

for immunization in immunizations:
    print(immunization)
    
    results = defaultdict(list)
    for i_D, D in enumerate(Ds):
        for N0 in tqdm(N0s):
            Nstar = 1e10
            r_array = np.linspace(0, L, Nx+1)
            integral_in = np.zeros_like(tGrid)
            for i_t in range(1, Nt+1):
                C[i_t,:] = (N0/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
                #C[C<1]=0
                integral_in[i_t] = np.sum(C[i_t,xGrid>r0]*r_array[xGrid>r0]**2*dx)*(4*np.pi)
            
            popt, pcov = curve_fit(f = my_linear_func, xdata = tGrid[integral_in>Nstar/2][:5], ydata = np.log(integral_in)[integral_in>Nstar/2][:5])
            tstar = tGrid[integral_in>(Nstar)][0]
            rstar = np.sqrt(2*D*tstar)
            results['N0'].append(N0)
            results['D'].append(D)
            results['lamA'].append(popt[1])
            results['r0'].append(r0)
            results['t'].append(tstar)
            results['alpha'].append(rstar/r0)
            #print(popt[1])

    results_db = pd.DataFrame(results)
    fig_lam, ax_lam = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_t, ax_t = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    alpha_array = np.logspace(np.log10(0.03), np.log10(.07), 50)
    markers = ['o', 's', '^', 'D', '*']
    for i_D, D in enumerate(Ds):
        data_D = results_db.loc[results_db['D']==D]
        ax_lam.plot(data_D['alpha']**2, (data_D['lamA'])*r0**2/data_D['D'], lw = 5, ls = ''
            , marker = markers[i_D], ms = 6, color = my_purple, label = r'$%.0e$'%(D))
        
        ax_t.plot(data_D['alpha']**2, 2*data_D['t']*data_D['D']/r0**2, lw = 5, ls = ''
            , marker = markers[i_D], ms = 6, color = my_purple, label = r'$%.0e$'%(D))
    
    ax_lam.plot(alpha_array, 1.1*alpha_array**(-2), color = my_purple2, lw = 2)
    #ax_lam.plot(alpha_array, .5*alpha_array**(-5), color = my_purple2, lw = 2)
    my_plot_layout(ax = ax_lam, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_lam.set_xlim(left = 0, right = 7)
    #ax_lam.set_ylim(bottom = 1e0, top = 2e11)
    ax_lam.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$D$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig_lam.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_' + immunization + '_asymptotics.pdf')

    ax_t.plot(alpha_array, alpha_array, color = my_purple2, lw = 2)
    # ax_t.plot(alpha_array, 0.05*alpha_array**(-3), color = my_purple2, lw = 2)
    my_plot_layout(ax = ax_t, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_t.set_xlim(left = 0, right = 7)
    #ax_t.set_ylim(bottom = 1e0, top = 2e11)
    ax_t.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$D$')#, labels = [r'$%.1e$'%i for i in N0s])
    fig_t.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_' + immunization + '_asymptotics.pdf')
    plt.close()
    





