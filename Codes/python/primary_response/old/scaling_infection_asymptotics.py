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
Nt = 40000  # Number of time steps
lamA = 6

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
Ds = np.logspace(np.log10(0.03*D0), np.log10(3000*D0), 10)
Nstars = [1e10/100, 1e10/1]
r0s = np.linspace(0.1, 1, 3)

for immunization in immunizations:
    print(immunization)    
    results = defaultdict(list)
    for r0 in r0s:
        for i_D, D in enumerate(tqdm(Ds)):
            Nstar = 1e10
            r_array = np.linspace(0, L, Nx+1)
            growth_rate = lamA
            integral_in = np.zeros_like(tGrid)
            for i_t in range(1, Nt+1):
                C[i_t,:] = np.exp(lamA*i_t*dt)*(1/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
                #C[C<1]=0
                integral_in[i_t] = np.sum(C[i_t,xGrid>r0]*r_array[xGrid>r0]**2*dx)*(4*np.pi)

            popt, pcov = curve_fit(f = my_linear_func, xdata = tGrid[integral_in>Nstar/2][:10], ydata = np.log(integral_in)[integral_in>Nstar/2][:10])
            tstar = tGrid[integral_in>Nstar][0]
            rstar = np.sqrt(2*D*tstar)
            results['N'].append(Nstar)
            results['D'].append(D)
            results['lamA'].append(popt[1])
            results['r0'].append(r0)
            results['t'].append(tstar)
            results['alpha'].append(rstar/r0)
            #print(popt[1])

    results_db = pd.DataFrame(results)
    # print(results_db)
    fig_lam, ax_lam = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig_t, ax_t = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    alpha_array = np.logspace(np.log10(0.04), np.log10(30), 50)
    
    markers = ['o', 's', '^', 'D', '*']
    for i_r0, r0 in enumerate(r0s):
        data_N = results_db.loc[results_db['r0']==r0]
        alpha_D_r_array = np.logspace(np.min(r0**(1/2)*np.log10(data_N['alpha']/data_N['D']**(1/4))), np.log10(4), 100)
       
        ax_lam.plot(data_N['alpha'], (data_N['lamA']-lamA)*r0**2/data_N['D'], lw = 5, ls = ''
            , marker = markers[i_r0], ms = 6, color = my_purple, label = r'$%.1f$'%(r0))
        
        ax_t.plot(data_N['alpha'], (data_N['t']-((1/lamA)*np.log(Nstar)))*lamA, lw = 5, ls = ''
            , marker = markers[i_r0], ms = 6, color = my_purple, label = r'$%.1f$'%(r0))
    

    ax_lam.plot(alpha_array, .5*alpha_array**(-4), color = my_purple2, lw = 2)
    ax_lam.plot(alpha_array, .5*alpha_array**(-5), color = my_purple2, lw = 2)
    my_plot_layout(ax = ax_lam, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_lam.set_xlim(left = 0, right = 7)
    #ax_lam.set_ylim(bottom = 1e0, top = 2e11)
    ax_lam.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in C0s])
    fig_lam.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/lam_' + immunization + '_asymptotics.pdf')


    ax_t.plot(alpha_array, 0.5*alpha_array**(-2), color = my_purple2, lw = 2)
    ax_t.plot(alpha_array, 0.5*alpha_array**(-3), color = my_purple2, lw = 2)
    my_plot_layout(ax = ax_t, xscale='log', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    #ax_t.set_xlim(left = 0, right = 7)
    #ax_t.set_ylim(bottom = 1e0, top = 2e11)
    ax_t.legend(fontsize = 20, title_fontsize = 24, loc = 0, title = r'$r_0$')#, labels = [r'$%.1e$'%i for i in C0s])
    fig_t.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/t_' + immunization + '_asymptotics.pdf')
    plt.close()

  



  