import sys
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'


# Define constants and parameters
D0 = 3*1e-8 * (60*60*24) # Diffusion coefficient [cm^2/days]
k_on = 1e6; #(M*days)^-1
Kstep = (0.5/60)/k_on;
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
r0 = 0.5 # cm

# Discretize the spatial and temporal dimensions
dx = L / Nx
dt = T / Nt

# Create spatial and temporal grids
xGrid = np.linspace(0, L, Nx + 1)
tGrid = np.linspace(0, T, Nt + 1)

# Initialize concentration matrix
C = np.zeros((Nt + 1, Nx + 1))

immunizations = ['infec']
colorsD = [my_purple, my_blue, my_green, my_red]
D0 = 3e-3
Ds = [0.2*D0, D0, 5*D0]
Nstars = [1e10/100, 1e10/1]
r0s = [0.5]

for immunization in immunizations:
    print(immunization)
    fig1, ax1 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    fig2, ax2 = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    ax2.plot(tGrid[:], np.exp(lamA*tGrid[:]), lw = 5, color = antigen_color, label = r'$0$')
    D0 = 3e-4
    for r0 in r0s:
        for i_D, D in enumerate(Ds):
            Nstar = ((1/N_A)*k_on*b0*(Kstar/Kstep)**(-p)/(lamA/(60*60*24)))**-1
            r_array = np.linspace(0, L, Nx+1)
            
            N0 = 1
            growth_rate = lamA
            #C[0,0] = N0
            NAeff = np.zeros_like(tGrid)
            integral_out = np.zeros_like(tGrid)
            for i_t in range(1, Nt):
                C[i_t,:] = np.exp(lamA*i_t*dt)*(N0/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
                #C[C<1]=0
                NAeff[i_t] = np.sum(C[i_t,xGrid>r0]*r_array[xGrid>r0]**2*dx)*(4*np.pi)
                # integral_out[i_t] = np.sum(C[i_t,xGrid>np.sqrt(2*D*i_t*dt)]*r_array[xGrid>np.sqrt(2*D*i_t*dt)]**2*dx)*(4*np.pi)
            lamAeff_array = np.diff(np.log(NAeff))/np.diff(tGrid)
            R = 1-np.exp(-np.cumsum((1/N_A)*k_on*(60*60*24)*b0*(Kstar/Kstep)**(-p)*NAeff*dt))
            tstar = tGrid[R<(0.5)][-1]
            lamAeff = lamAeff_array[R[:-1]<0.5][-1]
            rstar = np.sqrt(D*tstar)
            alpha = r0/rstar
            print('D=%.0e'%D, 'alpha = %.2f'%alpha)

            if i_D == 1:
                for i_t in [100, 200, 300, 400, 500]:
                    ax1.plot(xGrid/r0, C[i_t, :], lw = 5, color = colorsD[i_D], alpha = i_t/(500), label = r'$%.1f$'%(i_t*dt))
                    if i_D == 1:
                        ax1.fill_between((xGrid/r0)[xGrid>r0], np.ones_like(C[i_t, xGrid>r0]), C[i_t, xGrid>r0], lw = 5, color = colorsD[i_D], alpha = .1)
                        ax1.hlines(N0*0.08*np.exp(i_t*dt*growth_rate)/(2*np.pi*D*dt*i_t)**(3/2), 0, 1.5*np.sqrt(2*D*i_t*dt)/r0, color = 'black', ls = '--')
            
            ax2.plot(tGrid[:], NAeff, lw = 4, ls = '-', label = r'$%.1f$'%(alpha), color = colorsD[i_D])
            # ax2.vlines(tstar, 1, 1e12, lw = 1, ls = '--', color = 'k')

    # Plot the concentration over space
    my_plot_layout(ax = ax1, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax1.set_xlim(left = 0, right = 4)
    ax1.set_ylim(bottom = 1e0, top = np.max([1e12, N0*100]))
    ax1.legend(fontsize = 20, title_fontsize = 24, title = r'$t \mathrm{ [days]}$', loc = 1)
    fig1.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/integral_space_' + immunization + '_log.pdf')

    # Plot the concentration over time
    my_plot_layout(ax = ax2, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )
    ax2.set_xlim(left = 0, right = 8)
    ax2.set_ylim(bottom = 1e0, top = 1e12)
    ax2.legend(fontsize = 20, title_fontsize = 24, title = r'$\alpha$', loc = 4)
    fig2.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/integral_time_' + immunization + '_log.pdf')

    
        





