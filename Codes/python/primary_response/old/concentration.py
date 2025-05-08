import sys
sys.path.append('../../lib/')
from functions_2 import*
plt.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'


# Define constants and parameters
D = 3*1e-8 * (60*60*24) # Diffusion coefficient [cm^2/days]
k_on = 1e6*24*3600; #(M*days)^-1
b0 = 10**5
u = .000  # Velocity
L = 10 # Length of the domain [cm]
T = 15  # Total simulation time [days]
Nx = 10000  # Number of spatial grid points
Nt = 1500  # Number of time steps
lamA = 6

# Discretize the spatial and temporal dimensions
dx = L / Nx
dt = T / Nt

# Create spatial and temporal grids
xGrid = np.linspace(0, L, Nx + 1)
tGrid = np.linspace(0, T, Nt + 1)

# Initialize concentration matrix
C = np.zeros((Nt + 1, Nx + 1))
C_mask = np.zeros((Nx + 1, 1))
C_mask[0, 0] = 1
# Set initial condition (e.g., a Gaussian profile)
#C[0] = np.exp(-((xGrid - L / 2)**2 / (2 * 1**2)))
growthRate = 6.0  # Exponential growth rate

immunizations = ['vacc', 'infec']

for immunization in immunizations:

    r_array = np.linspace(0, L, Nx+1)
    if immunization == 'vacc':
        C0 = 1e12
        C[0,0] = C0
        for i_t in range(1, Nt):
            C[i_t,:] = (C0/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
        C[C<1]=0
        print('integral = ', np.sum(C[10,xGrid<np.sqrt(2*D*10*dt)]*r_array[xGrid<np.sqrt(2*D*10*dt)]**2*dx)*(4*np.pi)/C0)
        print('integral = ', np.sum(C[100,xGrid<np.sqrt(2*D*100*dt)]*r_array[xGrid<np.sqrt(2*D*100*dt)]**2*dx)*(4*np.pi)/C0)
    else:
        C0 = 1
        C[0,0] = C0
        for i_t in range(1, Nt):
            C[i_t,:] = np.exp(lamA*i_t*dt)*(C0/(4*np.pi*D*i_t*dt)**(3/2))*np.exp(-r_array**2/(4*D*i_t*dt))
        C[C<1]=0
        print('integral = ', np.sum(C[10,xGrid<np.sqrt(2*D*10*dt)]*r_array[xGrid<np.sqrt(2*D*10*dt)]**2*dx)*(4*np.pi)/C0)
        print('integral = ', np.sum(C[100,xGrid<np.sqrt(2*D*100*dt)]*r_array[xGrid<np.sqrt(2*D*100*dt)]**2*dx)/(np.exp(lamA*100*dt))*(4*np.pi))

    x0 = int(0/2)
    r1 = np.sqrt((2*D/lamA)*np.log((lamA*N_A)/(0.8*k_on*b0)))
    rho1 = (lamA*N_A/(k_on*b0))
    print('r1 = %.2f'%r1)
    # Plot the concentration over time
    fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    if immunization == 'vacc':
        points = [x0 + 2*int(r1/dx), x0 + 3*int(r1/dx), x0 + 4*int(r1/dx)]
    else:
        points = [x0 + int(0.5*r1/dx), x0 + 1*int(r1/dx), x0 + 2*int(r1/dx), x0 + 3*int(r1/dx)]#, x0 + 4*int(r1/dx)]
    for i, x in enumerate(points):
        if immunization == 'vacc':
            ax.plot(tGrid[:], C[:,x][:]/10**0, label = r'$%.1f$'%((xGrid[x])/r1), lw = 5, color = antigen_color)
            lamAeff = 4*D*(np.log(C0/(rho1))**2)/(xGrid[x])**2# - 4*D*(np.log(C0/(rho1)))/(2*(xGrid[x]-5.0)**2)
            #t1 = xGrid[x]**2 / (4*D*np.log(C0))
            #t1 = tGrid[C[:,x]<((lamAeff*N_A/(k_on*b0)))][-int(30*(.3*i+1))]
            #exp_growth = np.exp(lamAeff*(tGrid - t1))
            print('lam_A^eff = %.1f'%lamAeff)
            #ax.plot(tGrid[exp_growth<(rho1)], rho1*exp_growth[exp_growth<(rho1)], lw = 2, color = ax.get_lines()[-1].get_color(), ls = 'dashed')
            #ax.hlines(rho1, 1, 7, lw = 2, color = antigen_color, ls = 'dotted')
            print('%.1e'%rho1)
        else:
            ax.plot(tGrid[:], C[:,x][:]/10**0, label = r'$%.1f$'%((xGrid[x])/r1), lw = 5, color = antigen_color, alpha = (4*int(r1/dx)-x)/(4*int(r1/dx)-int(0.5*r1/dx)))
            #t1 = xGrid[x]**2 / (4*D*np.log(C0))
            #t1 = tGrid[C[:,x]<1][-int(5*(1.5*i+1))]
            #exp_growth = np.exp(lamA*(tGrid - t1))
            print('lam_A = %.1f'%lamA)
            #ax.plot(tGrid[exp_growth<1e8], exp_growth[exp_growth<1e8], lw = 2, color = ax.get_lines()[-1].get_color(), ls = 'dashed')

    if immunization == 'vacc':
        # ax.plot(tGrid[:], 1e-0*np.exp(growthRate * 1 * tGrid), color = 'k', ls = '--')
        ax.set_xlim(left = 0, right = 5)
        ax.set_ylim(bottom = 2e3, top = 2e12)
    else:
        # ax.plot(tGrid[:], 1e-4*np.exp(growthRate * tGrid), color = 'k', ls = '--')
        ax.set_xlim(left = 0, right = 7)
        ax.set_ylim(bottom = 1e1, top = 2e11)
    
    my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )  
    ax.legend(fontsize = 22, title_fontsize = 24, title = r"$r/r_1$", loc = 4)
    fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/antigen_' + immunization + '_log2.pdf')

    my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )  
    ax.legend(fontsize = 22, title_fontsize = 24, title = r"$r/r_1 \hspace{2mm} (cm)$", loc = 0)
    fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/antigen_' + immunization + '2.pdf')
    plt.close()

    # Plot the concentration over time at x0
    fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
    ax.plot(tGrid[:], C[:,x0][:], label = r'$%.2f$'%(xGrid[x0]), lw = 2)
    # plt.plot(tGrid[:], 1e-40*np.exp((u**3*(tGrid[:])**2)/(xGrid[x]*D))/np.sqrt(2*np.pi*xGrid[x]/u) )
    #plt.plot(tGrid[:], 1e-10*np.exp(growthRate*tGrid[:]))
    if immunization == 'vacc':
         # ax.plot(tGrid[:], 1e-0*np.exp(growthRate * 1 * tGrid), color = 'k', ls = '--')
        ax.set_xlim(left = 0, right = 4)
        ax.set_ylim(bottom = C0/1000, top = C0)
    else:
        #t1 = xGrid[x]**2 / (4*D*np.log(C0))
        t1 = 0
        exp_growth = np.exp(lamA*(tGrid - t1))
        print(t1, lamA)
        ax.plot(tGrid[exp_growth<1e4], exp_growth[exp_growth<1e4], lw = 2, color = ax.get_lines()[-1].get_color(), ls = 'dashed')
        # ax.plot(tGrid[:], 1e-4*np.exp(growthRate * tGrid), color = 'k', ls = '--')
        ax.set_xlim(left = 0, right = 9)
        ax.set_ylim(bottom = 1e0, top = 1e12)

    my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )  
    ax.legend(fontsize = 22, title_fontsize = 24, title = r'$x \hspace{2mm} (cm)$', loc = 0)
    fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/antigen_' + immunization + '_02.pdf')
    plt.close()


    # Plot the concentration in space
    fig, ax = plt.subplots(figsize=(5*1.62,5), gridspec_kw={'left':0.12, 'right':.98, 'bottom':.15, 'top': 0.94})
    for t in [50, 100, 150, 200, 250, 300, 350, 400]:
        if t%100==0:
            ax.plot((xGrid[:])/r1, C[t,:]/10**0, label = r'$%d$'%(tGrid[t]), lw = 5, color = antigen_color, alpha = t/(1000) + t/(100)*.1)
        else:
            ax.plot((xGrid[:])/r1, C[t,:]/10**0, lw = 5, color = antigen_color, alpha = t/(1000) + t/(100)*.1)
        # plt.plot(tGrid[:], 1e-40*np.exp((u**3*(tGrid[:])**2)/(xGrid[x]*D))/np.sqrt(2*np.pi*xGrid[x]/u) )
    #plt.plot(tGrid[:], 1e-10*np.exp(growthRate*tGrid[:]))
    if immunization == 'vacc':
        ax.set_xlim(left = 0, right = 4)
        ax.set_ylim(bottom = 2e3, top = 2e17)
        # ax.plot(tGrid[:], 1e-0*np.exp(growthRate * 1 * tGrid), color = 'k', ls = '--')
    else:
        # ax.plot(tGrid[:], 1e-4*np.exp(growthRate * tGrid), color = 'k', ls = '--')
        ax.set_xlim(left = 0, right = 4)
        ax.set_ylim(bottom = 1e1, top = 2e11)

    my_plot_layout(ax = ax, xscale='linear', yscale= 'log', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )  
    ax.legend(fontsize = 22, title_fontsize = 24, title = r'$t \hspace{2mm} (\mathrm{days})$', loc = 0)
    fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/antigen_profile_' + immunization + '_log2.pdf')
    plt.close()

    my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 30, x_fontsize=30, y_fontsize=30 )  
    ax.legend(fontsize = 22, title_fontsize = 24, title = r'$t \hspace{2mm} (days)$', loc = 0)
    fig.savefig('../../../Figures/primary_immune_response/1_Dynamics/antigen_concentration/antigen_profile_' + immunization + '2.pdf')
    plt.close()


