import sys
sys.path.append('../../lib/')
from funcs import*
plt.rcParams['text.usetex'] = True

def simulate_poisson_gaussian(g, mu, N_mut, T, DGs0):
    event_times = [0]

    DDGs = np.zeros((g, N_mut))
    # DDGs[:, 0] = DGs0

    # Initialize time to zero
    t = 0
    my_pos = np.random.randint(1*int(N_mut/5), 4*int(N_mut/5))
    for mutation in range(N_mut-1):

        # Draw next event times for all processes
        dt = np.random.exponential(1/(mu*g))
        t_next = t + dt
        # if t_next > T:
        #     break
        event_times.append(t_next)
        
        # Identify which process had the event
        if mutation == my_pos:
            DDGs[:, mutation+1] = DDGs[:, mutation]
            DDGs[0, mutation+1]+= 1
        else:
            r = np.random.random()
            array_e = np.linspace(0, 1, g+1)
            e_mut = np.searchsorted(array_e, r, side = 'left')
                
            DDGs[:, mutation+1] = DDGs[:, mutation]
            DDGs[e_mut-1, mutation+1]+= abs(np.random.normal(1, 0.5))
        
        t = t_next  # Update time
        
    return event_times, DDGs, my_pos


project = 'memory_response'
subproject = 'multi-epitope'
experiment = 0
output_plot = '../../../Figures/'+project+'/'+subproject+'/'+str(experiment)
os.makedirs(output_plot, exist_ok=True)

# Parameters
n_ens = 10000
gs = [1, 5]  # Number of Poisson processes
mu = 1.0  # Poisson rate
T = 12  # Total simulation time
theta = 1.8  # Values of theta to compare
gamma = 0.6
my_colors = [my_red, my_cyan, my_green, my_green, my_brown]
my_colors2 = ['darkred', 'lightcoral', 'darkorange', 'crimson', 'peru']
alpha = 1e-10
depth = 6

print(theta*gamma)
for i_g, g in enumerate(tqdm(gs)):
    N_mut = T*g
    fig, ax =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figtimes, axtimes =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figTittau, axTittau =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figXmut, axXmut =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figXtau, axXtau =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figXTit, axXTit =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figX, axX =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figHtau, axHtau =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figHTit, axHTit =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})
    figHX, axHX =  plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.9, 'bottom':.1, 'top': 0.9})

    times = np.zeros(N_mut)
    Tdrop = np.zeros(N_mut)
    Xinvitro = np.zeros(N_mut)
    Xinvivo = np.zeros(N_mut)
    NA_peak = np.zeros(N_mut)
    Hinvitro = np.zeros(N_mut)
    Hinvivo = np.zeros(N_mut)

    for i_ens in range(n_ens):
        
        # Initialize epsilon
        DGs0 = np.random.normal(1.5, .5, g)  # Gaussian distributed random numbers
        DGs0 = np.sort(DGs0) - np.min(DGs0)  # Order and shift so first element is zero
        DGs0 = np.zeros(g)

        # Run the simulation
        event_times, DDGs, my_pos = simulate_poisson_gaussian(g, mu, N_mut, T, DGs0)
        DDGsDGs0 = DDGs - DGs0[:, np.newaxis]
        gammaDDGsDGs0 = gamma*DDGs - DGs0[:, np.newaxis]
        
        # IN VITRO RESPONSE
        Xinvitro_i = np.sum(np.exp(-DDGsDGs0)/np.sum(np.exp(-DGs0)), axis = 0)
        Hinvitro_i = 1 - np.product(1-(1/(1+(alpha*1e3*1e8*np.exp(-DDGsDGs0)/np.sum(np.exp(-DGs0)))**-1)), axis = 0)
        # Hinvitro_i = 1/(1+(alpha*1e3*1e9*Xinvitro_i)**-1)
        
        # To be use as biophysical relevant axis.
        Tdrop_i = Xinvitro_i/Xinvitro_i[0]

        # IN VIVO RESPONSE
        Xinvivo_i = np.sum(np.exp(-theta*gammaDDGsDGs0)/np.sum(np.exp(-theta*DGs0)), axis = 0)
        NA_peak_i = np.exp(6*2)*(alpha*1e3*1e7*Xinvivo_i)**(-3)
        Hinvivo_i = 1/(1+NA_peak_i/1e12)
        
        times+=event_times
        Tdrop+=np.log(Tdrop_i)
        Xinvitro+=np.log(Xinvitro_i)
        Xinvivo+=np.log(Xinvivo_i)
        NA_peak+=np.log(NA_peak_i)
        Hinvitro+=np.log(Hinvitro_i)
        Hinvivo+=np.log(Hinvivo_i)

        if i_ens%200==0: 

            axtimes.step(event_times, np.array(range(N_mut))/g, color = 'k', alpha = .5, lw = .5)

            axTittau.step(event_times, -np.log(Tdrop_i), color = my_blue2, alpha = .5, lw = .5)
            
            # axXmut.step(event_times, np.log(Xinvitro_i), color = my_blue2, alpha = .5, lw = .5)
            axXmut.step(range(N_mut), np.log(Xinvivo_i) + depth, color = my_green, alpha = .5, lw = .5)
            axXmut.vlines(range(N_mut)[my_pos], np.log(Xinvivo_i[my_pos+1])+depth, np.log(Xinvivo_i[my_pos])+depth, color = my_green, zorder = 20)
            axXmut.scatter(range(N_mut)[my_pos], np.log(Xinvivo_i[my_pos])+depth, color = my_green, s = 5, marker = 'D', zorder = 20)

            axXtau.step(event_times, np.log(Xinvitro_i), color = my_red, alpha = .5, lw = .5)
            # axXtau.step(event_times, np.log(Xinvivo_i) + depth, color = my_green, alpha = .5, lw = .5)
            # axXtau.vlines(event_times[my_pos], np.log(Xinvivo_i[my_pos+1])+depth, np.log(Xinvivo_i[my_pos])+depth, color = my_green, zorder = 20)
            # axXtau.scatter(event_times[my_pos], np.log(Xinvivo_i[my_pos])+depth, color = my_green, s = 5, marker = 'D', zorder = 20)
            
            axXTit.step(-np.log(Tdrop_i), np.log(Xinvitro_i) + depth, color = my_blue2, alpha = .5, lw = .5)
            axXTit.step(-np.log(Tdrop_i), np.log(Xinvivo_i)  + depth , color = my_green, alpha = .5, lw = .5)
            axXTit.vlines(-np.log(Tdrop_i)[my_pos], np.log(Xinvivo_i[my_pos+1])+depth, np.log(Xinvivo_i[my_pos])+depth, color = my_green, zorder = 20)
            axXTit.scatter(-np.log(Tdrop_i)[my_pos], np.log(Xinvivo_i[my_pos])+depth, color = my_green, s = 5, marker = 'D', zorder = 20)

            axX.step(-np.log(Xinvitro_i), -np.log(Xinvitro_i), color = my_blue2, alpha = .5, lw = .5)
            axX.step(-np.log(Xinvitro_i), -np.log(Xinvivo_i), color = my_green, alpha = .5, lw = .5)

            # axHtau.step(event_times, Hinvitro_i, color = my_blue2, alpha = .5, lw = .5)
            # axHtau.vlines(event_times[my_pos], Hinvitro_i[my_pos], Hinvitro_i[my_pos+1], color = my_blue2, zorder = 20)
            # axHtau.scatter(event_times[my_pos], Hinvitro_i[my_pos], color = my_blue2, marker = '*', zorder = 20)
            
            axHtau.step(event_times, 1 - Hinvivo_i, color = my_purple, alpha = .2, lw = .4)        
            axHtau.vlines(event_times[my_pos], 1 - Hinvivo_i[my_pos], 1 - Hinvivo_i[my_pos+1], color = my_purple, zorder = 20)
            axHtau.scatter(event_times[my_pos], 1 - Hinvivo_i[my_pos], color = my_purple, s = 5, marker = 'D', zorder = 20)

            axHX.step(-np.log(Xinvitro_i), Hinvitro_i, color = my_blue2, alpha = .5, lw = .5)
            axHX.step(-np.log(Xinvivo_i), Hinvivo_i, color = my_green, alpha = .5, lw = .5)

            axHTit.step(-np.log(Tdrop_i), Hinvitro_i, color = my_blue2, alpha = .5, lw = .5)
            axHTit.step(-np.log(Tdrop_i), Hinvivo_i, color = my_green, alpha = .5, lw = .5)

        if g==1:
            if i_ens in [1, 2, 3, 4, 5]:
                ax.step(event_times, DDGs[0,:], color = 'grey')
                # ax.vlines(event_times[my_pos], DDGs[0,:][my_pos], DDGs[0,:][my_pos+1], color = 'k', zorder = 20)
                # ax.scatter(event_times[my_pos], DDGs[0,:][my_pos], color = 'k', marker = '*', zorder = 20)
        else:
            if i_ens == 0:
                for e in range(g):
                    ax.step(event_times, DDGs[e,:], color = my_colors2[e])
                    # ax.vlines(event_times[my_pos], DDGs[0,:][my_pos], DDGs[0,:][my_pos+1], color = 'k', zorder = 20)
                    # ax.scatter(event_times[my_pos], DDGs[0,:][my_pos], color = 'k', marker = '*', zorder = 20)

    times = times/n_ens
    Tdrop = np.exp(Tdrop/n_ens)
    Xinvitro = np.exp(Xinvitro/n_ens)
    Xinvivo = np.exp(Xinvivo/n_ens)
    NA_peak = np.exp(NA_peak/n_ens)
    Hinvitro = np.exp(Hinvitro/n_ens)
    Hinvivo = np.exp(Hinvivo/n_ens)

    axTittau.plot(times, -np.log(Tdrop), color = my_blue2, alpha = 1, label = r'$\textit{in vitro}$')

    axXmut.plot(range(N_mut), np.log(Xinvitro) + 6, color = my_blue2, alpha = 1, label = r'$\textit{in vitro}$')
    # axXmut.plot(range(N_mut), np.log(Xinvivo) + depth, color = my_green, alpha = 1, label = r'$\textit{in vivo}$')
    axXmut.plot(range(N_mut), -0.25*np.array(range(N_mut)) + depth, color = my_blue2, alpha = 1, ls = ':')
    axXmut.plot(range(N_mut), -np.array(range(N_mut)) + depth, color = my_blue2, alpha = 1, ls = '--')
    # axXmut.plot(times, -((theta*gamma))*times + 6, color = my_green, alpha = 1, ls = '--')
    axXmut.hlines(0, 0, (T - 1.5)*g, ls = '--', color = 'k')
    
    axXtau.plot(times, np.log(Xinvitro), color = my_red, lw = 2, alpha = 1, label = r'$%d$'%g)
    # axXtau.plot(times, 0.3*np.log(Xinvivo) + depth, color = my_green, alpha = 1, label = r'$\textit{in vivo}$')
    axXtau.plot(times, -times, color = 'grey', alpha = 1, lw = 2, ls = '--', label = r'$1$')
    # axXtau.plot(times, 0.8*times, color = my_blue2, alpha = 1, ls = ':')
    # axXtau.plot(times, -((theta*gamma))*times + 6, color = my_green, alpha = 1, ls = '--')
    axXtau.hlines(0, 0, T - 1.5, ls = '--', color = 'k')

    # axXTit.plot(-np.log(Tdrop), np.log(Xinvitro) + depth, color = my_blue2, alpha = 1, label = r'$\textit{in vitro}$')
    axXTit.plot(-np.log(Tdrop), np.log(Xinvivo) + depth, color = my_green, alpha = 1, label = r'$\textit{in vivo}$')
    axXTit.hlines(0, 0, T - 1.5, ls = '--', color = 'k')

    axX.plot(-np.log(Xinvitro), -np.log(Xinvitro), color = my_blue2, alpha = 1, label = r'$\textit{in vitro}$')
    axX.plot(-np.log(Xinvitro), -np.log(Xinvivo), color = my_green, alpha = 1, label = r'$\textit{in vivo}$')
    axX.plot(-np.log(Xinvitro), -((theta*gamma))*np.log(Xinvitro), color = my_green, alpha = 1, ls = '--')

    axHtau.plot(times, 1 - Hinvitro, color = my_blue2, alpha = 1, label = r'$\textit{in vitro}$', ls = '--')
    axHtau.plot(times, 1 - Hinvivo, color = my_purple, alpha = 1, label = r'$\textit{in vivo}$')
    axHtau.hlines(0.5, 0, T - 1.5, ls = '--', color = 'k')

    axHX.plot(-np.log(Xinvitro), Hinvitro, color = my_blue2, alpha = 1, label = r'$\textit{in vitro}$')
    axHX.plot(-np.log(Xinvivo), Hinvivo, color = my_green, alpha = 1, label = r'$\textit{in vivo}$')
    axHX.hlines(0.5, 0, T - 1.5, ls = '--', color = 'k')

    axHTit.plot(-np.log(Tdrop), Hinvitro, color = my_blue2, alpha = 1, label = r'$\textit{in vitro}$')
    axHTit.plot(-np.log(Tdrop), Hinvivo, color = my_green, alpha = 1, label = r'$\textit{in vivo}$')
    axHTit.hlines(0.5, 0, T - 1.5, ls = '--', color = 'k')

    my_plot_layout(ax = axtimes, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axtimes.set_xlabel('Time')
    # axtimes.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axtimes.set_ylim(bottom = -0.4, top = T-1)
    axtimes.set_xlim(left = -0.4, right = T-1)
    # axtimes.set_title(r'Trajectory of $Z$ over Time')
    # axtimes.legend(fontsize = 16)
    figtimes.savefig(output_plot + '/times_g-%d.pdf'%(g))

    my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # ax.set_xlabel('Time')
    # ax.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    ax.set_ylim(bottom = -0.4, top = T-1)
    ax.set_xlim(left = -0.4, right = T-1)
    # ax.set_title(r'Trajectory of $Z$ over Time')
    # ax.legend(fontsize = 16)
    fig.savefig(output_plot + '/DDG_g-%d.pdf'%(g))

    my_plot_layout(ax = axTittau, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axTittau.set_xlabel('Time')
    # axTittau.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axTittau.set_ylim(bottom = 1e-2, top = T-1)
    axTittau.set_xlim(left = -0.4, right = T-1)
    # axTittau.set_title(r'Trajectory of $Z$ over Time')
    axTittau.legend(fontsize = 16)
    figTittau.savefig(output_plot + '/Tit_tau_g-%d.pdf'%(g))

    my_plot_layout(ax = axXmut, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axXmut.set_xlabel('Time')
    # axXmut.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axXmut.set_ylim(bottom = -8 + depth + 0.2, top = 0 + depth + 0.2)
    axXmut.set_xlim(left = -0.4, right = (T-1)*g)
    # axXmut.set_xticks([])
    # axXmut.set_title(r'Trajectory of $Z$ over Time')
    axXmut.legend(fontsize = 16)
    figXmut.savefig(output_plot + '/X_mut_g-%d.pdf'%(g))

    my_plot_layout(ax = axXtau, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axXtau.set_xlabel('Time')
    # axXtau.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    # axXtau.set_ylim(bottom = -8 + depth + 0.2, top = 0 + depth + 0.2)
    axXtau.set_xlim(left = -0.4, right = T-1)
    # axXtau.set_xticks([])
    # axXtau.set_title(r'Trajectory of $Z$ over Time')
    axXtau.legend(fontsize = 18, title = r'$g$', title_fontsize = 20)
    figXtau.savefig(output_plot + '/X_tau_g-%d.pdf'%(g))

    my_plot_layout(ax = axXTit, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axXTit.set_xlabel('Time')
    # axXTit.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axXTit.set_ylim(bottom = -8 + depth + 0.2, top = 0 + depth + 0.2)
    axXTit.set_xlim(left = -0.4, right = T-1)
    axXTit.set_xticks([])
    # axXTit.set_title(r'Trajectory of $Z$ over Time')
    axXTit.legend(fontsize = 16)
    figXTit.savefig(output_plot + '/X_Tit_g-%d.pdf'%(g))

    my_plot_layout(ax = axX, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axX.set_xlabel('Time')
    # axX.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axX.set_ylim(bottom = -0.4, top = T-1)
    axX.set_xlim(left = -0.4, right = T-1)
    # axX.set_title(r'Trajectory of $Z$ over Time')
    axX.legend(fontsize = 16)
    figX.savefig(output_plot + '/X_X_g-%d.pdf'%(g))

    my_plot_layout(ax = axHtau, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axHtau.set_xlabel('Time')
    # axHtau.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axHtau.set_ylim(top = 1.05, bottom = 1e-2)
    axHtau.set_xlim(left = -0.4, right = T-1)
    axHtau.set_xticks([])
    axHtau.set_yticks([1, 0.5, 0])
    # axHtau.set_title(r'Trajectory of $Z$ over Time')
    axHtau.legend(fontsize = 18)
    figHtau.savefig(output_plot + '/H_tau_g-%d.pdf'%(g))

    axHtau.plot(times + 4 + g - 2, np.exp(-1*times), color = my_blue2, alpha = 1, ls = '--')
    axHtau.plot(times + 7 + g - 2, np.exp(-3*theta*gamma*times), color = my_green, alpha = 1, ls = '--')

    my_plot_layout(ax = axHtau, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axHtau.set_xlabel('Time')
    # axHtau.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axHtau.set_ylim(top = 1.05, bottom = 1e-2)
    axHtau.set_xlim(left = -0.4, right = T-1)
    # axHtau.set_title(r'Trajectory of $Z$ over Time')
    axHtau.legend(fontsize = 16)
    figHtau.savefig(output_plot + '/H_tau_g-%d_log.pdf'%(g))

    my_plot_layout(ax = axHX, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axHX.set_xlabel('Time')
    # axHX.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axHX.set_ylim(top = 1.05, bottom = 1e-2)
    axHX.set_xlim(left = -0.4, right = T-1)
    # axHX.set_title(r'Trajectory of $Z$ over Time')
    axHX.legend(fontsize = 16)
    figHX.savefig(output_plot + '/H_X_g-%d.pdf'%(g))

    axHX.plot(-np.log(Xinvitro) + 4, Xinvitro**(1), color = my_blue2, alpha = 1, ls = '--')
    axHX.plot(-np.log(Xinvivo) + 4, Xinvivo**(3), color = my_green, alpha = 1, ls = '--')

    my_plot_layout(ax = axHX, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axHX.set_xlabel('Time')
    # axHX.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axHX.set_ylim(top = 1.05, bottom = 1e-2)
    axHX.set_xlim(left = -0.4, right = T-1)
    # axHX.set_title(r'Trajectory of $Z$ over Time')
    axHX.legend(fontsize = 16)
    figHX.savefig(output_plot + '/H_X_g-%d_log.pdf'%(g))

    my_plot_layout(ax = axHTit, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axHTit.set_xlabel('Time')
    # axHTit.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axHTit.set_ylim(top = 1.05, bottom = 1e-2)
    axHTit.set_xlim(left = -0.4, right = T-1)
    # axHTit.set_title(r'Trajectory of $Z$ over Time')
    axHTit.legend(fontsize = 16)
    figHTit.savefig(output_plot + '/H_Tit_g-%d.pdf'%(g))

    axHTit.plot(-np.log(Tdrop) + 4, Tdrop**(1), color = my_blue2, alpha = 1, ls = '--')
    axHTit.plot(-np.log(Tdrop) + 7, Tdrop**(3*theta*gamma), color = my_green, alpha = 1, ls = '--')

    my_plot_layout(ax = axHTit, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
    # axHTit.set_xlabel('Time')
    # axHTit.set_ylabel(r'$Z = \sum_{i=1}^{n}(\exp(\theta * \epsilon_i))$')
    axHTit.set_ylim(top = 1.2, bottom = 1e-2)
    axHTit.set_xlim(left = -0.4, right = T-1)
    # axHTit.set_title(r'Trajectory of $Z$ over Time')
    axHTit.legend(fontsize = 16)
    figHTit.savefig(output_plot + '/H_Tit_g-%d_log.pdf'%(g))