import sys
sys.path.append('../../lib/')
from funcs import*
import random

def gillespie_algorithm(initial_state, reactions, propensities, parameters, max_time):
    """
    Implements the Gillespie Stochastic Simulation Algorithm.
    
    """
    
    # Initialize
    state = np.array(initial_state, dtype=int)
    time = 0.0
    times = [time]
    states = [state.copy()]
    
    while time < max_time:
        
        # Compute reaction propensities
        a = np.array([prop(state, parameters) for prop in propensities])
        a0 = np.sum(a)
        
        if a0 == 0:
            break  # No reactions can occur
        
        # Generate two random numbers
        r1, r2 = random.random(), random.random()
        
        # Compute time step
        tau = -np.log(r1) / a0
        
        # Select reaction to fire
        cumulative_a = np.cumsum(a)
        reaction_index = np.searchsorted(cumulative_a, r2 * a0)
        
        # Update system state
        state = reactions[reaction_index](state)
        time += tau
        
        # Record state and time
        times.append(time)
        states.append(state.copy())
    
    return times, states

# Example usage
if __name__ == "__main__":
    # State: p, pMHC
    
    # Define reactions:
    
    def reaction1(state): # Import of peptide
        new_state = state.copy()
        new_state[0] += 1  # p increases
        return new_state
    
    def reaction2(state): # Export of peptide
        new_state = state.copy()
        new_state[0] -= 1  # p decreases
        return new_state

    def reaction3(state): # p-MHC binding
        new_state = state.copy()
        new_state[0] -= 1  # p decreases
        new_state[1] += 1  # pMHC increases
        return new_state

    def reaction4(state): # p-MHC unbinding
        new_state = state.copy()
        new_state[0] += 1  # p increases
        new_state[1] -= 1  # pMHC decreases
        return new_state

    def reaction5(state): # p-MHC exportation
        new_state = state.copy()
        new_state[1] -= 1  # pMHC decreases
        new_state[2] += 1  # pMHCp increases
        return new_state

    def reaction6(state): # pMHCp unbinding
        new_state = state.copy()
        new_state[2] -= 1  # pMHCp decreases
        return new_state
    
    # Define propensities
    def propensity1(state, parameters):
        return kin
    
    def propensity2(state, parameters):
        return kout * state[0]

    def propensity3(state, parameters):
        return kon_rhoM * state[0]
    
    def propensity4(state, parameters):
        return koff * state[1]

    def propensity5(state, parameters):
        return kstep * state[1]

    def propensity6(state, parameters):
        return koff * state[2]

    #Parameters
    kin = 10
    kout = 0.1
    kon_rhoM = 10
    ksteps = np.logspace(-2, 0, 3)
    koffs = np.logspace(-4, 5, 20)
    ptotal = kin/kout

    for kstep in ksteps:
        fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
        pfrees = []
        pMHCs = []
        pMHCps = []
        for koff in koffs:
            parameters = [kin, kout, kon_rhoM, koff, kstep]

            # Initial state [A, B]
            initial_state = [10, 10, 0]
            
            # Run Gillespie simulation
            times, states = gillespie_algorithm(initial_state, [reaction1, reaction2, reaction3, reaction4, reaction5, reaction6],
             [propensity1, propensity2, propensity3, propensity4, propensity5, propensity6], parameters, max_time=400)
            
            # # Print results
            # for t, s in zip(times, states):
            #     print(f"Time: {t:.4f}, State: {s}")

            pfree_t = [states[i][0] for i in range(len(times))]
            pMHC_t = [states[i][1] for i in range(len(times))]
            pMHCp_t = [states[i][2] for i in range(len(times))]

            pfrees.append(np.mean(pfree_t[int(len(times)/2):])/1)
            pMHCs.append(np.mean(pMHC_t[int(len(times)/2):])/1)
            pMHCps.append(np.mean(pMHCp_t[int(len(times)/2):])/1)

            # lastplot = ax.plot(times, pMHC, alpha = .8, label = r'$%.4f$'%koff)
            # ax.plot(times, pfree, ls = '--', color = lastplot[-1].get_color(), alpha = .8)
            # ax.hlines(1*kin/(kout+rhoM*kon-((rhoM*kon*koff)/(koff+kstep))), 0, times[-1], color = lastplot[-1].get_color())
            # ax.hlines(1*kin/(kout*(koff+kstep)/(rhoM*kon) + kstep), 0, times[-1], color = lastplot[-1].get_color())
            # ax.hlines(np.mean(pfree), 0, times[-1], color = lastplot[-1].get_color())
            # ax.hlines(np.mean(pMHC[int(len(times)/2):]), 0, times[-1], ls = ':', color = lastplot[-1].get_color())

        ax.plot(koffs, pfrees, color = my_blue, ls = '', marker = 'o', alpha = .8, label = r'$p$')
        ax.plot(koffs, pMHCs, color = my_green, ls = '', marker = 'o', alpha = .8, label = r'$pMHC$')
        ax.plot(koffs, pMHCps, color = my_red, ls = '', marker = 'o', alpha = .8, label = r'$pMHCp$')

        ax.plot(koffs, 1*kin/(kout+kon_rhoM-((kon_rhoM*koffs)/(koffs+kstep)))/1, color = my_blue, ls = '--', alpha = .8)
        ax.plot(koffs, 1*kin/(kout*(koffs+kstep)/(kon_rhoM) + kstep)/1, color = my_green, ls = '--', alpha = .8)
        ax.plot(koffs, (kstep/koffs)*kin/(kout*(koffs+kstep)/(kon_rhoM) + kstep)/1, color = my_red, ls = '--', alpha = .8)
        
        my_plot_layout(ax = ax, xscale='log', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30, bottom = 0.1, top = 8e3)
        ax.legend(fontsize = 20, loc = 1)
        fig.savefig('../../../Figures/presentation_MHC/SS_kstep-%.2f.pdf'%kstep)
