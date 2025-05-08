import sys
sys.path.append('../../lib/')
from funcs import*


# Parameters
l = 16  # Length of the sequence
L = 10**8  # Total sequences to simulate
energy_matrix = np.random.random((l, 20))  # Example energy matrix
beta = .5  # Initial inverse temperature
num_samples = int(10**6)  # Number of MCMC samples
energy_model = 'TCRen'
N_epi = 1
antigen = [15, 0, 0, 6, 10, 18, 1, 16, 12, 0, 13, 10, 5, 19, 14, 6]
# Calculate motif
motif = get_motif(antigen, energy_model, '../../')*1.2
E_ms = np.zeros(N_epi)

for epi in range(N_epi):
    E_m = -3
    # Normalize motif
    for i in range(l):
        E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
        motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
    print(E_m)
    E_ms[epi] = E_m
    # Calculate Q0, Es, dE, betas
    Es, dE, Q0, betas = calculate_Q0(0.01, 50, 400000, motif[:, epi*l:(epi+1)*l], E_m, l)

fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

start_time = time.time()
samples0 = np.array([calculate_energy(motif, np.random.randint(0, 20, size=l)) + E_ms[0] for i in range(int(L))])
end_time = time.time()
print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")
samples0 = samples0[samples0<-9]
ax.hist(samples0, density = False, bins = np.linspace(-26, -6, 20), alpha = .8, label = 0)

for i, k in enumerate([1, 10, 100, 1000]):
    # Run MCMC
    start_time = time.time()
    samples2 = mcmc_sampling(l, motif, beta, int(num_samples/k))
    end_time = time.time()
    print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")
    # Extract low-energy sequences
    # low_energy_samples = sorted(samples2, key=lambda x: x[1])[:10]  # Top 100 lowest energy
    samples2 = np.array([samples2[i][1] + E_ms[0] for i in range(len(samples2))])
    samples2 = samples2[samples2<-9]
    # samples2 = samples2[samples2>np.min(samples2)]
    # samples2 = samples2[samples2>-18]
    ax.hist(samples2, density = False, bins = np.linspace(-26, -6, 20), alpha = .4, label = i+1)


ax.legend()
ax.plot(Es[:-1], Q0*L)
ax.set_ylim(bottom = 1e-4)
ax.set_yscale('log')
plt.show()
