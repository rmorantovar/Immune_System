import numpy as np
import matplotlib.pyplot as plt

# Define state space (2 states for simplicity)
n_states = 10

# Define base free energies (Delta G)
np.random.seed(0)
delta_G = np.linspace(1, 4, n_states)

# Define perturbation delta (increase in energy)
delta = np.linspace(0.1, 1.0, n_states)[::-1]  # higher delta for lower energy states

# Define l and l' (correlated but different)
l = np.exp(-delta_G)                 # biased toward low energy
l = l / np.sum(l)

l_prime = np.exp(-0.5 * delta_G)     # flatter than l, still correlated
l_prime = l_prime / np.sum(l_prime)

# Define temperatures
beta = 1.5

# Compute W and Z before and after perturbation
W = np.sum(l * np.exp(-delta_G))
W_new = np.sum(l * np.exp(-delta_G - delta))

Z = np.sum(l_prime * np.exp(-beta * delta_G))
Z_new = np.sum(l_prime * np.exp(-beta * (delta_G + 0.5*delta)))

# Compute relative decays
decay_W = W_new / W
decay_Z = (Z_new / Z)**(1) #this is the new exponent to translate into viral load

# Plotting
labels = ['W', 'Z']
decays = [decay_W, decay_Z]

plt.bar(labels, decays, color=['grey', 'black'])
plt.ylim(0, 1.1)
plt.title("Relative Decay After Perturbation")
plt.ylabel("Fraction Remaining")
plt.grid(axis='y')
plt.text(0, decay_W + 0.02, f"{decay_W:.2f}", ha='center')
plt.text(1, decay_Z + 0.02, f"{decay_Z:.2f}", ha='center')
plt.show()
