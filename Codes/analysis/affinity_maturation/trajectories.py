import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Number of trajectories in the ensemble.
num_trajectories = 10

average_fitnesses = []
for i in range(num_trajectories):
    # Load the DataFrame from the CSV file.
    df = pd.read_csv(f"../../out/affinity_maturation/trajectory_{i}.csv")
    
    # Compute the average fitness at each time step.
    average_fitness = df.groupby('time')['fitness'].mean()
    average_fitnesses.append(average_fitness)

# Convert to a numpy array for easier manipulation.
average_fitnesses = np.array(average_fitnesses)

# Plot each trajectory.
for i in range(num_trajectories):
    plt.plot(average_fitnesses[i], label=f"Trajectory {i}")

# Plot the average trajectory.
mean_trajectory = average_fitnesses.mean(axis=0)
plt.plot(mean_trajectory, label="Mean trajectory", color="black", linewidth=2)

plt.xlabel("Time step")
plt.ylabel("Average fitness")
plt.legend()
plt.show()
