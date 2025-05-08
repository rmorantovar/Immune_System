import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity


def generate_random_points(N, M, dimensions):
    # Generate N random points in the d-dimensional hypercube
    points_set_1 = np.random.rand(N, dimensions)
    
    # Generate M random points in the same d-dimensional hypercube
    points_set_2 = np.random.rand(M, dimensions)
    
    return points_set_1, points_set_2

def calculate_distances(points_set_1, points_set_2):
    # Calculate distances between all points in both sets
    distances = np.linalg.norm(points_set_1[:, None] - points_set_2, axis=2)
    return distances

# Define parameters
N_A = 100  # Number of points in the first set
N_B = 10  # Number of points in the second set
dimensions = 5  # Dimensionality of the space

# Generate random points in the hypercube
points_set_A, points_set_B = generate_random_points(N_A, N_B, dimensions)
# Calculate distances between points in the sets
distances_matrix = calculate_distances(points_set_A, points_set_B)
distances_matrix_A = calculate_distances(points_set_A, points_set_A)



# Assuming 'distances_matrix' is your NxM distance matrix from earlier

# Calculate cosine similarity between rows of the matrix
cos_sim = cosine_similarity(distances_matrix)

# Display the cosine similarity matrix
print("Cosine Similarity Matrix:")
print(cos_sim, distances_matrix_A)


# # Plotting the heatmap
# plt.imshow(distances_matrix, cmap='viridis', interpolation='nearest')
# plt.colorbar(label='Distance')
# plt.xticks(np.arange(N_B), labels=np.arange(1, N_B + 1))
# plt.yticks(np.arange(N_A), labels=np.arange(1, N_A + 1))
# plt.tight_layout()
# plt.show()


# Plotting the heatmap
plt.scatter(distances_matrix_A, 1-cos_sim)
plt.tight_layout()
plt.show()