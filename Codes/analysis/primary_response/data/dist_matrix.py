from Bio import AlignIO, SeqIO
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
import os

# Load aligned sequences
aligned_fasta = "CDR3_aligned.fasta"  # Make sure this file exists from Step 1
alignment = AlignIO.read(aligned_fasta, "fasta")
sequences = [str(record.seq) for record in alignment]
seq_ids = [record.id for record in alignment]

# Convert gaps to a standard character to avoid issues
sequences = [seq.replace("-", "X") for seq in sequences]

# Compute pairwise Hamming distances
def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2)) / len(s1)

dist_matrix = np.array([[hamming_distance(s1, s2) for s2 in sequences] for s1 in sequences])

# Convert to condensed format for clustering
condensed_dist_matrix = squareform(dist_matrix)

# Perform hierarchical clustering
Z = linkage(condensed_dist_matrix, method="ward")  # Ward's method minimizes cluster variance

# Set the expected number of clusters (adjustable)
num_clusters = 5  # Adjust as needed
clusters = fcluster(Z, num_clusters, criterion="maxclust")

# # ✅ Print cluster assignments
# print("Cluster assignments:")
# for seq_id, cluster in zip(seq_ids, clusters):
#     print(f"{seq_id} → Cluster {cluster}")

# Save clusters into separate FASTA files
output_folder = "clustered_alignments"
os.makedirs(output_folder, exist_ok=True)

# Store sequences by cluster
cluster_dict = {i: [] for i in range(1, num_clusters + 1)}

for record, cluster_id in zip(alignment, clusters):
    cluster_dict[cluster_id].append(record)

# ✅ Save each cluster to a separate FASTA file
for cluster_id, records in cluster_dict.items():
    output_fasta = os.path.join(output_folder, f"cluster_{cluster_id}.fasta")
    with open(output_fasta, "w") as out_f:
        SeqIO.write(records, out_f, "fasta")
    print(f"✅ Saved Cluster {cluster_id} with {len(records)} sequences in '{output_fasta}'")
