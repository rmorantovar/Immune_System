import sys
sys.path.append('../../lib/')
from funcs import*
import random

# Parameters
N = 100          # Number of sequences
L = 15          # Length of each segment
total_len = 3*L
aa_list = list("ACDEFGHIKLMNPQRSTVWY")
mutations_per_generation = (1, 3)  # Range of mutations per new sequence

def mutate_sequence(seq, num_mutations):
    seq = list(seq)
    positions = random.sample(range(len(seq)), num_mutations)
    for pos in positions:
        original = seq[pos]
        new_aa = random.choice([aa for aa in aa_list if aa != original])
        seq[pos] = new_aa
    return ''.join(seq)
    

# Step 1: Create the root (ancestral) sequence
root_seq = ''.join(random.choices(aa_list, k=total_len))

# Step 2: Create a hierarchy of related sequences
sequences = [root_seq]
while len(sequences) < N:
    # Pick a random existing sequence to mutate
    parent = random.choice(sequences)
    num_muts = random.randint(*mutations_per_generation)
    child = mutate_sequence(parent, num_muts)
    sequences.append(child)

# Step 3: Output final sequences
M = np.zeros((N, N))
for i, seq in enumerate(sequences):
    # print(f">Seq{i+1}\n{seq}")
    for j, seq2 in enumerate(sequences[i:]):
        # print(i, j+i, hamming_distance(seq, seq2))
        M[i, j+i] = hamming_distance(seq, seq2)
        M[j+i, i] = M[i, j+i]

print(M)




