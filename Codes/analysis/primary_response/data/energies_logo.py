import numpy as np
import pandas as pd
from Bio import AlignIO

my_list = ['CGG', 'OVA', 'NP-OVA', 'HA']
for my_string in my_list:
	# Load the MSA
	alignment_file = "CDR3_aligned_"+my_string+".fasta"  # Change to your file
	alignment = AlignIO.read(alignment_file, "fasta")

	# Convert alignment to a list of sequences
	sequences = [str(record.seq) for record in alignment]
	seq_ids = [record.id for record in alignment]

	# Define amino acids
	amino_acids = "ACDEFGHIKLMNPQRSTVWY"
	# amino_acids = "CTGA"

	# Compute frequency matrix
	alignment_length = len(sequences[0])
	freq_matrix = np.zeros((len(amino_acids), alignment_length))

	for seq in sequences:
	    for i, aa in enumerate(seq):
	        if aa in amino_acids:
	            freq_matrix[amino_acids.index(aa), i] += 1

	# Normalize to probability
	freq_matrix /= len(sequences)

	# Assume uniform background probabilities (1/20 for amino acids)
	background = np.full(len(amino_acids), 1 / len(amino_acids))

	# Compute PWM (log-odds scores)
	pwm = np.log2((freq_matrix + 1e-9) / background[:, None])  # Avoid log(0)

	# Convert PWM to DataFrame
	df_pwm = pd.DataFrame(pwm, index=list(amino_acids), columns=[f"Pos {i+1}" for i in range(alignment_length)])

	# Function to compute PWM score for a sequence
	def compute_pwm_score(sequence, pwm, amino_acids):
	    score = 0
	    for i, aa in enumerate(sequence):
	        if aa in amino_acids:
	            score += pwm.loc[aa, f"Pos {i+1}"]
	    return score

	# Example ΔG values (hypothetical, in kcal/mol)
	delta_g = np.random.uniform(-2, 2, size=(len(amino_acids), alignment_length))  # Replace with real values if available

	# Function to compute ΔG-based motif energy
	def compute_motif_energy(sequence, pwm, delta_g, amino_acids):
	    energy = 0
	    for i, aa in enumerate(sequence):
	        if aa in amino_acids:
	            energy += pwm.loc[aa, f"Pos {i+1}"] * delta_g[amino_acids.index(aa), i]
	    return -energy  # Negative sign follows ΔG convention

	# Compute scores and energies for all sequences
	sequence_scores = []
	for seq_id, seq in zip(seq_ids, sequences):
	    pwm_score = compute_pwm_score(seq, df_pwm, amino_acids)
	    motif_energy = compute_motif_energy(seq, df_pwm, delta_g, amino_acids)
	    sequence_scores.append({"Sequence ID": seq_id, "Sequence": seq, "PWM Score": pwm_score, "Motif Energy (ΔG)": motif_energy})

	# Convert to DataFrame
	df_scores = pd.DataFrame(sequence_scores)

	# Save to CSV file
	file_path = "/Users/robertomorantovar/Dropbox/Research/Immune_system/primary_response/data/victora_2016/"  # Change this to your actual file
	output_csv = file_path + "CDR3_aligned_"+my_string+"_motif_scores.csv"
	df_scores.to_csv(output_csv, index=False)

	print(f"✅ Scores saved to {output_csv}")
