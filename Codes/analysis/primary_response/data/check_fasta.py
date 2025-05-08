from Bio import SeqIO

fasta_file = "sequences.fasta"

valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY*X-")  # Standard amino acids + '*' (stop) + 'X' (unknown)

for record in SeqIO.parse(fasta_file, "fasta"):
    invalid_chars = set(record.seq.upper()) - valid_amino_acids
    if invalid_chars:
        print(f"⚠️ Warning: Sequence {record.id} contains invalid characters: {invalid_chars}")
    if len(record.seq) < 10:  # Short sequences might cause issues
        print(f"⚠️ Warning: Sequence {record.id} is very short ({len(record.seq)} amino acids).")
