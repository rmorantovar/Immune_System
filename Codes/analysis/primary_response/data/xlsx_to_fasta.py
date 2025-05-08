import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Load Excel file (modify 'sequences.xlsx' and sheet name if needed)
file_path = "/Users/robertomorantovar/Dropbox/Research/Immune_system/primary_response/data/victora_2020/mmc1.xlsx"  # Change this to your actual file
df = pd.read_excel(file_path, sheet_name="Photoactivation CGG", header = 0, skiprows = 1)  # Adjust sheet name if necessary

# Assuming sequences are in a column named 'Sequence' and IDs in 'ID' (adjust if needed)
sequence_column = "CDR3:"
id_column = "Sequence ID"

# Remove duplicate sequences while keeping the first occurrence
df_unique = df.drop_duplicates(subset=[sequence_column])

# Check the first few rows after removing duplicates
print(f"✅ Unique sequences found: {len(df_unique)}")
print(df_unique.head())

# Create a list of SeqRecord objects for FASTA format
records = []
for index, row in df_unique.iterrows():
    sequence = str(row[sequence_column])  # Convert to string in case of NaN values
    seq_id = str(row[id_column]) if id_column in df.columns else f"seq_{index+1}"
    records.append(SeqRecord(Seq(sequence), id=seq_id, description=""))

# Save to FASTA file
fasta_file = "CDR3.fasta"
with open(fasta_file, "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")

print(f"FASTA file '{fasta_file}' created successfully!")
