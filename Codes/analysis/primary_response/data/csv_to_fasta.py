import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/primary_response/data/victora_2020"

my_list = ['CGG', 'OVA', 'NP-OVA', 'HA']
my_list = ['all', 'larger', 'larger2']
for my_string in my_list:
    # Load Excel file (modify 'sequences.xlsx' and sheet name if needed)
    file_path = root_dir + "/data_cluster_"+my_string+".csv"  # Change this to your actual file
    df = pd.read_csv(file_path, header = 0)  # Adjust sheet name if necessary

    print(df)
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
    fasta_file = root_dir + "/CDR3_"+my_string+".fasta"
    with open(fasta_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    print(f"FASTA file '{fasta_file}' created successfully!")
