import subprocess

input_fasta = "CDR3.fasta"  # Your input file
output_fasta = "clusters.fasta"  # Output clustered file

try:
    subprocess.run([
        "cd-hit",
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", "0.4",  # 🔹 Lower similarity threshold (50%) to reduce cluster count
        "-n", "2",  # 🔹 Keeps sensitivity for short sequences
        "-d", "0",  # 🔹 Avoids unnecessary description truncation
        "-M", "16000",  # 🔹 Uses more memory for better clustering
        "-T", "8"  # 🔹 Use 8 threads (modify based on your CPU)
    ], check=True)
    print(f"✅ Clustering successful! Output saved in '{output_fasta}'")
except subprocess.CalledProcessError as e:
    print(f"❌ Error running CD-HIT: {e}")
