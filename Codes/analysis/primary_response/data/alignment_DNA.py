import subprocess

input_fasta = "Sequence_all.fasta"
output_aln = "Sequence_aligned_all.fasta"

# Run Clustal Omega with --force to allow overwriting
try:
    subprocess.run(["clustalo", "-i", input_fasta, "-o", output_aln, "--auto", "-v", "--force"], check=True)
    print(f"✅ Alignment successful! Output saved in '{output_aln}'")
except subprocess.CalledProcessError as e:
    print(f"❌ Error running Clustal Omega: {e}")
