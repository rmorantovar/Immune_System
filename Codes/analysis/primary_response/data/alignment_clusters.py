import glob
import subprocess

clustered_files = glob.glob("clustered_alignments/*.fasta")

for cluster_fasta in clustered_files:
    output_aln = cluster_fasta.replace(".fasta", "_aligned.fasta")

    try:
        subprocess.run(["mafft", "--auto", cluster_fasta], stdout=open(output_aln, "w"), check=True)
        print(f"✅ Alignment successful for {cluster_fasta}! Output: {output_aln}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error aligning {cluster_fasta}: {e}")
