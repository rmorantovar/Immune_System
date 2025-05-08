import glob
import subprocess

aligned_files = glob.glob("clustered_alignments/*_aligned.fasta")

for aln_file in aligned_files:
    output_logo = aln_file.replace("_aligned.fasta", ".pdf")

    try:
        subprocess.run([
            "weblogo", "-f", aln_file, "-o", output_logo, "--format", "PDF", "--title", aln_file
        ], check=True)
        print(f"✅ Sequence logo saved: {output_logo}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error generating logo for {aln_file}: {e}")
