import subprocess

my_list = ['CGG', 'OVA', 'NP-OVA', 'HA']
for my_string in my_list:
    input_fasta = "CDR3_"+my_string+".fasta"
    output_aln = "CDR3_aligned_"+my_string+".fasta"

    try:
        subprocess.run(["mafft", "--auto", input_fasta], stdout=open(output_aln, "w"), check=True)
        print(f"✅ MAFFT alignment successful! Output saved in '{output_aln}'")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running MAFFT: {e}")
