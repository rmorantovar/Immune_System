import subprocess
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/primary_response/data/victora_2020"

my_list = ['CGG', 'OVA', 'NP-OVA', 'HA']
my_list = ['all', 'larger', 'larger2']
for my_string in my_list:
    input_fasta = root_dir + "/CDR3_"+my_string+".fasta"
    output_aln = root_dir + "/CDR3_aligned_"+my_string+".fasta"

    try:
        subprocess.run(["mafft", "--auto", input_fasta], stdout=open(output_aln, "w"), check=True)
        print(f"✅ MAFFT alignment successful! Output saved in '{output_aln}'")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running MAFFT: {e}")
