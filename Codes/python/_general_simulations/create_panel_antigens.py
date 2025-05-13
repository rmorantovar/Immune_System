import sys
sys.path.append('../../lib/')
from funcs import*


parser = argparse.ArgumentParser(description="Generate random antigens and save properties to a CSV file.")
parser.add_argument('--N_ant', type=int, default=80, help="Number of antigens.")
parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
parser.add_argument('--N_evo', type=int, default = 1)
parser.add_argument('--N_epi', type=int, default = 3)
parser.add_argument('--L0', type=int, default=10**8, help="Number of random antigens.")
parser.add_argument('--l', type=int, default=16, help="Length of the antigens.")
parser.add_argument('--t_lim', type=float, default=8., help="Threshold for activation time.")
parser.add_argument('--E_lim', type=float, default=-11., help="Threshold for the sum of entries.")
parser.add_argument('--E_m', type=float, default=-24, help="Threshold for the sum of entries.")
parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
parser.add_argument('--p', type=float, default=4.0, help="# steps.")
parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
parser.add_argument('--lamA', type=float, default=6., help="Antigen growth rate.")
parser.add_argument('--lamB', type=float, default=2., help="Antigen growth rate.")
parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
parser.add_argument('--random_antigen', type=int, default=0)
# parser.add_argument('--antigen', type=str, default='TACNSEYPNTTRAKCGRWYR')
parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR')
parser.add_argument('--energy_model', type=str, default = 'TCRen')
parser.add_argument('--seqs', type=int, default = 1)
parser.add_argument('--pro', type=str, default='epitope_complexity', help="project.")
parser.add_argument('--subpro', type=str, default='minimal_ID', help="subproject.")
parser.add_argument('--exp', type=int, default=0, help="experiment.")
args = parser.parse_args()


# Retrieving arguments
N_ant = args.N_ant
N_ens = args.N_ens
N_inf = args.N_inf
N_evo = args.N_evo
N_epi = args.N_epi
L0 = args.L0
l = args.l
E_lim = args.E_lim
t_lim = args.t_lim
E_m = args.E_m
chunk_size = args.chunk_size
p = args.p
k_step = args.k_step
lamA = args.lamA
lamB = args.lamB
n_jobs = args.n_jobs
random_antigen = args.random_antigen
antigen = args.antigen
energy_model = args.energy_model

seqs = args.seqs

dT = 0.05
C = 1e4
T = 12
time_array = np.linspace(0, T, int((T-0)/dT))

Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)


project = 'epitope_complexity'
subproject = 'epistasis'
experiment = args.exp
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"

output_dir = root_dir + pars_dir_1 + pars_dir_2
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'antigens.csv')

mutations_per_generation = (1, 2)  # Range of mutations per new sequence

def mutate_sequence(seq, num_mutations):
    seq = list(seq)
    mutable_positions = [i for i in range(len(seq))]
    for i in range(N_epi):
        mutable_positions.remove(i*l)
    positions = random.sample(mutable_positions, num_mutations)
    for pos in positions:
        # print('epitope:', pos//l)
        original = seq[pos]
        new_aa = random.choice([aa for aa in Alphabet if aa != original])
        seq[pos] = new_aa
    return ''.join(seq)

# Step 1: Create the root (ancestral) sequence
root_seq = ''.join(random.choices(Alphabet, k=N_epi*l))

# Step 2: Create a hierarchy of related antigens
antigens = [root_seq]
while len(antigens) < N_ant:
    # Pick a random existing sequence to mutate
    parent = random.choice(antigens)
    num_muts = random.randint(*mutations_per_generation)
    # print(num_muts)
    child = mutate_sequence(parent, num_muts)
    antigens.append(child)

# Step 3: Output final antigens
M = np.zeros((N_ant, N_ant))
for i, seq in enumerate(antigens):
    # print(f">Seq{i+1}\n{seq}")
    for j, seq2 in enumerate(antigens[i:]):
        # print(hamming_distance(seq, seq2))
        M[i, j+i] = hamming_distance(seq, seq2)
        M[j+i, i] = M[i, j+i]

print(M)
antigens_int = []
for antigen in antigens:
	antigens_int.append(from_aa_to_i_Alphabet(Alphabet, antigen))
antigens_df = pd.DataFrame([])
antigens_df['antigen'] = antigens_int
antigens_df.to_csv(output_file, index=False)

