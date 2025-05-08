import sys
sys.path.append('../../lib/')
from funcs import*

parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
parser.add_argument('--N_ant', type=int, default=10, help="Number of antigens.")
parser.add_argument('--N_ens', type=int, default=10, help="Number of times to execute the process.")
parser.add_argument('--N_inf', type=int, default=2, help="Number of infections.")
parser.add_argument('--N_evo', type=int, default = -1)
parser.add_argument('--N_epi', type=int, default = 1)
parser.add_argument('--L0', type=int, default=10**8, help="Number of random sequences.")
parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
parser.add_argument('--t_lim', type=float, default=8., help="Threshold for activation time.")
parser.add_argument('--E_lim', type=float, default=-11., help="Threshold for the sum of entries.")
parser.add_argument('--E_m', type=float, default=-23, help="Threshold for the sum of entries.")
parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
parser.add_argument('--p', type=float, default=4, help="# steps.")
parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
parser.add_argument('--lamA', type=float, default=6., help="Antigen growth rate.")
parser.add_argument('--lamB', type=float, default=2., help="Antigen growth rate.")
parser.add_argument('--n_jobs', type=int, default=-1, help="n_jobs.")
parser.add_argument('--random_antigen', type=int, default=0)
# parser.add_argument('--antigen', type=str, default='TACNSEYPNTTRAKCGRWYR')
parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR')
parser.add_argument('--energy_model', type=str, default = 'TCRen')
parser.add_argument('--seqs', type=int, default = 1)
parser.add_argument('--pro', type=str, default='memory_response', help="project.")
parser.add_argument('--subpro', type=str, default='Z_dynamics', help="subproject.")
parser.add_argument('--exp', type=int, default=1, help="experiment.")
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

Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)

project = args.pro
subproject = args.subpro
experiment = args.exp
root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}_N_evo-R"

output_dir = root_dir + pars_dir_1 + pars_dir_2
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'antigens.csv')

antigens = []
for a in range(N_ant):
	antigen = np.random.choice(range(20), l*N_epi)
	# antigen = ''.join(np.random.choice(Alphabet[:-1], l*N_epi))
	# antigen = ("".join([antigen[epi*l:(epi+1)*l] + '-' for epi in range(N_epi)]))[:-1]
	antigens.append(list(antigen))

antigens_df = pd.DataFrame([])
antigens_df['antigen'] = antigens
antigens_df.to_csv(output_file, index=False)


