import sys
sys.path.append('../../my_lib/')
from funcs import*
# from classes import*
#from functions_2 import*
from matplotlib.colors import LogNorm

def main():
	# Set up command-line argument parser
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")

	# Define command-line arguments with explanations
	parser.add_argument('--N_ant', type=int, default=100, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default=0, help="Evolution count.")
	parser.add_argument('--N_epi', type=int, default=3, help="Number of epitopes.")
	parser.add_argument('--L0', type=int, default=10**7, help="Number of random sequences.")
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--t_lim', type=float, default=8.0, help="Activation time threshold.")
	parser.add_argument('--E_lim', type=float, default=-6.0, help="Threshold for the sum of entries.")
	parser.add_argument('--E_m', type=float, default=-24, help="Energy threshold.")
	parser.add_argument('--chunk_size', type=int, default=1000000, help="Size of each chunk.")
	parser.add_argument('--p', type=float, default=4.0, help="Number of steps.")
	parser.add_argument('--k_step', type=float, default=720, help="Step rate.")
	parser.add_argument('--lamA', type=float, default=6.0, help="Antigen growth rate A.")
	parser.add_argument('--lamB', type=float, default=2.0, help="Antigen growth rate B.")
	parser.add_argument('--n_jobs', type=int, default=-1, help="Number of jobs for parallel processing.")
	parser.add_argument('--random_antigen', type=int, default=0, help="Random antigen flag.")
	parser.add_argument('--antigen', type=str, default='TACNSYPNTAKCRWYR', help="Antigen sequence.")
	parser.add_argument('--energy_model', type=str, default='TCRen', help="Energy model.")
	parser.add_argument('--seqs', type=int, default=1, help="Number of sequences.")
	parser.add_argument('--one_WT', type=int, default=0, help="Single WT flag.")
	parser.add_argument('--secondary', type=int, default=0, help="Secondary infection flag.")
	parser.add_argument('--secondary_all', type=int, default=1, help="Secondary all infections flag.")
	parser.add_argument('--add_mutant', type=int, default=0, help="add mutant.")
	parser.add_argument('--pro', type=str, default='epitope_complexity', help="Project name.")
	parser.add_argument('--subpro', type=str, default='epistasis', help="Subproject name.")
	parser.add_argument('--exp', type=int, default=0, help="Experiment ID.")
	parser.add_argument('--potency_all', type = int, default = 1)

	args = parser.parse_args()

	# ------------ PARAMETERS AND INPUTS ------------
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
	if L0>=1e6:
		chunk_size = args.chunk_size
	else:
		chunk_size = args.L0
	p = args.p
	k_step = args.k_step
	lamA = args.lamA
	lamB = args.lamB
	n_jobs = args.n_jobs
	random_antigen = args.random_antigen
	antigen = args.antigen
	energy_model = args.energy_model
	
	seqs = args.seqs
	one_WT = args.one_WT
	secondary = args.secondary
	secondary_all = args.secondary_all
	potency_all = args.potency_all
	add_mutant = args.add_mutant

	if N_evo == -1:
		N_evo = 'R'

	dT = 0.001
	C = 1e4
	T = 12
	time_array = np.linspace(0, T, int((T-0)/dT))
	Alphabet = np.loadtxt('../../in/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)

	project = args.pro
	subproject = args.subpro
	experiment = args.exp
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	pars_dir_1 = f"/L0-{int(L0/10**int(np.log10(L0)))}e{int(np.log10(L0))}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}"#_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv", converters={"antigen": literal_eval})
	antigens = antigens_data['antigen']
	# print(antigens)

	# ------------ ------------ ------------

	start_time = time.time()
	print('Starting simulation ...')

	if one_WT:
		WTs = [antigens.iloc[0]]
	else:
		WTs = list(antigens)
	
	colors = [my_blue, my_red, my_green, my_cyan]
	lss = ['-', '--']
	lws = [1, 2]

	effects_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/DDG_epis.csv')
	effects0_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/DDG_epis_epi.csv')
	
	output_plot = '/Users/robertomorantovar/Dropbox/My_Documents/Science/Projects/Immune_System/_Repository/Figures/'+project+'/'+subproject
	os.makedirs(output_plot, exist_ok=True)

	mutations = ['_0_3_8', '_16_13_11', '_32_7_8']
	
	for mutation in mutations:
		# for kappa in effects_df['kappa'].unique():
		for alpha in effects_df['alpha'].unique():
			fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.10, 'right':.95, 'bottom':.1, 'top': 0.95})
			
			# ax.hist(effects_df[(effects_df['mut'].str.endswith(mutation)) & (effects_df['kappa'] == kappa)]['DDG'], bins = 10, density = True, color = my_red, label = r'$\kappa$ fixed')
			ax.hist(effects_df[(effects_df['mut'].str.endswith(mutation)) & (effects_df['alpha'] == alpha)]['DDG'], bins = 10, density = True, color = my_green, label = r'$\alpha$ fixed')
			ax.hist(effects0_df[(effects0_df['mut'].str.endswith(mutation)) & (effects0_df['DDG'] != 0.0)]['DDG'], bins = 10, density = True, alpha = .5, color = 'k')

			my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
			# ax.set_xticks([])
			# ax.set_yticks([])
			# ax.set_ylim(bottom = -0.3, top = 7.3)
			# ax.set_xlim(left = -0.3, right = 7.3)
			ax.legend(fontsize = 16, loc = 0, title = r'$p_m$', title_fontsize = 18)
			# fig.savefig(output_plot + '/hist_DDG_kappa_' + str(kappa) + mutation + '.pdf')
			fig.savefig(output_plot + '/hist_DDG_alpha_' + str(alpha) + mutation + '.pdf')
			plt.close()

	# Print Final execution time
	end_time = time.time()
	print(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()