import sys
sys.path.append('../../lib/')
from functions_memory import*
# from classes import*
#from functions_2 import*

def main():
	parser = argparse.ArgumentParser(description="Generate random sequences and save properties to a CSV file.")
	parser.add_argument('--antigen', type=str, default='', help="Antigen sequence. Epitopes are separated by -")
	parser.add_argument('--N_ant', type=int, default=1, help="Number of antigens.")
	parser.add_argument('--N_ens', type=int, default=1, help="Number of times to execute the process.")
	parser.add_argument('--N_inf', type=int, default=1, help="Number of infections.")
	parser.add_argument('--N_evo', type=int, default = 0)
	parser.add_argument('--N_epi', type=int, default = 1)
	parser.add_argument('--l', type=int, default=16, help="Length of the sequences.")
	parser.add_argument('--new', type=int, default=0, help="run Z values again.")
	args = parser.parse_args()

	# Parameters -----------------------------------------------------

	N_ant = args.N_ant
	N_ens = args.N_ens
	N_inf = args.N_inf
	N_evo = args.N_evo
	N_epi = args.N_epi
	l = args.l
	new = args.new

	if N_evo == -1:
		N_evo = 'R'

	L0 = int(1e6)  # Number of random sequences
	E_lim = -11.  # Threshold for the sum of entries
	t_lim = 8.  # Threshold for the sum of entries
	chunk_size = 1e6  # Size of each chunk
	p = 3
	k_step = 720
	n_jobs = -1

	lamA = 6.0
	lamB = 3 * np.log(2) #(days)^-1
	lamB = 2.
	dT = 0.05
	C = 1e4
	time_array = np.linspace(0, 15, int((10-0)/dT))
	colors_inf = plt.cm.jet(np.linspace(0,1,N_inf))
	colors_mut = [my_blue, my_red]

	#----------------------------------------------------------------
	energy_model = 'TCRen'
	#energy_model = 'MJ2'
	# antigen = args.antigen
	# epitopes = antigen.split('-')
	# l=len(epitopes[0])
	
	project = 'memory'
	subproject = 'multi-epitope'
	subproject = 'PS'
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}"
	pars_dir_1 = f"/L0-{L0}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	pars_dir_2 = f"/N_ant-{N_ant}_N_ens-{N_ens}_N_epi-{N_epi}_N_evo-{N_evo}"
	antigens_data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/antigens.csv")
	antigens = antigens_data['antigen']

	N_ant = len(antigens)
	
	if new:
		sera1 = defaultdict(list)
		sera2 = defaultdict(list)
		for a1, antigen1 in enumerate(tqdm(antigens)):
			try:
				data1 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/activated_repertoire.csv')
				data2 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/%d'%(a1+1) + '/activated_repertoire.csv')
				sera1 = run_essay_homo(data1, 'xxxx', 0, antigen1, a1+1, antigen1, a1+1, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, time_array, sera1)
				sera2 = run_essay_homo(data2, antigen1, a1+1, antigen1, a1+1, antigen1, a1+1, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, time_array, sera2)
			except FileNotFoundError:
				print(f'skipping {antigen1}')
				continue					
		sera1_df = pd.DataFrame(sera1)
		sera1_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Zt12_1.csv', index = False)
		sera2_df = pd.DataFrame(sera2)
		sera2_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Zt12_2.csv', index = False)

	fig, ax = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	fig1, ax1 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	fig2, ax2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

	sera1_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Zt12_1.csv')
	sera2_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/Zt12_2.csv')
	sera21_df = sera2_df.copy()
	sera21_df['Z12'] = sera2_df['Z12']/sera1_df['Z12']
	sera21_df['t12'] = sera2_df['t12'] - sera1_df['t12']
 
	sns.scatterplot(data = sera21_df, x = 't12', y = 'Z12', hue = 'current id', ax = ax)
	sns.scatterplot(data = sera1_df, x = 't12', y = 'Z12', hue = 'current id', ax = ax1)
	sns.scatterplot(data = sera2_df, x = 't12', y = 'Z12', hue = 'current id', ax = ax2)

	my_plot_layout(ax = ax, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# ax.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# ax.set_xlim(left = 4, right = 14)
	# ax.set_ylim(bottom = 1e8, top = 1.5e11)
	# ax.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig.savefig('../../../Figures/memory/DeltaZ/Zt12_'+energy_model+'_N_ant-'+str(N_ant)+'_N_ens-'+str(N_ens)+'.pdf')

	my_plot_layout(ax = ax1, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# ax1.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	ax1.set_xlim(left = 4, right = 14)
	ax1.set_ylim(bottom = 1e8, top = 1.5e11)
	# ax1.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax1.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig1.savefig('../../../Figures/memory/DeltaZ/Zt12_1_'+energy_model+'_N_ant-'+str(N_ant)+'_N_ens-'+str(N_ens)+'.pdf')

	my_plot_layout(ax = ax2, xscale='linear', yscale= 'log', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# ax2.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	ax2.set_xlim(left = 4, right = 14)
	ax2.set_ylim(bottom = 1e8, top = 1.5e11)
	# ax2.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax2.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig2.savefig('../../../Figures/memory/DeltaZ/Zt12_2_'+energy_model+'_N_ant-'+str(N_ant)+'_N_ens-'+str(N_ens)+'.pdf')


if __name__ == "__main__":
    main()

