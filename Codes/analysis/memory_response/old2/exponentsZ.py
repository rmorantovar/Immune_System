import sys
sys.path.append('../../lib/')
from functions_memory import*
# from classes import*
#from functions_2 import*
from matplotlib.colors import LogNorm

def main():

	project = 'memory_response'
	subproject = 'multi-epitope'
	subproject = 'Panel'
	experiment = 0
	root_dir = f"/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/{experiment}"
	
	p = 3
	k_step = 720
	lamA = 6.0
	lamB = 3 * np.log(2) #(days)^-1
	lamB = 2.

	fig1, ax1 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	fig2, ax2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})
	# fig3, ax3 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.96})

	df1 = pd.DataFrame(columns=['beta','exp_Z'])
	df2 = pd.DataFrame(columns=['beta','exp_Z'])
	# df = pd.DataFrame(columns=['beta','alpha1','alpha2'])

	pars_dir_1 = f"/L0-{int(1e6)}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	# for N_ens in [100]:
	# 	for N_evo in ['R']:	
	# 		pars_dir_2 = f"/N_ant-{100}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
	# 		data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])
			
	# 		data1 = data.loc[data['past id']==0]
	# 		df1 = pd.concat([df1, data1])
	# 		sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

	# 		data2 = data.loc[data['past id']!=0]
	# 		df2 = pd.concat([df2, data2])
	# 		sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)
	# for N_ens in [4, 5]:
	# 	for N_evo in [0, 1, 2, 3, 'R']:	
	# 		pars_dir_2 = f"/N_ant-{10}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
	# 		data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])

	# 		data1 = data.loc[data['past id']==0]
	# 		df1 = pd.concat([df1, data1])
	# 		sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

	# 		data2 = data.loc[data['past id']!=0]
	# 		df2 = pd.concat([df2, data2])
	# 		sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)

	# for L0 in [int(1e6), int(1e7)]:
	# 	pars_dir_1 = f"/L0-{L0}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	# 	for N_ens in [10]:
	# 		for N_evo in ['R']:	
	# 			pars_dir_2 = f"/N_ant-{100}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
	# 			data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])
				
	# 			data1 = data.loc[data['past id']==0]
	# 			df1 = pd.concat([df1, data1])
	# 			sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

	# 			data2 = data.loc[data['past id']!=0]
	# 			df2 = pd.concat([df2, data2])
	# 			sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)
	for N_ens in [40]:
		for N_evo in [0, 1, 2, 3, 'R']:	
			pars_dir_2 = f"/N_ant-{10}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
			data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])

			data1 = data.loc[data['past id']==0]
			df1 = pd.concat([df1, data1])
			sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

			data2 = data.loc[data['past id']!=0]
			df2 = pd.concat([df2, data2])
			sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)

	pars_dir_1 = f"/L0-{int(1e7)}_p-{p}_k_step-{k_step}_lamA-{lamA}_lamB-{lamB}"
	for N_ens in [40]:
		for N_evo in [0, 1, 2, 3, 'R']:	
			pars_dir_2 = f"/N_ant-{10}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
			data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])

			data1 = data.loc[data['past id']==0]
			df1 = pd.concat([df1, data1])
			sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

			data2 = data.loc[data['past id']!=0]
			df2 = pd.concat([df2, data2])
			sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)
	# for N_ens in [10]:
	# 	for N_evo in ['R']:	
	# 		pars_dir_2 = f"/N_ant-{20}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
	# 		data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])
			
	# 		data1 = data.loc[data['past id']==0]
	# 		df1 = pd.concat([df1, data1])
	# 		sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

	# 		data2 = data.loc[data['past id']!=0]
	# 		df2 = pd.concat([df2, data2])
	# 		sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)
	for N_ens in [200]:
		for N_evo in ['R']:	
			pars_dir_2 = f"/N_ant-{100}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
			data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])
			
			data1 = data.loc[data['past id']==0]
			df1 = pd.concat([df1, data1])
			sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

			data2 = data.loc[data['past id']!=0]
			df2 = pd.concat([df2, data2])
			sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)
	for N_ens in [100]:
		for N_evo in ['R']:	
			pars_dir_2 = f"/N_ant-{200}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
			data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])
			
			data1 = data.loc[data['past id']==0]
			df1 = pd.concat([df1, data1])
			sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

			data2 = data.loc[data['past id']!=0]
			df2 = pd.concat([df2, data2])
			sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)
	for N_ens in [20]:
		for N_evo in ['R']:	
			pars_dir_2 = f"/N_ant-{1000}_N_ens-{N_ens}_N_epi-{1}_N_evo-{N_evo}"
			data = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + "/Z_time.csv", usecols = ['past id', 'beta', 'exp_Z'])
			
			data1 = data.loc[data['past id']==0]
			df1 = pd.concat([df1, data1])
			sns.scatterplot(data = data1, x = 'beta',  y = 'exp_Z', ax = ax1, alpha = 0.4)

			data2 = data.loc[data['past id']!=0]
			df2 = pd.concat([df2, data2])
			sns.scatterplot(data = data2, x = 'beta',  y = 'exp_Z', ax = ax2, alpha = 0.4)

	
	v = lamA/p

	bins = np.linspace(np.min(df1['beta']*.9), np.max(df1['beta']*1.01), 18)
	cuts = pd.cut(df1['beta'], bins = bins, labels = (bins[:-1]+bins[1:])/2)
	df_avg1 = df1.groupby([cuts]).exp_Z.mean().reset_index(name='exp_Z_mean')
	sns.scatterplot(data = df_avg1, x = 'beta',  y = 'exp_Z_mean', ax = ax1, s = 40, color = 'black')
	ax1.hlines(lamB, 0.6, 2, color = 'k', ls = '--')
	beta_array = np.linspace(0.6, np.max(df1['beta']/0.9), 50)
	beta_array = np.linspace(0.6, 3.2, 100)
	sigma = beta_array - 1 - lamB/v
	t_mes = 1.5
	Z = np.exp(lamB*t_mes)*(np.exp(v*sigma*t_mes) - 1)/(sigma)
	Zprime = v*np.exp((lamB + v*sigma)*t_mes) + (lamB*np.exp(lamB*t_mes)*(-1+np.exp(v*sigma*t_mes)))/(sigma)
	ax1.plot(beta_array, lamB + v*sigma, color = 'k', ls = '--')
	ax1.plot(beta_array, Zprime/Z, color = 'k', ls = '-')
	
	
	bins = np.linspace(np.min(df2['beta']*.9), np.max(df2['beta']*1.01), 18)
	cuts = pd.cut(df2['beta'], bins = bins, labels = bins[:-1]+ np.diff(bins)[0]/2)
	df_avg2 = df2.groupby([cuts]).exp_Z.mean().reset_index(name='exp_Z_mean')
	sns.scatterplot(data = df_avg2, x = 'beta',  y = 'exp_Z_mean', ax = ax2, s = 40, color = 'black')
	ax2.hlines(lamB, 0.9, 3, color = 'k', ls = '--')
	beta_array = np.linspace(0.9/0.9, np.max(df2['beta']/0.9), 50)
	sigma = beta_array - 1 - 2*lamB/v
	t_mes = 2.
	Z = np.exp(lamB*t_mes)*(np.exp(v*sigma*t_mes) - 1)/(sigma)
	Zprime = v*np.exp((lamB + v*sigma)*t_mes) + (lamB*np.exp(lamB*t_mes)*(-1+np.exp(v*sigma*t_mes)))/(sigma)
	ax2.plot(beta_array, (lamB + v*sigma), color = 'k', ls = '--')
	ax2.plot(beta_array, (Zprime/Z), color = 'k', ls = '-')

	#sns.scatterplot(data = df, x = 'alpha1',  y = 'alpha2', hue = 'beta', ax = ax3)
	
	my_plot_layout(ax = ax1, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# ax1.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# ax1.set_xlim(left = 2e-6, right = 2e2)
	ax1.set_ylim(bottom = 1.8, top = 4.4)
	# ax1.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax1.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig1.savefig('../../../Figures/memory_response/DeltaZ/Z_exponents1_Panel_'+'TCRen'+'.pdf')
	
	
	my_plot_layout(ax = ax2, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# ax2.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# ax2.set_xlim(left = 2e-6, right = 2e2)
	ax2.set_ylim(bottom = 1.8, top = 5.4)
	# ax2.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# ax2.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	fig2.savefig('../../../Figures/memory_response/DeltaZ/Z_exponents2_Panel_'+'TCRen'+'.pdf')
	
	
	# my_plot_layout(ax = ax3, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
	# # ax3.legend(fontsize = 18, title_fontsize = 22, loc = 2, title = r'$\mathrm{Infection}$')
	# # ax3.set_xlim(left = 2e-6, right = 2e2)
	# # ax3.set_ylim(bottom = 2e-6, top = 2e2)
	# # ax3.set_xticks(np.linspace(0, time_array[-1]*(N_inf-1), int(N_inf/8)))
	# # ax3.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
	# fig3.savefig('../../../Figures/memory/DeltaZ/Z_exponents12_Panel_'+'TCRen'+'.pdf')



if __name__ == "__main__":
    main()

