import sys
sys.path.append('../../lib/')
from functions_memory import*
# from classes import*
#from functions_2 import*
from matplotlib.colors import LogNorm

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
        for a1, antigen1 in enumerate(tqdm(antigens)):
            data1 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/activated_repertoire.csv')
            DDG_df = data1[['seq']].copy()
            for a11, antigen11 in enumerate(antigens):
                DDG_df = DDG_distributions(data1['E'], DDG_df, 'xxxx', 0, antigen1, a1+1, antigen11, a11+1, energy_model, N_epi, l)
                
            DDG_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/DDGs.csv', index = False)

            for a2, antigen2 in enumerate(antigens):
                data2 = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/%d'%(a2+1) + '/activated_repertoire.csv')
                DDG_df = data2[['seq']].copy()
                for a22, antigen22 in enumerate(antigens):
                    DDG_df = DDG_distributions(data2['E'], DDG_df, antigen1, a1+1, antigen2, a2+1, antigen22, a22+1, energy_model, N_epi, l)

                DDG_df.to_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/%d'%(a2+1) + '/DDGs.csv', index = False)
    

    for a1, antigen1 in enumerate(antigens):
        fig1, ax1 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.92})
        DDG_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a1+1) + '/DDGs.csv')
        for a2, antigen2 in enumerate(antigens):
            d = []
            for i in range(len(antigen1)):
                if antigen1[i]==antigen2[i]:
                    d.append(0)
                else:   
                    d.append(1)
            if sum(d)==1:
                ax1.hist(DDG_df[antigen2], bins = np.linspace(-2, 7, 20), alpha = .5, label = antigen2)

        my_plot_layout(ax = ax1, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
        ax1.legend(fontsize = 10, title_fontsize = 12, loc = 0, title = r'$\mathrm{Test}$')
        ax1.set_title('Currrent:' + antigen1,fontsize = 12)
        # ax1.set_xlim(left = 2e-6, right = 2e2)
        # ax1.set_ylim(bottom = 2e-6, top = 2e2)
        # ax1.set_xticks(range(-1, 6))
        # ax1.set_yticks(range(-1, 6))
        # ax1.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
        fig1.savefig('../../../Figures/memory/DeltaZ/DDG_distribution_1_'+energy_model+'_N_ens-'+str(N_ens)+'_N_evo-'+ str(N_evo) + '_sera-' + str(a1+1) + '.pdf')
        plt.close() 

    for a, antigen in enumerate(antigens):
        for a1, antigen1 in enumerate(antigens):
            fig2, ax2 = plt.subplots(figsize=(5*1.62, 5), gridspec_kw={'left':0.12, 'right':.9, 'bottom':.1, 'top': 0.92})
            DDG_df = pd.read_csv(root_dir + pars_dir_1 + pars_dir_2 + '/%d'%(a+1) + '/%d'%(a1+1) + '/DDGs.csv')
            for a2, antigen2 in enumerate(antigens):
                d = []
                for i in range(len(antigen1)):
                    if antigen1[i]==antigen2[i]:
                        d.append(0)
                    else:   
                        d.append(1) 
                if sum(d)==1:
                    ax2.hist(DDG_df[antigen2], bins = np.linspace(-2, 7, 20), alpha = .5, label = antigen2)
            
            my_plot_layout(ax = ax2, xscale='linear', yscale= 'linear', ticks_labelsize= 24, x_fontsize=30, y_fontsize=30)
            ax2.set_title('Past:' + antigen + ' ; Currrent:' + antigen1, fontsize = 12)
            ax2.legend(fontsize = 10, title_fontsize = 12, loc = 0, title = r'$\mathrm{Test}$')
            # ax2.set_xlim(left = 2e-6, right = 2e2)
            # ax2.set_ylim(bottom = 2e-6, top = 2e2)
            # ax2.set_xticks(range(-1, 6))
            # ax2.set_yticks(range(-1, 6))
            # ax2.set_xticklabels([r'$%d$'%i for i in range(1, N_inf+1, 8)])
            fig2.savefig('../../../Figures/memory/DeltaZ/DDG_distribution_2_'+energy_model+'_N_ens-'+str(N_ens)+'_N_evo-'+ str(N_evo) + '_past-' + str(a+1) + '_sera-' + str(a1+1) + '.pdf')
            plt.close()       
   
if __name__ == "__main__":
    main()

