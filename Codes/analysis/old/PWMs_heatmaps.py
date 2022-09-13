import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
from Immuno_models import*
import scipy.special as sc
import pickle
from matplotlib import style
from scipy.optimize import curve_fit

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'

RT = 0.593

#Matrix = 'BLOSUM62'
Matrix = 'MJ2'
#Matrix = 'MM'


if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    M2_list = M2.tolist()
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    
if(Matrix == 'MM'):
    M2 = (np.loadtxt(Text_files_path + Matrix + '.txt')+1)*e0
    M2_list = M2.tolist()
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()

if(Matrix == 'BLOSUM62'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,25))
    M2_list = M2.tolist()
    Alphabet = np.array(['A', 'R'  ,'N' , 'D'  ,'C' , 'Q'  ,'E'  ,'G'  ,'H' , 'I'  ,'L'  ,'K'  ,'M' , 'F' , 'P' , 'S'  ,'T' , 'W' , 'Y' , 'V' , 'B' , 'Z'  ,'X',  '*'])
    Alphabet_list = Alphabet.tolist()
    
L_alphabet = len(Alphabet)
#----------------------------------------------------------------------



seq_names = ['H1', 'H3', 'MHC']
seqs = ['TFSDYWMNWV', 'GSYYGMDYWG', Alphabet[np.random.randint(0, len(Alphabet), 9)]]
file_names= ['PWM_Adams_etal_2016_1.pkl', 'PWM_Adams_etal_2016_2.pkl', 'Collesano_2022.csv']
file_exts = ['pkl', 'pkl', 'csv']
alphabet_data_Adams = np.array(['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'S', 'T', 'N', 'Q', 'C', 'P', 'H', 'K', 'R', 'D', 'E'])
alphabet_data_Collesano = np.array(['R','H','K','D','E','S','T', 'N', 'Q', 'C','G', 'P','A','V','I','L', 'M', 'F','Y','W'])
alphabets_data = np.array([alphabet_data_Adams, alphabet_data_Adams, alphabet_data_Collesano], dtype = object)
for j, seq_name in enumerate(seq_names[:]):

    fig, ax = plt.subplots(1, 3, figsize=(40,8))

    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=Matrix, ax = ax[0])

    antigen = seqs[j]

    L = len(seqs[j])

    antigen_list = [i for i in antigen]

    contributions = np.zeros(shape = (1,20))
    antigen_seq = np.array([], dtype = int)
    for i, aa in enumerate(antigen_list):
        index = Alphabet_list.index(aa)
        antigen_seq = np.append(antigen_seq, int(index))
    PWM = M2[:,antigen_seq]
    #loading PWM from file
    if(file_exts[j]=='pkl'):
        file_PWM = open(Text_files_path + 'Data/' + file_names[j], 'rb')
        log10Kd_WT = -8.91169118
        PWM_data = np.log(10**(pickle.load(file_PWM) - log10Kd_WT))*RT
        #PWM_data = (pickle.load(file_PWM)) - logKd_WT
        file_PWM.close()
    #loading PWM from Collesano 2022
    if(file_exts[j]=='csv'):
        PWM_data = pd.read_csv(Text_files_path + 'Data/'+ file_names[j]).to_numpy()
        PWM_data = PWM_data*np.log(10)
    #WT_seq = np.array(['T', 'F', 'S', 'D', 'Y', 'W', 'M', 'N', 'W', 'V'])
    #log10Kd_OPT = log10Kd_WT - ((PWM_H1_data[0, 2]) + (PWM_H1_data[15, 3]) + (PWM_H3_data[1, 1]) + (PWM_H3_data[9, 2]) + (PWM_H3_data[19, 6]) + (PWM_H3_data[4, 8]))

    #PWMs from data with the same order of aas as my MJ matrix
    PWM_data_2 = np.ones_like(PWM_data)

    for i in np.arange(L_alphabet):
        PWM_data_2[i,:] = PWM_data[np.where(Alphabet==alphabets_data[j][i])[0][0], :]
    #Change values by the minimum
    for i in np.arange(L):
        PWM[:,i]-=np.min(PWM[:,i], axis=0)
        PWM_data[:,i]-=np.min(PWM_data[:,i], axis=0)
        PWM_data_2[:,i]-=np.min(PWM_data_2[:,i], axis=0)

    PWM_list = PWM.tolist()
    PWM_data_list = PWM_data.tolist()
    PWM_data_2_list = PWM_data_2.tolist()

    plot_PWM(PWM=PWM_list, Alphabet=Alphabet, sequence = antigen, title=r'PWM', ax = ax[1])
    plot_PWM(PWM= PWM_data_2_list, Alphabet=Alphabet, sequence = seqs[j], title=seq_name , ax = ax[2])

    min_E = np.sum([np.min(PWM[:,i]) for i in range(len(PWM[0,:]))])
    avg_E = np.sum([np.mean(PWM[:,i]) for i in range(len(PWM[0,:]))])
    var_E = np.sum([np.var(PWM[:,i]) for i in range(len(PWM[0,:]))])
    max_E = np.sum([np.max(PWM[:,i]) for i in range(len(PWM[0,:]))])

    min_E_data = np.sum([np.min(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    avg_E_data = np.sum([np.mean(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    var_E_data = np.sum([np.var(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    max_E_data = np.sum([np.max(PWM_data[:,i]) for i in range(len(PWM_data[0,:]))])
    print(np.exp(min_E_data), np.exp(avg_E_data))
    #----------------------------------------------------------------------

    #----------------------------------------------------------------------

    # linear_PWM = np.reshape(PWM, (L*20,1))
    # linear_PWM_data = np.reshape(PWM_data, (L*20,1))
    # ax[3].hist(linear_PWM, bins = 20, alpha = .5, align='left', label = 'MJ')
    # ax[3].hist(linear_PWM_data, bins = 20, alpha = .5, align='left', label = 'Adams etal 2016')
    # my_plot_layout(ax=ax[3], xscale = 'linear')
    # ax[3].legend(loc = 0, fontsize = 20, title = r'Model', title_fontsize = 22)

    fig.savefig('../../Figures/6_Data/PWM_model_with_data_%s.pdf'%(seq_name))

