import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
from Immuno_models import*
import scipy.special as sc
import pickle
from matplotlib import style
from scipy.optimize import curve_fit

Text_files_path = '/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/'



#Matrix = 'BLOSUM62'
Matrix = 'MJ2'
#Matrix = 'MM'

fig, ax = plt.subplots(1, 4, figsize=(40,8))

if(Matrix == 'MJ2'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,21))
    M2_list = M2.tolist()
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'MJ Model', ax = ax[0])
if(Matrix == 'MM'):
    M2 = (np.loadtxt(Text_files_path + Matrix + '.txt')+1)*e0
    M2_list = M2.tolist()
    Alphabet = np.loadtxt(Text_files_path + 'Alphabet.txt', dtype=bytes, delimiter='\t').astype(str)
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'MJ Model', ax = ax[0])
if(Matrix == 'BLOSUM62'):
    M2 = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 1, usecols=range(1,25))
    M2_list = M2.tolist()
    Alphabet = np.array(['A', 'R'  ,'N' , 'D'  ,'C' , 'Q'  ,'E'  ,'G'  ,'H' , 'I'  ,'L'  ,'K'  ,'M' , 'F' , 'P' , 'S'  ,'T' , 'W' , 'Y' , 'V' , 'B' , 'Z'  ,'X',  '*'])
    Alphabet_list = Alphabet.tolist()
    plot_energy_matrix(Energy_Matrix=M2_list, Alphabet=Alphabet, title=r'BLOSUM62 Model', ax = ax[0])
L_alphabet = len(Alphabet)
#----------------------------------------------------------------------

#antigen_seq = np.random.randint(0, len(Alphabet), L)
#antigen = Alphabet[antigen_seq]
antigen = 'MYMFQALNSL'
L = len(antigen)
antigen_list = [i for i in antigen]
contributions = np.zeros(shape = (1,20))
antigen_seq = np.array([], dtype = int)
for i, aa in enumerate(antigen_list):
    index = Alphabet_list.index(aa)
    antigen_seq = np.append(antigen_seq, int(index))
PWM = M2[:,antigen_seq]
#loading PWM from Adams et.al.
file_PWM_H1 = open(Text_files_path + '/Data/PWM_Adams_etal_2016_1.pkl', 'rb')
file_PWM_H3 = open(Text_files_path + '/Data/PWM_Adams_etal_2016_2.pkl', 'rb')

logKd_WT = -8.91169118

PWM_H1_data = (pickle.load(file_PWM_H1)) - logKd_WT
PWM_H3_data = (pickle.load(file_PWM_H3)) - logKd_WT
file_PWM_H1.close()
file_PWM_H3.close()


alphabet_data = np.array(['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'S', 'T', 'N', 'Q', 'C', 'P', 'H', 'K', 'R', 'D', 'E'])
WT_seq = np.array(['T', 'F', 'S', 'D', 'Y', 'W', 'M', 'N', 'W', 'V'])

logKd_OPT = logKd_WT - ((PWM_H1_data[0, 2]) + (PWM_H1_data[15, 3]) + (PWM_H3_data[1, 1]) + (PWM_H3_data[9, 2]) + (PWM_H3_data[19, 6]) + (PWM_H3_data[4, 8]))

#PWMs from data with the same order of aas as my MJ matrix
PWM_H1_data_2 = np.ones_like(PWM_H1_data)
PWM_H3_data_2 = np.ones_like(PWM_H3_data)

for i in np.arange(L_alphabet):
    PWM_H1_data_2[i,:] = PWM_H1_data[np.where(Alphabet==alphabet_data[i])[0][0], :]
    PWM_H3_data_2[i,:] = PWM_H3_data[np.where(Alphabet==alphabet_data[i])[0][0], :]
#Change values by the minimum
for i in np.arange(L):
    PWM[:,i]-=np.min(PWM[:,i], axis=0)
    PWM_H1_data[:,i]-=np.min(PWM_H1_data[:,i], axis=0)
    PWM_H1_data_2[:,i]-=np.min(PWM_H1_data_2[:,i], axis=0)
    PWM_H3_data[:,i]-=np.min(PWM_H3_data[:,i], axis=0)
    PWM_H3_data_2[:,i]-=np.min(PWM_H3_data_2[:,i], axis=0)

PWM_list = PWM.tolist()
PWM_H1_data_list = PWM_H1_data.tolist()
PWM_H1_data_2_list = PWM_H1_data_2.tolist()
PWM_H3_data_list = PWM_H3_data.tolist()
PWM_H3_data_2_list = PWM_H3_data_2.tolist()

plot_PWM(PWM=PWM_list, Alphabet=Alphabet, sequence = antigen, title=r'PWM', ax = ax[1])
plot_PWM(PWM=PWM_H1_data_2_list, Alphabet=Alphabet, sequence = np.array(['#']*L), title=r'CDRH1', ax = ax[2])
plot_PWM(PWM=PWM_H3_data_2_list, Alphabet=Alphabet, sequence = np.array(['#']*L), title=r'CDRH3', ax = ax[3])

min_E = np.sum([np.min(PWM[:,i]) for i in range(len(PWM[0,:]))])
avg_E = np.sum([np.mean(PWM[:,i]) for i in range(len(PWM[0,:]))])
var_E = np.sum([np.var(PWM[:,i]) for i in range(len(PWM[0,:]))])
max_E = np.sum([np.max(PWM[:,i]) for i in range(len(PWM[0,:]))])

min_E_data = np.sum([np.min(PWM_H3_data[:,i]) for i in range(len(PWM_H3_data[0,:]))])
avg_E_data = np.sum([np.mean(PWM_H3_data[:,i]) for i in range(len(PWM_H3_data[0,:]))])
var_E_data = np.sum([np.var(PWM_H3_data[:,i]) for i in range(len(PWM_H3_data[0,:]))])
max_E_data = np.sum([np.max(PWM_H3_data[:,i]) for i in range(len(PWM_H3_data[0,:]))])
print(np.exp(min_E_data), np.exp(avg_E_data))
#----------------------------------------------------------------------

#----------------------------------------------------------------------

# linear_PWM = np.reshape(PWM, (L*20,1))
# linear_PWM_data = np.reshape(PWM_data, (L*20,1))
# ax[3].hist(linear_PWM, bins = 20, alpha = .5, align='left', label = 'MJ')
# ax[3].hist(linear_PWM_data, bins = 20, alpha = .5, align='left', label = 'Adams etal 2016')
# my_plot_layout(ax=ax[3], xscale = 'linear')
# ax[3].legend(loc = 0, fontsize = 20, title = r'Model', title_fontsize = 22)

fig.savefig('../../Figures/6_Data/PWM_model_with_data.pdf')

