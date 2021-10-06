import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import pickle
import time
from tqdm import tqdm

path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Python/ensemble/"

growth_models = ['exponential', 'linear']
growth_models = ['exponential']
N_ensemble  = 500

alphas = [1]
betas = [1]
for beta in betas:
	for alpha in alphas:
		for i, growth_model in enumerate(growth_models):
			print('Growth model: ',growth_model)
			file_clone_sizes = open(path+'clone_size_data_alpha-%.2f_beta-%.2f_'%(alpha, beta)+growth_model+'.pkl','wb')
			clone_sizes = np.array([])

			for n in tqdm(np.arange(N_ensemble)):
				
				my_response = Immune_response(L=15, N=2000, alpha = alpha, beta=beta, antigen_str = 'FMLFMAVFVMTSWYC', text_files_path=path, energy_model = 'MJ',
					growth_model = growth_model)
				my_response.run(T = 20)
				clone_sizes = np.append(clone_sizes, my_response.B_cells_Tseries[my_response.activation_status,-1])

			pickle.dump(clone_sizes, file_clone_sizes)
			file_clone_sizes.close()

