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
N_ensemble  = 1000
energy_model = "MM"

alphas = [1]
betas = [0.5, 1, 2]
ds = [10]
Ns = [2e3]
L = 15
e0 = 4

if(energy_model=="MM"):

	for d in ds:
		for beta in betas:
			for alpha in alphas:
				for N in Ns:
					for i, growth_model in enumerate(growth_models):
						print('Growth model: ',growth_model, '; Energy model: ',energy_model)
						file_clone_sizes = open(path+'clone_size_data_d-%d_beta-%.2f_N-%.1e_'%(d, beta, N)+growth_model+'_'+energy_model+'.pkl','wb')
						clone_sizes = np.array([])

						for n in tqdm(np.arange(N_ensemble)):
							
							#my_response = Immune_response(L=15, N=2000, alpha = alpha, beta=beta, antigen_str = 'FMLFMAVFVMTSWYC', text_files_path=path, energy_model = 'MJ', growth_model = growth_model)
							my_response = Immune_response(L=L, N=int(N), alpha = alpha, beta=beta,  text_files_path=path, energy_model = energy_model, d = d, e0=e0, bcells_filter = True)
							my_response.run(T = 25)
							clone_sizes = np.append(clone_sizes, my_response.B_cells_Tseries[my_response.activation_status,-1])

						pickle.dump(clone_sizes, file_clone_sizes)
						file_clone_sizes.close()

if(energy_model=="MJ"):
	d = 20
	for beta in betas:
		for alpha in alphas:
			for N in Ns:
				for i, growth_model in enumerate(growth_models):
					print('Growth model: ',growth_model)
					file_clone_sizes = open(path+'clone_size_data_d-%d_beta-%.2f_N-%.1e_'%(d, beta, N)+growth_model+'_'+energy_model+'.pkl','wb')
					clone_sizes = np.array([])

					for n in tqdm(np.arange(N_ensemble)):
						
						#my_response = Immune_response(L=15, N=2000, alpha = alpha, beta=beta, antigen_str = 'FMLFMAVFVMTSWYC', text_files_path=path, energy_model = 'MJ', growth_model = growth_model)
						my_response = Immune_response(L=L, N=int(N), alpha = alpha, beta=beta,  text_files_path=path, energy_model = energy_model, d = d, e0=e0, bcells_filter = True)
						my_response.run(T = 25)
						clone_sizes = np.append(clone_sizes, my_response.B_cells_Tseries[my_response.activation_status,-1])

					pickle.dump(clone_sizes, file_clone_sizes)
					file_clone_sizes.close()


