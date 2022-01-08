import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy



path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Single_trajectory/"

#my_response = Immune_response(L=15, N=5000, alpha = 1, beta=0.5, antigen_str = 'FMLFMAVFVMTSWYC', text_files_path=path, energy_model = 'MJ')
my_response = Immune_response(L=15, N=int(2e3), alpha = 1, beta=1,  text_files_path=path, energy_model = 'MM', d = 10, e0=4, bcells_filter = True)
my_response.run(T = 25)

#print(my_response.E_matrix)
plt.hist(my_response.energies)
#plt.ylim(0, 20)
plt.yscale('log')
plt.show()


#plt.plot(my_response.time, my_response.antigen_Tseries)
#plt.yscale('log')
#plt.show()
if(np.sum(np.where(my_response.activation_status==True))):	
	plt.stackplot(my_response.time, my_response.B_cells_Tseries[my_response.activation_status]/np.sum(my_response.B_cells_Tseries[my_response.activation_status], axis = 0))
	plt.show()