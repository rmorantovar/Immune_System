import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy



path = "../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Single_trajectory/"
my_response = Immune_response(L=15, N=5000, alpha = 1, beta=0.5, antigen_str = 'FMLFMAVFVMTSWYC', text_files_path=path, energy_model = 'MJ')
my_response.run(T = 20)

plt.plot(my_response.time, my_response.antigen_Tseries)
plt.show()
plt.stackplot(my_response.time, my_response.B_cells_Tseries[my_response.activation_status]/np.sum(my_response.B_cells_Tseries[my_response.activation_status], axis = 0))
plt.show()