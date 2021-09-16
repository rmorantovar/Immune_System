import sys
sys.path.append('../library/')
from Immuno_models import *
import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import pickle
import time
from tqdm import tqdm

def neutralization(a, b, x):
	return 1/(1+((np.exp(a*(x-1.2)))/(b)))

N_epitopes = np.arange(1,10,2)

N_ensemble = 1000

T = 8

fig, ax = plt.subplots(figsize=(12, 8), )

for l in np.arange(len(N_epitopes)):
	time = np.linspace(0, T, T*N_epitopes[l]+1)
	cross_protection = np.zeros_like(time)
	for n in np.arange(N_ensemble):
		mutations = np.random.randint(low=1, high=N_epitopes[l]+1, size = T*N_epitopes[l])
		cross_protection_temp = np.array([np.product([1-neutralization(2.5, 10, np.count_nonzero(mutations[0:t]==i)) for i in N_epitopes[0: l+1]]) for t in np.arange(len(time))])
		cross_protection = cross_protection + cross_protection_temp

	if(N_epitopes[l]==1):
		cross_protection_1 = np.array([(1-neutralization(2.5, 10, np.count_nonzero(mutations[0:t]==1))) for t in np.arange(len(time))])
		time_1 = time


	plot1 = ax.plot(time, cross_protection/N_ensemble, label = '%d'%(N_epitopes[l]))
	ax.plot(time_1, cross_protection_1**N_epitopes[l], linestyle = 'dashed', color = plot1[0].get_color())


my_plot_layout(ax)
ax.legend()
plt.show()
	