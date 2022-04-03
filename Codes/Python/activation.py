import sys
sys.path.append('../library/')
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from Immuno_models import*
import scipy.special as sc
import pickle

Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/"


# Paramters:
N_engagements = 1e7
k_off_array = np.logspace(-5, 4, 19)
k_on = 1e6
m_off = 0
m_on_array = np.array([1e-2, 1,  1e2])
#l_off = k_off
l_on = 0
gamma = 1

#States
AB = 1
ABp = 0
A_B = 0
Activation = 0

#The coda simulates the different competitive processes depending on the current state:
for j, m_on in enumerate(m_on_array):
	file_p_a = open(Text_files_path+'recognition_gamma-%.0e_m_on-%.0e.pkl'%(gamma, m_on),'wb')
	p_a_array = np.array([])
	for i, k_off in enumerate(k_off_array):
		p_a = 0
		for n in np.arange(N_engagements):
			AB = 1
			ABp = 0
			A_B = 0
			Activation = 0
			while((1-A_B)*(1-Activation)):
				r = np.random.random() 
				propensities = np.array([k_off*AB, m_on*AB, k_off*ABp, gamma*ABp])
				cumsum = propensities.cumsum()
				alpha = propensities.sum()

				transitionId = np.searchsorted(cumsum,r*alpha)%2
				
				A_B = AB*(1 - transitionId) + ABp*(1 - transitionId)
				Activation = ABp*transitionId
				ABp = AB*transitionId
				AB = 0

			p_a+=Activation

		p_a_array = np.append(p_a_array, p_a/N_engagements)

	pickle.dump(np.array([k_off_array, p_a_array]), file_p_a)
	file_p_a.close()


