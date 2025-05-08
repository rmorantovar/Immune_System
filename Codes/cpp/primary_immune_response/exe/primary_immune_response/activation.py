import sys
sys.path.append('../library/')
from Immuno_models import*
plt.rcParams['text.usetex'] = True


Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/"


# Paramters:
N_engagements = 1e5
k_off_array = np.logspace(-5, 2, 19)
k_on = 1e6
k_pr_array = np.array([0.1, 1, 3]) #min^-1
#l_off = k_off
l_on = 0
k_act = 0.5 #hour^-1
k_act = k_act * 24

#States
AB = 1
ABp = 0
A_B = 0
Activation = 0

#The code simulates the different competitive processes depending on the current state:
for j, k_pr in enumerate(k_pr_array):
	file_p_a = open(Text_files_path+'recognition_k_act-%.0e_k_pr-%.0e.pkl'%(k_act/24, k_pr),'wb')
	p_a_array = np.array([])
	for i, k_off in enumerate(k_off_array):
		p_a = 0
		for n in tqdm(np.arange(N_engagements)):
			AB = 1
			ABp = 0
			A_B = 0
			Activation = 0
			while((1-A_B)*(1-Activation)):
				r = np.random.random() 
				propensities = np.array([k_off*AB, k_pr*AB, k_off*ABp, k_act*ABp])
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


