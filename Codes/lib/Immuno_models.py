import numpy as np
import matplotlib.pyplot as plt
#from Bio import Phylo
from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
import scipy.special as sc
import seaborn as sns
import pickle
import json
from scipy.optimize import curve_fit


class Sequence():
	"""docstring for Sequence"""
	def __init__(self, seq_id, master_sequence, energy_parent, complementary_sequence, Energy_Matrix, parent,  parent_id = 0, Master_Seq = False):
		super(Sequence, self).__init__()

		#A given Sequence object has the following attributes
		self.id = seq_id
		self.Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
		self.parent = parent
		self.parent_id = parent_id
		self.master_sequence = master_sequence
		self.complementary_sequence = complementary_sequence
		self.active = False
		self.clone_size = 1
		self.energy = 0
		self.energy_parent = energy_parent
		self.delta_energy = 0
		self.sequence = parent
		self.pos_mut = 0
		self.tree_position = 1 # 0 = internal ; 1 = external
		self.hamming_distance = 0
		if not(Master_Seq):

			self.mutate_sequence(Energy_Matrix = Energy_Matrix)
		if (Master_Seq):
			self.energy = self.energy_parent

	def calculate_delta_energy(self, Energy_Matrix, old_letter, new_letter):

		M = Energy_Matrix

		list_comp_seq = list(self.complementary_sequence)

		pos_new_letter = np.where(np.isin(self.Alphabet,new_letter))[0][0]
		pos_old_letter = np.where(np.isin(self.Alphabet,old_letter))[0][0]

		pos_comp_seq = np.where(np.isin(self.Alphabet,list_comp_seq[self.pos_mut]))[0][0]

		self.delta_energy = M[pos_comp_seq][pos_new_letter]-M[pos_comp_seq][pos_old_letter]
		self.energy = self.energy_parent + self.delta_energy


	def mutate_sequence(self, Energy_Matrix):
		""" This function will create a new mutations and give an energy value to the new sequence according to the Energy_Matrix. """
		self.pos_mut = np.random.randint(9)
		list_seq = list(self.sequence)
		old_letter = self.sequence[self.pos_mut]
		self.Alphabet.remove(old_letter)
		new_letter = np.random.choice(self.Alphabet)
		list_seq[self.pos_mut] = new_letter
		self.sequence = "".join(list_seq)
		self.hamming_distance = hamming_distance(self.master_sequence, self.sequence)
		self.Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
		self.calculate_delta_energy(Energy_Matrix = Energy_Matrix, old_letter = old_letter, new_letter = new_letter)
	

#----------------- Models -----------------

#Inside the classes for models, there are some built-in functions for generating plots of individual simulations.

class Deterministic_simulation():
	"""docstring for Deterministic_simulation"""
	def __init__(self, Sequences, n_linages, T, U, nu, beta, gamma, energy_translation, initial_time, dt):
		super(Deterministic_simulation, self).__init__()
		self.n_linages = n_linages
		self.Sequences = Sequences
		self.T = T
		self.U = U
		self.nu = nu #Linages growth rate
		self.beta = beta #Antigen growth rate
		self.gamma = gamma #Antigen clearance rate
		self.energy_translation = energy_translation
		self.dt = dt
		self.N_A = 6.02214076e23

		self.initial_time  = initial_time
		self.linages_time_series = np.ones(shape =(n_linages, int((self.T-self.initial_time)/self.dt)))
		self.activation_time_series = np.zeros(shape=(n_linages, int((self.T-self.initial_time)/self.dt)))
		self.active_linages = 0
		self.antigen_time_series = np.ones(int((self.T-self.initial_time)/self.dt)) 
		self.time_series = np.linspace(self.initial_time, self.T, int((self.T-self.initial_time)/self.dt))

	def ODE(self):
		#ANTIGEN
		self.antigen_time_series[0] = np.exp(self.beta*self.initial_time)
		for t in range(1,int((self.T-self.initial_time)/self.dt)):
			N_t_active = np.sum(self.linages_time_series[(np.sum(self.activation_time_series, axis=1)!=0),:], axis = 0) - np.sum(self.linages_time_series[(np.sum(self.activation_time_series, axis=1)!=0),:], axis = 0)[0]
			self.antigen_time_series[t] = self.antigen_time_series[t-1] + (self.beta*self.antigen_time_series[t-1] - self.gamma*self.antigen_time_series[t-1]*N_t_active[t-1])*self.dt
			if(self.antigen_time_series[t]<(1)):
				self.antigen_time_series[t] = 0
			#BCELL LINAGES 
			for i in range(self.n_linages):
				self.linages_time_series[i,t] = self.linages_time_series[i,t-1] + self.nu*self.linages_time_series[i,t-1]*self.Sequences[i].active*self.dt
				f = (self.antigen_time_series[t]/self.N_A)/((self.antigen_time_series[t]/self.N_A)+np.exp(self.energy_translation+self.Sequences[i].energy)) 
				if(f > 0.5):
					self.Sequences[i].active = 1
				if(self.Sequences[i].active == 1):
					self.activation_time_series[i, t] = 1

	def plot_antigen_time(self, color, ax):

		ax.plot(self.time_series, self.antigen_time_series/self.N_A, linewidth  = 4, color = color)
		ax.set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'Antigen $\rho$ [M]', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.set_ylim(np.exp(self.beta*self.initial_time)/self.N_A, (np.max(self.antigen_time_series)/self.N_A)*5)
		#handles, labels = ax.get_legend_handles_labels()
		#ax.legend(np.concatenate(([Line2D([0], [0], color='tab:red', linewidth=4, linestyle='solid', ms = 8)],handles)),np.concatenate(([r'$n_b(r, \rho)$'],labels)), loc = 0, fontsize = 20)

	def plot_prob_binding(self, ax):
		rho_array = np.logspace(int(np.log10(max(self.antigen_time_series)/self.N_A)) - 4, np.log10(max(self.antigen_time_series)/self.N_A)	, 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		energies  = np.array([self.Sequences[i].energy for i in range(int(len(self.Sequences)))])
		energies_array = np.linspace(np.min(energies),np.max(energies),100)
		for i, rho in enumerate(rho_array): 
			ax.plot(energies_array, (1/(1+np.exp(self.energy_translation + energies_array - np.log(rho)))), linewidth  = 4, color = colors[i], label = r'$\rho \approx 10^{%.0d}$'%(np.log10(rho)))
		ax.set_yscale('log')
		ax.set_xlabel(r'Energy $\epsilon$', fontsize = 20)
		ax.set_ylabel(r'Probability of binding $p_b$', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def stackplot_linages_time(self, ax, antigen = False, time = True):
		colors = []
		for i in self.Sequences:
			if(i.active==True):
				colors.append('indianred')
			else:
				colors.append('indigo')
		if(time):
			ax.stackplot(self.time_series, self.linages_time_series/np.sum(self.linages_time_series, axis = 0), colors=colors, alpha = 0.9);
			ax.set_xlabel(r'Time $t$', fontsize = 20)
		if(antigen):
			ax.stackplot(self.antigen_time_series, self.linages_time_series/np.sum(self.linages_time_series, axis = 0), colors=colors, alpha = 0.9);
			ax.set_xlabel(r'Antigen $\rho$ [M]', fontsize = 20)
			ax.set_xscale('log')
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'B cell Linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		#ax.legend(loc = 0, fontsize = 20)

	def hist_sequences_hamming_distance(self, Sequences, ax):

		rho_array = np.logspace(np.log10(1/self.N_A), np.log10(max(self.antigen_time_series)/self.N_A), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		data_distances = ax.hist([Sequences[i].hamming_distance for i in range(int(len(Sequences)))], bins = range(10), align = 'left', label = r'$S(d)$', color = 'olive', alpha = 0.4)
		#ax.plot(data_distances[1][0:-1], sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), linewidth = 4 , color = 'lightsteelblue', alpha = 0.6)

		ax.hist([self.Sequences[i].hamming_distance for i in range(int(len(self.Sequences)))], bins = range(10), align = 'left', label = r'$US(d)$', color = 'indigo', alpha = 0.6)
		#ax.plot(data_distances[1][0:-1], self.U*sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), linewidth = 4 , color = 'indigo', alpha = 0.6)

		ax.hist([self.Sequences[i].hamming_distance for i in range(int(len(self.Sequences))) if self.Sequences[i].active], bins = range(10), align = 'left', label = r'Activated Linages', color = 'indianred', alpha = 0.8)
		#for i, rho in enumerate(rho_array):
		#	ax.plot(data_distances[1][0:-1], self.U*sc.comb(9, data_distances[1][0:-1].astype(int))*((20-1)**data_distances[1][0:-1])*(1/(1+np.exp(self.master_Sequence_energy + data_distances[1][0:-1] - np.log(rho)))) , color = colors[i], linestyle = 'dashed', linewidth = 3)

		ax.set_ylim(0.1, 2e5)    
		ax.set_yscale('log')
		ax.set_xlabel(r'Hamming Distance $d$', fontsize = 20)
		ax.set_ylabel(r'Number of linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def hist_sequences_energy(self, Sequences, n_bins, ax):

		rho_array = np.logspace(np.log10(1/self.N_A), np.log10(max(self.antigen_time_series)/self.N_A), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		energies = np.array([Sequences[i].energy for i in range(int(len(Sequences)))])
		data_energies = ax.hist(energies, bins = np.linspace(np.min(energies), np.max(energies), n_bins), align = 'left', label = r'$S(\epsilon)$', color = 'lightsteelblue', alpha = 0.5)

		sub_energies = np.array([self.Sequences[i].energy for i in range(int(len(self.Sequences)))])
		ax.hist(sub_energies , bins = np.linspace(np.min(energies), np.max(energies), n_bins), align = 'left', label = r'$US(\epsilon)$', color = 'indigo', alpha = 0.6)

		sub_energies_activated = np.array([self.Sequences[i].energy for i in range(int(len(self.Sequences))) if self.Sequences[i].active])
		ax.hist(sub_energies_activated, bins = np.linspace(np.min(energies), np.max(energies), n_bins), align = 'left', label = r'Activated Linages', color = 'indianred', alpha = 0.8)
  
		ax.set_yscale('log')
		ax.set_xlabel(r'Energy $\epsilon$', fontsize = 20)
		ax.set_ylabel(r'Number of linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def hist_sequences_k_D(self, Sequences, n_bins, ax):

		rho_array = np.logspace(np.log10(1/self.N_A), np.log10(max(self.antigen_time_series)/self.N_A), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		energies = np.array([Sequences[i].energy for i in range(int(len(Sequences)))])
		data_energies = ax.hist(np.exp(energies), bins = np.logspace(np.log10(np.exp(np.min(energies))), np.log10(np.exp(np.max(energies))), n_bins), align = 'mid', label = r'$S(\epsilon)$', color = 'olivedrab', alpha = 0.5)

		sub_energies = np.array([self.Sequences[i].energy for i in range(int(len(self.Sequences)))])
		ax.hist(np.exp(sub_energies), bins = np.logspace(np.log10(np.exp(np.min(energies))), np.log10(np.exp(np.max(energies))), n_bins), align = 'mid', label = r'$US(\epsilon)$', color = 'indigo', alpha = 0.6)

		sub_energies_activated = np.array([self.Sequences[i].energy for i in range(int(len(self.Sequences))) if self.Sequences[i].active])
		ax.hist(np.exp(sub_energies_activated), bins = np.logspace(np.log10(np.exp(np.min(energies))), np.log10(np.exp(np.max(energies))), n_bins), align = 'mid', label = r'Activated Linages', color = 'indianred', alpha = 0.8)
		
		ax.set_yscale('log')
		ax.set_xlabel(r'$k_D$', fontsize = 20)
		ax.set_ylabel(r'Number of linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def plot_k_largest_linages(self, k, ax):

		#Calculate array of the frequencies of the largest k linages
		Seq_sizes = self.linages_time_series[:,-1]
		k = k
		biggest_k_linages_sizes = np.sort(Seq_sizes)[-k:]
		Pos = np.array([i for i, j in enumerate(Seq_sizes) if np.isin(j,biggest_k_linages_sizes)])
		biggest_k_linages = self.linages_time_series[Pos,:]
		biggest_k_linages_freq = biggest_k_linages/np.sum(biggest_k_linages, axis = 0)

		ax.stackplot(self.time_series, biggest_k_linages_freq, alpha = 0.9);
		#ax[0].set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'k largest linages', fontsize = 20)
		ax.tick_params(labelsize = 20)

		return biggest_k_linages_freq

	def plot_entropy_k_largest_linages(self, k, biggest_k_linages_freq, ax):

		#Calculate entropy
		entropy = np.array([np.sum(-1*biggest_k_linages_freq[:,t]*np.log(biggest_k_linages_freq[:,t])) for t in range(int(len(self.time_series)))])
		ax.plot(self.time_series[::50], entropy[::50], marker = 'o', ms = 8, linestyle = '', linewidth = '4', color = 'olive', label = 'Simulation')
		#ax[1].set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'Entropy', fontsize = 20)
		ax.tick_params(labelsize = 20)

	def plot_normalized_entropy_k_largest_linages(self, k, biggest_k_linages_freq, ax):

		#Calculate entropy
		entropy = np.array([np.sum(-1*biggest_k_linages_freq[:,t]*np.log(biggest_k_linages_freq[:,t])) for t in range(int(len(self.time_series)))])
		normalized_entropy = entropy/np.log(len(biggest_k_linages_freq))
		ax.plot(self.time_series[::50], normalized_entropy[::50], marker = 'o', ms = 8, linestyle = '', linewidth = '4', color = 'olive', label = 'Simulation')
		#ax[1].set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'Entropy', fontsize = 20)
		ax.tick_params(labelsize = 20)

		return normalized_entropy

	def plot_entropy_drop(self, k_array, ax1, ax2):

		#____________ Create array for entropy drop. The entropy is always normalized by the maximum

		entropy_drop = np.array([0])

		for i, k in enumerate(k_array[1:]):

		    biggest_k_linages_freq = self.plot_k_largest_linages(k=k, ax=ax1[i,0])
		    normalized_entropy = self.plot_normalized_entropy_k_largest_linages(k=k, biggest_k_linages_freq = biggest_k_linages_freq, ax=ax1[i,1])
		    entropy_drop = np.append(entropy_drop, normalized_entropy[-1] - normalized_entropy[0])

		
		#____________ Plot entropy drop as a fucntion of k
		ax2.plot(k_array[1:], abs(entropy_drop)[1:], linestyle = '--', marker = '^', ms = 20, linewidth = 4, color = 'indianred', alpha = 0.8, label = 'Simulation')
		ax2.set_xlabel(r'$k$', fontsize = 20)
		ax2.set_ylabel(r'$\Delta S$', fontsize = 20)
		ax2.tick_params(labelsize = 22)
		ax2.set_xscale('log')
		ax2.legend(loc = 0, fontsize = 20)
		
class Stochastic_simulation():
	"""docstring for Stochastic_simulation"""
	def __init__(self, Sequences, n_linages, T, U, gamma, nu, R, beta, master_Sequence_energy):
		super(Stochastic_simulation, self).__init__()
		self.n_linages = n_linages
		self.Sequences = Sequences
		self.T = T
		self.U = U
		self.gamma = gamma
		self.nu = nu
		self.R = R
		self.beta = beta
		self.master_Sequence_energy = master_Sequence_energy
		self.N_A = 6.02214076e23

		self.linages_time_series = np.ones(shape =(n_linages, 1))
		self.activation_time_series = np.zeros(shape=(n_linages, 1))
		self.active_linages = 0
		self.antigen_time_series = np.array([20])
		self.time_series = np.array([0])
		self.probabilities = np.zeros((2*n_linages)+1)

	def calculate_probabilities(self):

		# Initialize with the event of antigen growth.
		self.probabilities[0] = self.beta*self.antigen_time_series[-1]

		# fill with Bcell activation and proliferation
		for i in range(1, self.n_linages+1):
			#Activation
			rho = self.antigen_time_series[-1]/self.N_A
			self.probabilities[(2*i)-1] = (rho/(rho+np.exp(self.master_Sequence_energy + self.Sequences[i-1].energy)))*(1-self.Sequences[i-1].active)
			#Proliferation
			self.probabilities[(2*i)] = self.nu*(self.linages_time_series[i-1,-1])*(self.Sequences[i-1].active)

	def gillespie_step(self):

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 1. Generate 2 random numbers uniformly distributed in (0,1)
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		r1 = np.random.rand()
		r2 = np.random.rand()

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 2. Calculate probabilities
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		self.calculate_probabilities()

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 3. Calculate alpha
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		probabilities_cumsum = np.cumsum(self.probabilities)
		alpha = sum(self.probabilities)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 4. Compute the time until the next event takes place
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		tau = (1/alpha)*np.log(float(1/r1))
		self.time_series = np.append(self.time_series, self.time_series[-1]+tau)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 5. Compute which event takes place
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		transitionIdx   = np.searchsorted(probabilities_cumsum,r2*alpha)
		transition_Type  = transitionIdx % 2
		transition_Agent  = int((transitionIdx+1)/2)

		#antigen profileration event
		if(transition_Agent == 0):
			self.antigen_time_series = np.append(self.antigen_time_series, self.antigen_time_series[-1]+100)
			self.linages_time_series = np.hstack((self.linages_time_series, self.linages_time_series[:,-1].reshape(self.n_linages, 1)))
			self.activation_time_series = np.hstack((self.activation_time_series, self.activation_time_series[:,-1].reshape(self.n_linages, 1)))
		else:
			#Activation event
			if(transition_Type==1):
				self.Sequences[transition_Agent-1].active = True
				self.antigen_time_series = np.append(self.antigen_time_series, self.antigen_time_series[-1])
				temp_array = np.copy(self.activation_time_series[:,-1]).reshape(self.n_linages, 1)
				temp_array[transition_Agent-1] =  1
				self.linages_time_series = np.hstack((self.linages_time_series, self.linages_time_series[:,-1].reshape(self.n_linages, 1)))
				self.activation_time_series = np.hstack((self.activation_time_series, temp_array))
				self.active_linages +=1
			#B cell prolifration event
			else:
				self.antigen_time_series = np.append(self.antigen_time_series, self.antigen_time_series[-1])
				temp_array = np.copy(self.linages_time_series[:,-1]).reshape(self.n_linages, 1)
				temp_array[transition_Agent-1] = temp_array[transition_Agent-1] + 1
				self.linages_time_series = np.hstack((self.linages_time_series, temp_array))
				self.activation_time_series = np.hstack((self.activation_time_series, self.activation_time_series[:,-1].reshape(self.n_linages, 1)))

	def Gillespie(self):

		while((self.time_series[-1] < self.T)):
		#while((self.antigen_time_series[-1]<9e2) and (self.active_linages < 15)):
			self.gillespie_step()

	def plot_antigen_time(self, ax):

		ax.plot(self.time_series, self.antigen_time_series/self.N_A, linewidth  = 4)
		ax.set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'Antigen $\rho$', fontsize = 20)
		ax.tick_params(labelsize = 20)
		#handles, labels = ax.get_legend_handles_labels()
		#ax.legend(np.concatenate(([Line2D([0], [0], color='tab:red', linewidth=4, linestyle='solid', ms = 8)],handles)),np.concatenate(([r'$n_b(r, \rho)$'],labels)), loc = 0, fontsize = 20)

	def plot_prob_binding(self, ax):
		rho_array = np.logspace(0, np.log10(max(self.antigen_time_series)/self.N_A), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		for i, rho in enumerate(rho_array): 
			ax.plot(np.linspace(-6,8,10), (1/(1+np.exp(self.master_Sequence_energy + np.linspace(-6,8,10) - np.log(rho)))), linewidth  = 4, color = colors[i], label = r'$\rho = %.0e$'%(rho))
		ax.set_yscale('log')
		ax.set_xlabel(r'Energy $\epsilon$', fontsize = 20)
		ax.set_ylabel(r'Probability of binding $p_b$', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def stackplot_linages_time(self, ax, antigen = False, time = True):
		colors = []
		for i in self.Sequences:
			if(i.active==True):
				colors.append('tab:red')
			else:
				colors.append('indigo')
		if(time):
			ax.stackplot(self.time_series, self.linages_time_series/np.sum(self.linages_time_series, axis = 0), colors=colors, alpha = 0.9);
			ax.set_xlabel(r'Time $t$', fontsize = 20)
		if(antigen):
			ax.stackplot(self.antigen_time_series, self.linages_time_series/np.sum(self.linages_time_series, axis = 0), colors=colors, alpha = 0.9);
			ax.set_xlabel(r'Antigen $\rho$', fontsize = 20)
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'B cell Linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		#ax.legend(loc = 0, fontsize = 20)

	def hist_sequences_hamming_distance(self, Sequences, ax):

		rho_array = np.logspace(0, np.log10(max(self.antigen_time_series)/self.N_A), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		data_distances = ax.hist([Sequences[i].hamming_distance for i in range(int(len(Sequences)))], bins = range(10), align = 'left', label = r'$S(d)$', color = 'lightsteelblue', alpha = 0.5)
		ax.plot(data_distances[1][0:-1], sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), linewidth = 4 , color = 'lightsteelblue', alpha = 0.6)

		ax.hist([self.Sequences[i].hamming_distance for i in range(int(len(self.Sequences)))], bins = range(10), align = 'left', label = r'$US(d)$', color = 'indigo', alpha = 0.6)
		ax.plot(data_distances[1][0:-1], self.U*sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), linewidth = 4 , color = 'indigo', alpha = 0.6)

		ax.hist([self.Sequences[i].hamming_distance for i in range(int(len(self.Sequences))) if self.Sequences[i].active], bins = range(10), align = 'left', label = r'Activated Linages', color = 'tab:red', alpha = 0.8)
		for i, rho in enumerate(rho_array):
			ax.plot(data_distances[1][0:-1], self.U*sc.comb(9, data_distances[1][0:-1].astype(int))*((20-1)**data_distances[1][0:-1])*(1/(1+np.exp(self.master_Sequence_energy + data_distances[1][0:-1] - np.log(rho)))) , color = colors[i], linestyle = 'dashed', linewidth = 3)

		ax.set_ylim(0.1, 2e5)    
		ax.set_yscale('log')
		ax.set_xlabel(r'Hamming Distance $d$', fontsize = 20)
		ax.set_ylabel(r'Number of linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def hist_sequences_energy(self, Sequences, bins, ax):

		rho_array = np.logspace(0, np.log10(max(self.antigen_time_series)/self.N_A), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		data_energies = ax.hist([Sequences[i].energy for i in range(int(len(Sequences)))], bins = bins, align = 'left', label = r'$S(\epsilon)$', color = 'lightsteelblue', alpha = 0.5)
		#ax.plot(data_energies[1][0:-1], sc.comb(9, data_energies[1][0:-1])*((20-1)**data_energies[1][0:-1]), linewidth = 4 , color = 'lightsteelblue', alpha = 0.6)

		ax.hist([self.Sequences[i].energy for i in range(int(len(self.Sequences)))], bins = bins, align = 'left', label = r'$US(\epsilon)$', color = 'indigo', alpha = 0.6)
		#ax.plot(data_energies[1][0:-1], self.U*sc.comb(9, data_energies[1][0:-1])*((20-1)**data_energies[1][0:-1]), linewidth = 4 , color = 'indigo', alpha = 0.6)

		ax.hist([self.Sequences[i].energy for i in range(int(len(self.Sequences))) if self.Sequences[i].active], bins = bins, align = 'left', label = r'Activated Linages', color = 'tab:red', alpha = 0.8)
		#for i, rho in enumerate(rho_array):
		#	ax.plot(data_energies[1][0:-1], self.U*sc.comb(9, data_energies[1][0:-1].astype(int))*((20-1)**data_energies[1][0:-1])*(1/(1+np.exp(self.master_Sequence_energy + data_energies[1][0:-1] - np.log(rho)))) , color = colors[i], linestyle = 'dashed', linewidth = 3)

		#ax.set_ylim(0.1, 2e5)    
		ax.set_yscale('log')
		ax.set_xlabel(r'Energy $\epsilon$', fontsize = 20)
		ax.set_ylabel(r'Number of linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def plot_k_largest_linages(self, k, ax):

		#Calculate array of the frequencies of the largest k linages
		Seq_states = [i.active for i in self.Sequences]
		Seq_sizes = self.linages_time_series[:,-1]
		k = k
		biggest_k_linages_sizes = np.sort(Seq_sizes)[-k:]
		Pos = np.array([i for i, j in enumerate(Seq_sizes) if np.isin(j,biggest_k_linages_sizes)])
		biggest_k_linages = self.linages_time_series[Pos,:]
		#for i in range(1,int(len(self.linages_time_series[0,:]))):
		#    biggest_k_linages = np.vstack((biggest_k_linages, self.linages_time_series[Pos,i]))
		#biggest_k_linages_freq = np.transpose(biggest_k_linages)/np.sum(np.transpose(biggest_k_linages), axis = 0)
		biggest_k_linages_freq = biggest_k_linages/np.sum(biggest_k_linages, axis = 0)

		ax.stackplot(self.time_series, biggest_k_linages_freq, alpha = 0.9);
		#ax[0].set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'k largest linages', fontsize = 20)
		ax.tick_params(labelsize = 20)

		return biggest_k_linages_freq

	def plot_entropy_k_largest_linages(self, k, biggest_k_linages_freq, ax):

		#Calculate entropy
		entropy = [np.sum(-1*biggest_k_linages_freq[:,t]*np.log(biggest_k_linages_freq[:,t])) for t in range(int(len(self.time_series)))]
		ax.plot(self.time_series, entropy, linewidth = '4', color = 'indigo')
		#ax[1].set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'Entropy', fontsize = 20)
		ax.tick_params(labelsize = 20)


#----------------- Functions -----------------

def hamming_distance(chaine1, chaine2):

    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_complementary_seq(sequence, Energy_Matrix):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],min(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def find_complementary_seq_2(sequence, Energy_Matrix):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],max(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def calculate_energy(Energy_Matrix, seq1, seq2):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq1 = list(seq1)
	list_seq2 = list(seq2)
	Energy = 0
	for i in range(9):
		pos_i = np.where(np.isin(Alphabet,list_seq1[i]))[0][0]
		pos_j = np.where(np.isin(Alphabet,list_seq2[i]))[0][0]
		Energy += M[pos_i][pos_j]
	return Energy

def my_linear_func(x, a, b):

    return a + b*x

def my_quadratic_func(x, a, b, c):

    return np.log(a)+np.log(np.sqrt(-b)) + b*(x-c)**2

def my_plot_layout(ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 24, xlabel = '', ylabel = '', title = '', x_fontsize=24, y_fontsize = 24, t_fontsize = 24):
    ax.tick_params(labelsize = ticks_labelsize)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize = x_fontsize)
    ax.set_ylabel(ylabel, fontsize = y_fontsize)
    ax.set_title(title, fontsize = t_fontsize)


#----------------- Plots -----------------

def plot_energy_matrix(Energy_Matrix, Alphabet, title, ax):

	M = Energy_Matrix
	
	Alphabet = Alphabet

	sns.heatmap(np.flip(M, axis = 0), ax = ax, cmap=plt.cm.viridis, center = 0, cbar = True)
	ax.set_title(title, fontsize = 22)
	ax.tick_params(labelsize = 20)
	ax.set_xticklabels(Alphabet)
	ax.set_yticklabels(np.flip(Alphabet));
	cbar = ax.collections[0].colorbar
	cbar.ax.tick_params(labelsize=18)

def plot_histogram_hamming_distance(Sequences, ax):

	distances = np.array([i.hamming_distance for i in Sequences])
	data_distances = np.histogram(distances, bins=range(int(max(distances))))

	#ax.plot(data_distances[1][0:-1], scipy.special.comb(9, data_distances[1][0:-1]), linewidth = 4 , label = 'Binary')
	ax.plot(data_distances[1][0:-1], sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), color = 'steelblue', linewidth = 4 , label = '20-Alphabet')
	#ax.plot(data_distances[1][0:-1], scipy.special.comb(9, data_distances[1][0:-1])*((4-1)**data_distances[1][0:-1]), linewidth = 4 , label = '4-Alphabet')

	#ax.plot(data_distances[1][0:-1], np.exp(4*data_distances[1][0:-1]), linewidth = 4, label = r'$e^{\lambda r}$')
	ax.plot(data_distances[1][0:-1], data_distances[0], linewidth = 4, label = 'Data', linestyle = '', marker = 'o')

	ax.set_yscale('log')
	#ax.set_ylim(1,1e8)
	ax.set_xlabel(r'Hamming Distance $r$', fontsize = 20)
	ax.set_ylabel(r'', fontsize = 20)
	ax.tick_params(labelsize = 20)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(np.concatenate(([],handles)),np.concatenate(([],labels)), loc = 2, fontsize = 20)

	return distances

def plot_histogram_energy(Sequences, normalization, bins, color, density, n_seq, ax):

	energies = np.array([i.energy for i in Sequences])
	data_energies = np.histogram(energies, bins=bins, density = density)

	ax.plot(data_energies[1][0:-1], data_energies[0]/normalization, linewidth = 4, color = color, linestyle = '', marker = 'o', label = 'N=%.e'%(n_seq))
	ax.set_yscale('log')
	#ax.set_ylim(1,1e10)
	ax.set_xlabel(r'Energy $r$', fontsize = 20)
	ax.set_ylabel(r'', fontsize = 20)
	ax.tick_params(labelsize = 20)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(np.concatenate(([],handles)),np.concatenate(([],labels)), loc = 0, fontsize = 20)

	return energies, data_energies

def plot_scatter_hamming_distance_energy(distances, energies, color, ax):

	ax.scatter(distances, energies, color = color, s = 15, marker = '^')
	ax.hlines(energies[0],0,9, linestyle = 'dashed', label = r'Master Seq. energy $\epsilon_0$')
	#ax.set_yscale('log')
	#ax.set_ylim(1,1e10)
	ax.set_xlabel(r'Hamming distance $d$', fontsize = 20)
	ax.set_ylabel(r'Energy $\epsilon$', fontsize = 20)
	ax.tick_params(labelsize = 20)
	ax.legend(loc = 0, fontsize = 20)

def plot_histogram_energy_subsampling(bins, n_linages, Sequences, sub_energies, fig, ax):

	energies, data_energies = plot_histogram_energy(Sequences = Sequences, bins = bins, ax = ax)

	#____________ fit and plot gaussian function
	popt, pcov = curve_fit(my_quadratic_func, data_energies[1][0:-1][data_energies[0]!=0], np.log(data_energies[0][data_energies[0]!=0]) , p0 = (5e5/np.sqrt(np.pi), -0.04, -30))
	print(r'integral a*$\pi$:',popt[0]*np.pi)
	r_array = np.linspace(np.min(energies)-5, np.max(energies)+5, 1000)
	ax.plot(r_array, np.exp(my_quadratic_func(r_array, *popt)), linestyle = '--', linewidth = 4, color = 'steelblue', alpha = 0.4)
	r_array2 = np.linspace(-80, 50, 10000)


	#____________ Plot histogram of sub_energies
	data_energies = np.histogram(sub_energies, bins=20, density = False)
	ax.plot(data_energies[1][0:-1], data_energies[0]/(2000), linewidth = 4, color = 'indigo', label = 'After samplig', linestyle = '', marker = 'o')
	ax.plot(r_array, (2e2/5e5)*np.exp(my_quadratic_func(r_array, *popt)), linestyle = '--', linewidth = 4, color = 'indigo', alpha = 0.4)
	ax.plot(r_array, (2e2/5e5)*popt[0]*np.sqrt(-popt[1])*(1+popt[1]*(r_array-popt[2])**2) )
	ax.set_ylim(1,8e4)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(np.concatenate(([],handles)),np.concatenate(([],labels)), loc = 2, fontsize = 20)

	#____________ Create inset with the integral of the gaussian function  
	left, bottom, width, height = [0.65, 0.67, 0.25, 0.2]
	ax2 = fig.add_axes([left, bottom, width, height])
	ax2.plot(r_array, np.cumsum((2e2/5e5)*np.exp(my_quadratic_func(r_array, *popt))*(np.max(energies)-np.min(energies))/1000), linewidth = 2, color = 'indigo')
	ax2.set_xlabel(r'Energy $r$', fontsize = 15)
	ax2.set_ylabel(r'', fontsize = 15)
	ax2.tick_params(labelsize = 15)
	ax2.yaxis.tick_right()
	ax2.set_yscale('log')

	return popt, pcov

#----------------- Generate files -----------------

		
def print_raw_file(Sequences, filename):

	file = open(filename, 'w+')
	for i in range(len(Sequences)):
		np.savetxt(file, np.array([Sequences[i].parent_id, Sequences[i].id]), fmt = '%d', delimiter = ' ', newline = ' ')
		file.write("\n")
	file.close()

def generate_newick_format(filename):

	file  = np.loadtxt(filename, dtype = int)
	n_f = '0()'
	print(file)


#----------------- Plots for ensemble averages -----------------

def plot_activation_rate_ensemble_deterministic(beta, b, nu, gamma, T, initial_time, eo, n_linages, n_left_tail, rho_min, rho_max, dt, popt, energies, comment, gaussian, exponential, ax):

	N_A = 6.02214076e23
	to = (eo+np.log(N_A))/beta
	#____________ Read and plot the activation of linages as function of antigen concetration
	t_new = np.linspace(initial_time, T, int((T-initial_time)/dt))
	activation_time_series = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System//Text_files/ensemble_deterministic_model_activation_time_series_"+comment+".pkl", "rb" ) )
	activation_time_series_var = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System//Text_files/ensemble_deterministic_model_activation_time_series_var_"+comment+".pkl", "rb" ) )
	ax.plot(np.exp(beta*t_new[1:][::50][np.where(activation_time_series[1:][::50]!=0)])/N_A, activation_time_series[1:][::50][np.where(activation_time_series[1:][::50]!=0)], linestyle = '', marker = '.', ms = 15, linewidth = 4, color = 'indigo', label = 'simulation')
	ax.fill_between(np.exp(beta*t_new[1:][::10][np.where(activation_time_series[1:][::10]!=0)])/N_A , activation_time_series[1:][::10][np.where(activation_time_series[1:][::10]!=0)] - np.sqrt(activation_time_series_var[1:][::10][np.where(activation_time_series[1:][::10]!=0)]), activation_time_series[1:][::10][np.where(activation_time_series[1:][::10]!=0)] + np.sqrt(activation_time_series_var[1:][::10][np.where(activation_time_series[1:][::10]!=0)]), linewidth = 4, color = 'indigo', alpha = 0.2)

	#____________ Plot the gaussian integral
	if(gaussian):
		r_array = np.linspace(np.min(energies), np.max(np.log(np.exp(t_new[1:][::10][np.where(activation_time_series[1:][::10]!=0)])/N_A)), 5000)
		ax.plot(np.exp(r_array), np.cumsum((2e2/5e5)*np.exp(my_quadratic_func(r_array, *popt))*(np.max(energies)-np.min(energies))/5000), linestyle = '--', ms = 15, linewidth = 4, color = 'violet', label = 'Gaussian integral')

	#____________ Plot the exponential integral
	if(exponential):
		rho_new = np.logspace(np.log10(rho_min), np.log10(rho_max), 50)
		alpha = n_linages/n_left_tail
		ax.plot(rho_new, (alpha/(b*beta))*(np.exp(-beta*b*to)*(rho_new*N_A)**(b)-1), linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = 'exponential model', color = 'indigo')


	ax.set_xlabel(r'Antigen concentration $[M]$', fontsize = 20)
	ax.set_ylabel(r'Activated linages', fontsize = 20)
	ax.tick_params(labelsize = 22)
	ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_ylim(0.1, max(activation_time_series[1:]*1.5))
	ax.legend(loc = 0, fontsize = 20)

def plot_size_distribution_ensemble_deterministic(beta, b, nu, gamma, T, eo, dt, n_bins, density, popt, comment, gaussian, exponential, ax):

	N_A = 6.02214076e23
	to = (eo+np.log(N_A))/beta

	#____________ Read and plot the distribution of clone sizes
	activated_linages_size = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_linage_sizes_"+comment+".pkl", "rb" ) )
	bins = np.logspace(0, np.log10(np.max(activated_linages_size)), n_bins)
	data_activated_linages_log = np.histogram(activated_linages_size, bins = bins, density = True)
	
	#Distribution
	ax.plot(data_activated_linages_log[1][:-1], 1 - np.cumsum(data_activated_linages_log[0]*np.diff(data_activated_linages_log[1])), marker = '.', ms = 15, linestyle = '', linewidth = 3, color = 'indigo', label = 'Simulation')
	
	n_array = np.logspace(0,np.log10(np.max(activated_linages_size)), 50)
	#____________ Plot the gaussian integral
	if(gaussian):
		ax.plot(n_array, 1-((len(activated_linages_size)/2e2)*(2e2/5e5)*np.exp(my_quadratic_func(np.log((np.exp(T)/N_A)/(n_array)), *popt))), linestyle = '--', ms = 20, linewidth = 4, color = 'violet', label = 'Gaussian model')

	#____________ Plot the exponential integral
	if(exponential):
		ax.plot(n_array, ((1/(beta*b))*(np.exp(beta*b*(T-to))*n_array**(-(b*beta)/nu)-1))/((((1/(beta*b))*(np.exp(beta*b*(T-to))*n_array**(-b*beta)-1)))[0]), linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = 'exponential model', color = 'indigo')

	ax.set_xlabel(r'Clone size $n_i$', fontsize = 20)
	ax.set_ylabel(r'counts', fontsize = 20)
	ax.tick_params(labelsize = 22)
	ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_xlim(.9,1000)
	ax.legend(loc = 0, fontsize = 20)

def plot_N_total_ensemble_deterministic(beta, b, nu, gamma, T, initial_time, eo, n_linages, n_left_tail, dt, popt, comment, gaussian, exponential, ax):

	N_A = 6.02214076e23
	to = (eo+np.log(N_A))/beta
	t_new = np.linspace(initial_time, T, int((T-initial_time)/dt))

	#____________ Read and plot the distribution of clone sizes
	N_total_sim = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_N_total_"+comment+".pkl", "rb" ) )

	#Total linage size
	ax.plot(t_new[::50], N_total_sim[::50], marker = '.', ms = 20, linestyle = '', linewidth = 3, color = 'indigo', label = 'Simulation')

	#____________ Plot the gaussian integral
	if(gaussian):
		N_total = np.array([np.cumsum(np.exp(t-np.linspace(0,t, 100))*(2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt))*(t/100))[-1] + 200 - np.cumsum((2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt))*(t/100))[-1] for t in t_new])
		ax.plot(t_new, N_total, linestyle = '--', linewidth = 4, color = 'violet', alpha = 0.4, label = 'gaussian model')

	#____________ Plot the exponential integral
	if(exponential):
		tau = np.exp(t_new-to)
		alpha = n_linages/n_left_tail
		N_total_exp = alpha*((tau**(beta*b)-tau**(nu))/(beta*b-nu)) + n_linages - (alpha/(beta*b))*(tau**(beta*b)-1)
		ax.plot(t_new, N_total_exp, linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = 'exponential model', color = 'indigo')

	ax.set_xlabel(r'Time $t$', fontsize = 20)
	ax.set_ylabel(r'size $N_{total}$', fontsize = 20)
	ax.tick_params(labelsize = 22)
	#ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_xlim(T/2,T)
	ax.legend(loc = 0, fontsize = 20)

def plot_entropy_ensemble_deterministic(beta, b, nu, gamma, T, initial_time, eo, n_linages, n_left_tail, dt, popt, comment, gaussian, exponential, ax):

	N_A = 6.02214076e23

	to = (eo+np.log(N_A))/beta
	t_new = np.linspace(initial_time, T, int((T-initial_time)/dt))
	#____________ Read and plot the distribution of clone sizes
	entropy = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_entropy_"+comment+".pkl", "rb" ) )

	#Total linage size
	ax.plot(t_new[::50], entropy[::50], marker = '.', ms = 20, linestyle = '', linewidth = 3, color = 'indigo', label = 'Simulation')

	#____________ Plot the gaussian integral
	if(gaussian):
		N_total = np.array([np.cumsum(np.exp(t-np.linspace(0,t, 100))*(2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt))*(t/100))[-1] + 200 - np.cumsum((2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt))*(t/100))[-1] for t in t_new])
		N_alpha_log_N_alpha = np.array([np.cumsum((np.exp(t-np.linspace(0,t, 100))*(t-np.linspace(0,t, 100))*(2e2/5e5)*np.exp(my_quadratic_func(np.log(np.exp(np.linspace(0,t, 100))/N_A), *popt)))*(t/100))[-1] for t in t_new])
		ax.plot(t_new, (-1/N_total)*(N_alpha_log_N_alpha) + np.log(N_total), linestyle = '--', linewidth = 4, color = 'violet', alpha = 0.4, label = 'gaussian model')

	#____________ Plot the exponential integral
	if(exponential):
		tau = np.exp(t_new-to)
		alpha = n_linages/n_left_tail
		N_total_exp = alpha*((tau**(beta*b)-tau**(nu))/(beta*b-nu)) + n_linages - (alpha/(beta*b))*(tau**(beta*b)-1)
		Entropy_exp = ((nu*alpha)/(N_total_exp*(beta*b-nu))) * (tau**(nu)*np.log(tau)  - ((tau**(beta*b)-tau**(nu))/(beta*b-nu))) + ((np.log(N_total_exp)*1)/(1))
		ax.plot(t_new, Entropy_exp, linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = 'exponential model', color = 'indigo')
		ax.hlines(1/(1-beta*b), initial_time, T, linewidth = 4, linestyle = 'dashed', alpha = 0.4, label = r'$S_{\infty}$')

	ax.set_xlabel(r'Time $t$', fontsize = 20)
	ax.set_ylabel(r'Entropy $S_{Linages}$', fontsize = 20)
	ax.tick_params(labelsize = 22)
	#ax.set_xscale('log')
	ax.set_yscale('log')
	#ax.set_xlim(T/2,T)
	ax.legend(loc = 0, fontsize = 20)

#----------------------------------------------------------------

def generate_Sequences(n_seq, Energy_Matrix, antigen_sequence, L, new_antigen = False):

	M = Energy_Matrix
	L = L

	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']

	antigen_sequence = antigen_sequence

	if (new_antigen):
		antigen_sequence = "".join(np.random.choice(Alphabet, L))

	print('Antigen Seq: ' + antigen_sequence + '\n')

	master_sequence = find_complementary_seq(sequence = antigen_sequence , Energy_Matrix =  M)
	print('Master Seq: ' + master_sequence + '\n')

	master_sequence_energy = calculate_energy(Energy_Matrix = M, seq1 = master_sequence, seq2 = antigen_sequence)

	print('Master Seq energy: ', master_sequence_energy, '\n')

	Master_Sequence = Sequence(seq_id = 0, parent = master_sequence, energy_parent = master_sequence_energy, Master_Seq=True, master_sequence = master_sequence, complementary_sequence = antigen_sequence, Energy_Matrix = M)
	Sequences = np.array([Master_Sequence])
	sequences = np.array([master_sequence])
	zero_date = datetime(2020, 1, 1)

	#----------UNCOMMENT TO WRITE OUTPUT FILES-----------
	#file = open('../../../../Dropbox/Research/Evolution_Immune_System/Text_files/file_parent_daughter2.txt', 'w+')
	#file_1 = open('../../../../Dropbox/Research/Evolution_Immune_System/Text_files/timedNodeFile2.txt', 'w+')
	#file_2 = open('../../../../Dropbox/Research/Evolution_Immune_System/Text_files/strainFile2.txt', 'w+')
	#np.savetxt(file_1, np.array(['Node', 'Parent', 'Time', 'SigClade']), fmt='%d')
	#file.write("\n")

	#np.savetxt(file, np.array([str(Master_Sequence.id)+'\t', str(Master_Sequence.parent_id)]), fmt='%s', delimiter='', newline='', header = 'node\t parent\n', comments='')
	#file.write("\n")
	#np.savetxt(file_1, np.array([str(Master_Sequence.id)+'\t', str(Master_Sequence.parent_id)+'\t', str(43830)+'\t' , str(0)]), fmt='%s', delimiter='', newline='', header = 'Node\tParent\tTime\tSigClade\n', comments='')
	#file_1.write("\n")

	n_seq = n_seq
	Energy = {'energy': {'0': Master_Sequence.energy}}
	Hamming = {'hamming': {'0': 0}}
	repeated = 0
	for i in range(1, n_seq):
	    succ = False 
	    while(succ == False):
	        parent = np.random.choice(Sequences)
	        new_seq = Sequence(seq_id = i, parent = parent.sequence, energy_parent = parent.energy,  parent_id = parent.id, master_sequence = master_sequence, complementary_sequence = antigen_sequence, Energy_Matrix = M)
	        if (np.isin(new_seq.sequence, sequences)):
	        	repeated +=1
	        #check if the new sequence is already in the tree. Here we can check for other conditions like that the energy is higher than the parent.
	        if not(np.isin(new_seq.sequence, sequences)):
	            parent.tree_position = 0    
	            Sequences = np.append(Sequences, new_seq)
	            sequences = np.append(sequences, new_seq.sequence)
	            Energy['energy'][str(new_seq.id)] = new_seq.energy
	            Hamming['hamming'][str(new_seq.id)] = new_seq.hamming_distance
	            succ = True

	    #----------UNCOMMENT TO WRITE OUTPUT FILES-----------
	    #np.savetxt(file, np.array([str(new_seq.id)+'\t', str(new_seq.parent_id)]), fmt='%s', delimiter='', newline='')
	    #file.write("\n")
	    #np.savetxt(file_1, np.array([str(new_seq.id)+'\t', str(new_seq.parent_id)+'\t', str(43830 + i)+'\t', str(0)]), fmt='%s', delimiter='', newline='')
	    #file_1.write("\n")

	print('Found %d repeatead sequences. \n'%(repeated))
	#----------UNCOMMENT TO WRITE OUTPUT FILES-----------
	#for i in range(1, n_seq):
	    #if(Sequences[i].tree_position==1):
	        #new_date = zero_date + timedelta(i)
	        #np.savetxt(file_2, np.array([str(Sequences[i].id)+'\t', 'A/Germany/'+str(i)+'/2020'+'\t', 'Germany'+'\t', 'EPI_ISL_'+str(i)+'\t', str(new_date.year)+'-'+str(new_date.month)+'-'+str(new_date.day)+'\t', '0']),fmt='%s', delimiter='', newline='')
	        #file_2.write("\n")
	#----------UNCOMMENT TO WRITE OUTPUT FILES-----------
	#with open("../../../../Dropbox/Research/Evolution_Immune_System/Text_files/energy.json", 'w') as of:
	#	json.dump(Energy, of, indent=True)
	#with open("../../../../Dropbox/Research/Evolution_Immune_System/Text_files/hamming.json", 'w') as of:
	#	json.dump(Hamming, of, indent=True)

	#file.close()
	#file_1.close()
	#file_2.close()

	return Sequences

def generate_Sequences_randomly(n_seq, Energy_Matrix, antigen_sequence, L, new_antigen = False):

	M = Energy_Matrix
	L = L
	n_seq = n_seq

	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']

	antigen_sequence = antigen_sequence

	if (new_antigen):
		antigen_sequence = "".join(np.random.choice(Alphabet, L))
	print('Antigen Seq: ' + antigen_sequence + '\n')

	master_sequence = find_complementary_seq(sequence = antigen_sequence , Energy_Matrix =  M)
	print('Master Seq: ' + master_sequence + '\n')

	master_sequence_energy = calculate_energy(Energy_Matrix = M, seq1 = master_sequence, seq2 = antigen_sequence)
	print('Master Seq energy: ', master_sequence_energy, '\n')

	Master_Sequence = Sequence(seq_id = 0, parent = master_sequence, energy_parent = master_sequence_energy, Master_Seq=True, master_sequence = master_sequence, complementary_sequence = antigen_sequence, Energy_Matrix = M)
	
	Sequences = np.array([Master_Sequence])

	#Create n_seq random sequences and calculate the energy with respect to the given antigen.
	for i in range(1, n_seq):

		new_sequence = "".join(np.random.choice(Alphabet, L))
		new_sequence_energy = calculate_energy(Energy_Matrix = M, seq1 = new_sequence, seq2 = antigen_sequence)
		new_Sequence = Sequence(seq_id = i, parent = new_sequence, energy_parent = new_sequence_energy,  parent_id = i, master_sequence = master_sequence, complementary_sequence = antigen_sequence, Energy_Matrix = M, Master_Seq = True)

		Sequences = np.append(Sequences, new_Sequence)

	return Sequences

def run_ensemble_linage_size_distribution(Sequences, n_linages, n_seq, nu, beta, T, master_Sequence_energy, n_sim, new = False):
	n_linages = n_linages
	n_seq = n_seq
	U = n_linages/n_seq
	nu = nu
	R=6
	beta = beta
	gamma = 1
	T = T
	master_Sequence_energy = master_Sequence_energy

	#_____ Choose one of the following____________________________________
	if(new):
		activated_linages_size_t = []
		final_antigen_concentration = []
	else:
		activated_linages_size_t = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/activated_linages_size_t.pkl", "rb" ) )
		final_antigen_concentration = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/final_antigen_concentration.pkl", "rb" ) )
	#_____________________________________________________________________

	for i in range(n_sim):
		if(i%int((n_sim/5))==0):
			print(i, '...')
		Sub_Sequences = np.random.choice(Sequences, n_linages)
		Model = Stochastic_simulation(Sequences = Sub_Sequences, n_linages=n_linages, T = T, U = U, gamma = gamma, nu = nu, R = R, beta = beta, master_Sequence_energy = master_Sequence_energy)
		Model.Gillespie()
		activated_linages_size_t = np.append(activated_linages_size_t, [Model.linages_time_series[i,-1] for i in range(n_linages) if Model.Sequences[i].active])
		final_antigen_concentration = np.append(final_antigen_concentration, Model.antigen_time_series[-1])

	pickle.dump(activated_linages_size_t, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/activated_linages_size_t.pkl", "wb" ) )
	pickle.dump(final_antigen_concentration, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/final_antigen_concentration.pkl", "wb" ) )

	print('Ensemble size:', len(activated_linages_size_t))

def run_ensemble_deterministic_model(Sequences, n_linages, n_seq, nu, beta, gamma, T, energy_translation, initial_time, dt, n_sim, comment = "", new = False):

	n_linages = n_linages
	n_seq = n_seq
	U = n_linages/n_seq
	nu = nu
	beta = beta
	gamma = gamma
	T = T
	initial_time = initial_time
	dt = dt
	energy_translation = energy_translation

	#_____ Choose one of the following____________________________________
	if(new):
		activated_linages_size = np.array([])
		m_bar = np.array([])
		activated_energies = np.array([])
	else:
		activated_linages_size = pickle.load( open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_linage_sizes_"+comment+".pkl", "rb" ) )
		m_bar = pickle.load( open("../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_m_bar_"+comment+".pkl", "rb" ) )
		activated_energies = pickle.load( open("../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_activated_energies_"+comment+".pkl", "rb" ) )
	#_____________________________________________________________________

	N_total = np.zeros(int((T-initial_time)/dt))
	activation_time_series = np.zeros(int((T-initial_time)/dt))
	activation_time_series_2 = np.zeros(int((T-initial_time)/dt))
	entropy = np.zeros(int((T-initial_time)/dt))


	for i in range(n_sim):
		if(i%int((n_sim/5))==0):
			print(i, '...')
		Sub_Sequences = np.random.choice(Sequences, n_linages)
		for i in range(n_linages):
			Sub_Sequences[i].active = False
		Model = Deterministic_simulation(Sequences = Sub_Sequences, n_linages=n_linages, T = T, U = U, nu = nu, beta = beta, gamma = gamma, energy_translation = energy_translation, initial_time = initial_time, dt = dt)
		Model.ODE()
		# Add the linage sizes to the statistics
		activated_linages_size = np.append(activated_linages_size, [Model.linages_time_series[i,-1] for i in range(n_linages) if Model.Sequences[i].active])
		# Add m_bar to the statistics
		m_bar = np.append(m_bar, np.sum(Model.activation_time_series, axis = 0)[-1])
		# Sum up total population size
		N_total_i = np.sum(Model.linages_time_series, axis=0)
		# Calculate linage frequencies and entropy
		linage_freqs = Model.linages_time_series/N_total_i
		entropy_i = np.array([np.sum(-1*linage_freqs[:,t]*np.log(linage_freqs[:,t])) for t in range(int(len(Model.time_series)))])
		# Sum up columns in array Model.activation_time_series for the activated linages
		activation_time_series_i = np.array([np.sum(Model.activation_time_series[:,i]) for i in range(int(len(Model.activation_time_series[0,:])))])
		# Add the energies of the activated linages
		activated_energies = np.append(activated_energies, [i.energy for i in Model.Sequences[np.where(np.sum(Model.activation_time_series, axis=1)!=0)]])
		# Sum the last arrays to the main arrays
		activation_time_series = activation_time_series + activation_time_series_i
		activation_time_series_2 = activation_time_series_2 + (activation_time_series_i)**2
		N_total = N_total + N_total_i
		entropy = entropy + entropy_i


	activation_time_series = activation_time_series/(n_sim)
	activation_time_series_2 = (activation_time_series_2)/(n_sim)
	activation_time_series_var = activation_time_series_2-activation_time_series**2
	N_total = N_total/(n_sim)
	entropy = entropy/(n_sim)
    
	pickle.dump(activated_linages_size, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_linage_sizes_"+comment+".pkl", "wb" ) )
	pickle.dump(activation_time_series, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_activation_time_series_"+comment+".pkl", "wb" ) )
	pickle.dump(activation_time_series_var, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_activation_time_series_var_"+comment+".pkl", "wb" ) )
	pickle.dump(N_total, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_N_total_"+comment+".pkl", "wb" ) )
	pickle.dump(entropy, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_entropy_"+comment+".pkl", "wb" ) )
	pickle.dump(m_bar, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_m_bar_"+comment+".pkl", "wb" ) )
	pickle.dump(activated_energies, open( "../../../../Dropbox/Research/Evolution_Immune_System/Text_files/ensemble_deterministic_model_activated_energies_"+comment+".pkl", "wb" ) )


	print('Ensemble size:', len(activated_linages_size))
