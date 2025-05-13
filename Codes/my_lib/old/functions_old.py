import struct
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.special as sc
import time
from scipy.integrate import odeint

if(sys.version_info[1]<= 7):
    import pickle5 as pickle
else:
    import pickle

import json
import warnings

from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
from tqdm import tqdm
from scipy.optimize import curve_fit
from matplotlib import style
from matplotlib.collections import LineCollection

# ----- CONSTANTS ------
N_A = 6.02214076e23
k_BT = 1.380649e-23*293

my_red = np.array((228,75,41))/256.
my_purple = np.array((125,64,119))/256.
my_purple2 = np.array((116,97,164))/256.
my_green = np.array((125,165,38))/256.
my_blue = np.array((76,109,166))/256.
my_gold = np.array((215,139,45))/256.
my_brown = np.array((182,90,36))/256.
my_blue2 = np.array((80,141,188))/256.
my_yellow = np.array((246,181,56))/256.
my_green2 = np.array((158,248,72))/256.
my_cyan = 'tab:cyan'

antigen_color = my_yellow



def get_data_ensemble_K(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_K_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		if os.path.exists(folder_path+'/energies_ensemble.csv'):
			data = pd.read_csv(folder_path + '/energies_ensemble.csv')
		else:
			with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
				dataraw = np.fromfile(f, dtype=dt)
				print('bin file has been read...')
				data = pd.DataFrame(dataraw)
			 	# Create an empty list to store the data
				#data = []
				# Loop until we reach the end of the file
				#while True:
				#	chunk = f.read(data_size)
				#	if not chunk:
				#		break
				#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_K_mf(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_mf.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_mf.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_K_mf_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_mf.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_mf.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_S_mf(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_S_mf.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_S_mf.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_S_mf_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_S_mf.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_S_mf.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_K_largest(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_largest.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_largest.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_K_largest_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_largest.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_largest.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_K_elite(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_elite.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_elite.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None, engine = 'c')
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_K_elite_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_elite.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_elite.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_K_aging(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_aging.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_aging.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_K_aging_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_aging.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_aging.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_K_L(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_L.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_L.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None, engine = 'c')
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_K_L_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_K_L.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_K_L.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_ranking(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_ranking_size_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_size.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_size.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_ranking_2(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_2.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_2.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_ranking_affinity_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_affinity.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_affinity.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_ranking_3(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_3.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_3.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_ranking_potency_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_potency.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_potency.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_ranking_4(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_4.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_4.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_ranking_1_affinity_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_1_affinity.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_1_affinity.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type

def get_data_ensemble_ranking_combined(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_combined.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_combined.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f)
			print('pkl file has been read...')
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			print('txt file has been read...')
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_data_ensemble_ranking_combined_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_ranking_combined.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_ranking_combined.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
	else:
		with open(folder_path+'/energies_ensemble.bin', 'rb') as f:
			dataraw = np.fromfile(f, dtype=dt)
			data = pd.DataFrame(dataraw)
			print('bin file has been read...')
		 	# Create an empty list to store the data
			#data = []
			# Loop until we reach the end of the file
			#while True:
			#	chunk = f.read(data_size)
			#	if not chunk:
			#		break
			#	data.append(struct.unpack("diidi", chunk))

	return data, return_data_type