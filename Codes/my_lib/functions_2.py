import sys
from functions_1 import*

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

# --- functions ----
def needleman_wunsch(x, y):
    """Run the Needleman-Wunsch algorithm on two sequences.

    x, y -- sequences.

    Code based on pseudocode in Section 3 of:

    Naveed, Tahir; Siddiqui, Imitaz Saeed; Ahmed, Shaftab.
    "Parallel Needleman-Wunsch Algorithm for Grid." n.d.
    https://upload.wikimedia.org/wikipedia/en/c/c4/ParallelNeedlemanAlgorithm.pdf
    """
    N, M = len(x), len(y)
    s = lambda a, b: int(a == b)

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1

    # Create tables F and Ptr
    F = {}
    Ptr = {}

    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j

    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] - 1,
            F[i, j - 1] - 1,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

    # Work backwards from (N - 1, M - 1) to (0, 0)
    # to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1

    return list(alignment)

def align_fast(x, y):
    """Align two sequences, maximizing the
    alignment score, using the Needleman-Wunsch
    algorithm.

    x, y -- sequences.
    """
    return needleman_wunsch(x, y)
    
delta_Nb = lambda t, tb, Nb, N, lambda_B, C: lambda_B*Nb*(1-(N/C))*np.heaviside(t-tb, 1)

#Quantity regime

t_E_0 = lambda E, theta, lambda_A, k_on, N_c, k_pr, N_r, L:-theta*E/lambda_A - np.log((lambda_A*N_A)/(k_on*N_c))/lambda_A + theta*np.log(k_pr/k_on)/lambda_A - np.log(N_r)/lambda_A + 1*L*np.log(20)*lambda_A

#Quality regime

E_t = lambda t, theta, lambda_A, k_on, N_c, k_pr:lambda_A*t/theta - np.log((lambda_A*N_A)/(k_on*N_c))/theta + np.log(k_pr/k_on) 

t_E = lambda E, theta, lambda_A, k_on, N_c, k_pr:theta*E/lambda_A + np.log((lambda_A*N_A)/(k_on*N_c))/lambda_A - theta*np.log(k_pr/k_on)/lambda_A


def deltaNb_func(Nb, t, n_act, lambda_B, C, activation_times):

	N  = np.sum(Nb)-n_act
	return lambda_B*Nb*(1-(N/C))*np.heaviside(t-activation_times, 1)

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = matplotlib.colors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def get_data_ensemble_times_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_times.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_times.pkl', 'rb')
		data = pickle.load(f)
		print('proccesed file has been read...')
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

def get_data_CSV(folder_path):
	return pd.read_csv(folder_path + 'output.csv')

def get_data_binary(folder_path, data_type):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_'+data_type+'.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		file = open(folder_path+'/processed_data_'+data_type+'.pkl', 'rb')
		data = pickle.load(file)
		print('Proccesed file has been read...')
		file.close()
		return data, return_data_type
	else:
		if os.path.exists(folder_path+'/energies_ensemble_active.bin'):
			with open(folder_path+'/energies_ensemble_active.bin', 'rb') as file:
				dataraw = np.fromfile(file, dtype=dt)
				data = pd.DataFrame(dataraw)
				print('.bin file has been read...')
			 	# Create an empty list to store the data
				#data = []
				# Loop until we reach the end of the file
				#while True:
				#	chunk = f.read(data_size)
				#	if not chunk:
				#		break
				#	data.append(struct.unpack("diidi", chunk))
			file.close()
			return data, return_data_type
		else:
			print(folder_path+'.bin file does not exist...')
			return [], 0

def get_data(folder_path, data_type):
    return_data_type = 0 #default for returning processed data
    # if os.path.exists(folder_path+'/processed_data_'+data_type+'.pkl'):
    try:
        file = open(folder_path+'/processed_data_'+data_type+'.pkl', 'rb')
        return_data_type = 1
        print('Data exists already and is proccesed.  Loading it ...')
        data = pickle.load(file)
        print('Proccesed file has been read...')
        file.close()
        return data, return_data_type
    # else:
    except FileNotFoundError:
        print(f'Processed data does not exist')
        return [], 0

def get_data_ensemble(folder_path):
	if os.path.exists(folder_path+'/energies_ensemble.pkl'):
		print('Object exists already, loading it ...')
		f = open(folder_path+'/energies_ensemble.pkl', 'rb')
		data = pickle.load(f) 
		return data
	else:
		print(f'Pickling data ...')
		data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
		f = open(folder_path+'/energies_ensemble.pkl', 'wb')
		pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
		return data

def get_data_ensemble_b(folder_path):
	#fmt = 'diidi'
	dt = np.dtype([('energy', 'd'), ('active', 'i4'), ('plasma', 'i4'), ('act_time', 'd'), ('i_ens', 'i4')])
	#data_size = struct.calcsize(fmt)
	with open(folder_path+'/energies_ensemble_active.bin', 'rb') as f:
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
	return data

def get_data_ensemble_times(folder_path):
	return_data_type = 0 #default for returning processed data
	if os.path.exists(folder_path+'/processed_data_times.pkl'):
		return_data_type = 1
		print('Data exists already and is proccesed.  Loading it ...')
		f = open(folder_path+'/processed_data_times.pkl', 'rb')
		data = pickle.load(f) 
	else:
		if os.path.exists(folder_path+'/energies_ensemble.pkl'):
			print('Data is not processed. Picke object exists already. Loading it ...')
			f = open(folder_path+'/energies_ensemble.pkl', 'rb')
			data = pickle.load(f) 
		else:
			print(f'Pickling data ...')
			data = pd.read_csv(folder_path+'/energies_ensemble.txt', sep = '\t', header=None)
			f = open(folder_path+'/energies_ensemble.pkl', 'wb')
			pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
			
	return data, return_data_type

def get_clones_sizes_C_old(n_act, time_array, activation_times, lambda_B, C, dT):
	clone_sizes = np.ones((n_act, len(time_array)))
	for i_t, t in enumerate(time_array[:-1]):
		for i in np.arange(0, n_act):
			tb = activation_times[i]
			Nb = clone_sizes[i, i_t]# * np.heaviside(tb-t)
			N = np.sum(clone_sizes[:, i_t]) - np.sum(clone_sizes[:, 0])
			clone_sizes[i, i_t+1] = Nb + delta_Nb(t, tb, Nb, N, lambda_B, C)*dT
	return clone_sizes

def get_clones_sizes_C(n_act, time_array, activation_times, lambda_B, C, dT):
	clone_sizes = np.ones((n_act, len(time_array)))
	for i_t, t in enumerate(time_array[:-1]):
		Nb = clone_sizes[:,i_t]
		N  = np.sum(Nb) - np.size(activation_times[activation_times>t])
		deltaNb = lambda_B*Nb*(1-(N/C))*np.heaviside(t-activation_times, 1)
		clone_sizes[:, i_t+1] = Nb + deltaNb*dT
		
	return clone_sizes

def get_clones_sizes_C_new(n_act, time_array, activation_times, lambda_B, C, dT):

	Nb0 = np.ones(n_act)

	return odeint(deltaNb_func, Nb0, time_array, args = (n_act, lambda_B, C, activation_times)).T

def get_motif(antigen, Matrix, Text_files_path):

	
	M = np.loadtxt(Text_files_path + Matrix + '.txt', skiprows= 0, usecols=range(0,20))
	Alphabet = np.loadtxt(Text_files_path+'/Alphabet_'+Matrix+'.txt', dtype=bytes, delimiter='\t').astype(str)
	Alphabet_list = Alphabet.tolist()

	antigen_list = [i for i in antigen]
	antigen_seq = np.array([], dtype = int)
	for i, aa in enumerate(antigen_list):
	    index = Alphabet_list.index(aa)
	    antigen_seq = np.append(antigen_seq, int(index))

	return M[:,antigen_seq], M, Alphabet

def get_motif_em(antigen, energy_model, M, Text_files_path, read_matrix = True):

	Alphabet = np.loadtxt(Text_files_path+'Input_files/Alphabet_'+energy_model+'.txt', dtype=bytes, delimiter='\t').astype(str)
	Alphabet_list = Alphabet.tolist()
	antigen_list = [i for i in antigen]

	if read_matrix:
		if energy_model == 'Gaussian':
			M = np.random.normal(0, 1, size=(20, 20))
		else:
			M = np.loadtxt(Text_files_path+'Input_files/' + energy_model + '.txt', skiprows= 0, usecols=range(0,20))
	else:
		M = M

	antigen_seq = np.array([], dtype = int)
	for i, aa in enumerate(antigen_list):
	    index = Alphabet_list.index(aa)
	    antigen_seq = np.append(antigen_seq, int(index))

	return M[:,antigen_seq], M, Alphabet

def get_motif_em_ids(antigen_ids, energy_model, M, Text_files_path, read_matrix = True):

	if read_matrix:
		if energy_model == 'Gaussian':
			M = np.random.normal(0, 1, size=(20, 20))
		else:
			M = np.loadtxt(Text_files_path+'Input_files/' + energy_model + '.txt', skiprows= 0, usecols=range(0,20))
	else:
		M = M

	return M[:,antigen_ids]

def get_repertoire_properties(betas, Q0, Es, dE, N_r):
	beta_r = betas[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
	E_r = Es[:-1][np.cumsum(Q0*dE)<(1/(N_r))][-1]
	Kd_r = np.exp(E_r)

	return beta_r, E_r, Kd_r

def get_proofreading_properties(betas, Q0, Es, dE, k_pr, k_on):
	E_pr = Es[:-1][Es[:-1]<np.log(k_pr/k_on)][-1]
	Kd_pr = np.exp(E_pr)
	beta_pr = betas[:-1][Es[:-1]<E_pr][-1]

	return beta_pr, E_pr, Kd_pr

def get_p_properties(betas, Q0, Es, dE, theta):

	beta_theta = betas[betas>theta][-1]
	E_theta = Es[betas>theta][-1]
	Kd_theta = np.exp(E_theta)

	return beta_theta, E_theta, Kd_theta

def hamming_distance(chaine1, chaine2):

    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_complementary_seq_min(sequence, Alphabet, Energy_Matrix):

	M = np.transpose(Energy_Matrix)
	#Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],min(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def find_complementary_seq_max(sequence, Alphabet, Energy_Matrix):

	M = Energy_Matrix
	#Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],max(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def calculate_energy(Alphabet, Energy_Matrix, seq1, seq2):

	M = Energy_Matrix
	#Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq1 = list(seq1)
	list_seq2 = list(seq2)
	Energy = 0
	for i in range(len(seq1)):
		pos_i = np.where(np.isin(Alphabet,list_seq1[i]))[0][0]
		pos_j = np.where(np.isin(Alphabet,list_seq2[i]))[0][0]
		Energy += M[pos_i][pos_j]
	return Energy

def calculate_energy_2(Energy_Matrix, seq):

	Energy = np.sum(Energy_Matrix[seq, range(int(len(seq)))])
	return Energy

def my_linear_func(x, a, b):

    return a + b*x

def my_quadratic_func(x, a, b, c):

    return np.log(a)+np.log(np.sqrt(-b)) + b*(x-c)**2

def my_log_logistic_func(x, a, b, c):

	#return (1.*1*np.exp(c*lambda_B*(x+a)**(b)))/(1.*1+(np.exp(c*lambda_B*(x+a)**(b))-1))
	return 1/(1+np.exp(-b*(x-c)**a))

def my_logistic_func(x, a, b, c):

	return np.log((((1.*C*np.exp(c*lambda_B*(x+a)**(b)))/(1.*C+(np.exp(c*lambda_B*(x+a)**(b))-1)))/Kd_r_renorm)/(C/Kd_r_renorm))

def binding_affinity(x, a, b):
	
    return 1/(1+((np.exp(a+b*x))))

def Z_PWM(PWM, T):
    Z = 1
    for i in range(len(PWM[0,:])):
        Z_i = 0
        for j in range(len(PWM[:,0])):
            Z_i = Z_i + np.exp((-PWM[j, i]/T))
        Z = Z*Z_i
    return Z

def Z_PWMj(PWM, beta):
    Z = 1
    for i in range(len(PWM[0,:])):
        Z_i = 0.j
        for j in range(len(PWM[:,0])):
            Z_i = Z_i + np.exp(-PWM[j, i]*beta)
        Z = Z*Z_i
    return Z

def Z_PWM_integral(T, lamda, k):
    return np.exp(-lamda*min_E)*(np.exp(max_E*(lamda-(1/T)))-np.exp(min_E*(lamda-(1/T))))/(lamda-(1/T))*k

def Z_PWM_integral2(T, lamda):
    return 2*np.exp(-lamda*min_E)*(np.exp(avg_E*(lamda-(1/T)))-np.exp(min_E*(lamda-(1/T))))/(lamda-(1/T))

def calculate_Q0(Tmin, Tmax, n_T, E_matrix, E_ms, L):

	Ts = np.linspace(Tmin, Tmax, n_T)
	betas = 1/Ts[:-1]
	F_PWM = -Ts*np.log(Z_PWM(E_matrix, Ts))
	Es = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts)) + E_ms
	dE = np.diff(Es)
	
	S = np.cumsum(betas[:-1]*dE)

	#Omega = 20**L
	Omega = np.sum(np.exp(S)*dE)
	#Omega = 2*np.sum(np.exp(S)*dE)
	
	#print('%.2e'%(20**L), '%.2e'%(np.sum(np.exp(S)*dE)))
	Q0 = np.exp(S)/Omega

	return Es, dE, Q0, betas

def calculate_Q0_2(Tmin, Tmax, n_T, E_matrix, E_ms, L):

	Ts = np.linspace(Tmin, Tmax, n_T)
	betas = 1/Ts[:-1]
	F_PWM = -Ts*np.log(Z_PWM(E_matrix, Ts))
	Es = F_PWM[:-1]-Ts[:-1]*(np.diff(F_PWM)/np.diff(Ts)) + E_ms
	dE = np.diff(Es)
	
	S1 = np.cumsum(betas[:-1]*dE)
	S2 = np.log(Z_PWM(E_matrix, Ts[:-2])) + betas[:-1]*(Es[:-1]-E_ms)
	S2 = -(np.diff(F_PWM)/np.diff(Ts))
	#Omega = 20**L
	Omega1 = np.sum(np.exp(S1)*dE)
	Omega2= np.sum(np.exp(S2[:-1])*dE)
	#Omega = 2*np.sum(np.exp(S)*dE)
	
	#print('%.2e'%(20**L), '%.2e'%(np.sum(np.exp(S)*dE)))
	Q01 = np.exp(S1)/Omega1
	Q02 = np.exp(S2)/Omega2

	return Es, dE, Q01, Q02, betas

def calculate_QR(Q0, k_on, k_act, rho_A, Es, q, lambda_A, N_c, dE):

	p_a = (1/(1 + (k_on*np.exp(Es[:-1])/k_act))**q)
	u_on = rho_A*k_on
	R = 1-np.exp(-u_on*p_a*N_c/lambda_A)
	QR = Q0*R

	return u_on, p_a, R, QR

def calculate_QR_t(Q0, k_on, k_act, E, rho_A_t, Es, q, lambda_A, N_c, dE):

	p_a = (1/(1 + (k_on*np.exp(Es[Es<E][-1])/k_act)**q) )
	u_on = rho_A_t*k_on
	R_t = 1-np.exp(-u_on*p_a*N_c/lambda_A)
	QR_t = R_t*Q0[Es[:-1]<E][-1]

	return u_on, p_a, R_t, QR_t

def get_t_act(time_array, N_r, Q0, k_on, k_pr, lambda_A, Es, dE, theta, N_c):
	u_on, p_a, P_act, Q_act = calculate_QR(Q0, k_on, k_pr, np.exp(lambda_A*time_array[0])/N_A, Es, theta, lambda_A, N_c, dE)
	m_bar = np.array([N_r*(1-np.sum(np.exp(-((p_a*(np.exp(lambda_A*t)/N_A*k_on*N_c))/lambda_A))*Q0*dE)) for t in time_array])
	t_act = time_array[m_bar>1][0] # Activation time	
	return t_act

def P_e_gaussian(avg_E, var_E, Es):

    return (2*np.pi*var_E)**(-0.5)*np.exp(-(Es-avg_E)**2/(2*var_E))

def P_min_e(N, avg_E, var_E, Es, dE):

    return (N*(1-np.cumsum(P_e_gaussian(avg_E, var_E, Es)*dE))**(N-1)*(P_e_gaussian(avg_E, var_E, Es)))

def P_min_e_Q0(N, Q0, dE):

    return (N*(1-np.cumsum(Q0*dE))**(N-1)*(Q0))

def apply_filter_C(clone_sizes, activation_times, energies, lim_size):
	filter_C = clone_sizes[:, -1] > lim_size
	n_C = np.sum(filter_C)
	clone_sizes_C = clone_sizes[filter_C, :]
	activation_times_C = activation_times[filter_C]
	energies_C = energies[filter_C]

	return clone_sizes_C, activation_times_C, energies_C, filter_C, n_C

def apply_filter_C_seqs(clone_sizes, activation_times, energies, seqs, lim_size):
	filter_C = clone_sizes[:, -1] > lim_size
	n_C = np.sum(filter_C)
	clone_sizes_C = clone_sizes[filter_C, :]
	activation_times_C = activation_times[filter_C]
	energies_C = energies[filter_C]
	seqs_C = seqs[filter_C]

	return clone_sizes_C, activation_times_C, energies_C, seqs_C, filter_C, n_C

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

	sns.heatmap(np.flip(M, axis = 0), ax = ax, cmap=plt.cm.seismic, center = 0, cbar = True)
	ax.set_title(title, fontsize = 22)
	ax.tick_params(labelsize = 20)
	ax.set_xticklabels(Alphabet)
	ax.set_yticklabels(np.flip(Alphabet));
	cbar = ax.collections[0].colorbar
	cbar.ax.tick_params(labelsize=18)

def plot_PWM(PWM, Alphabet, sequence, title, ax):

	M = PWM
	
	Alphabet = Alphabet

	sns.heatmap(np.flip(M, axis = 0), ax = ax, cmap=plt.cm.viridis, center = 0, cbar = True)
	ax.set_title(title, fontsize = 22)
	ax.tick_params(labelsize = 20)
	ax.set_xticklabels(sequence)
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
