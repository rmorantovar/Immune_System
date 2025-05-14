
from funcs_mini import*

# ------- functions -------

# Function to generate random sequences and compute properties

def generate_repertoire_Me(Alphabet, motif, Q0s, Ess, dEs, time_array, ensemble_id, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, infection, chunk_size, memory_clones, N_epi, DDE):
    lamA = lamA
    k_on = 1e6*24*3600; #(M*days)^-1
    b0 = 1e5
    b0_scaled = (b0 * 10 * k_on) / (lamA * N_A)
    K_step = k_step/k_on
    times = time_array
    exp_lamA_times = np.exp(lamA * times)
    properties = []
    for j in range(L0 // chunk_size):        
        for epi in range(N_epi):
            seqs_flat = np.random.rand(chunk_size)
            cum_Omega_0 = np.cumsum(Q0s[epi]*dEs[epi])[::]
            E_idx = np.searchsorted(cum_Omega_0,seqs_flat)-1
            Energies = Ess[epi][E_idx]
            Energies = Energies[Energies < np.min(Energies) + 4]
            factors = b0_scaled / (1 + (np.exp(Energies)/K_step))**p
            for i, factor in enumerate(factors):
                F1 = 1 - np.exp(-factor * (exp_lamA_times - 1))
                r1 = np.random.random()
                t1 = times[np.searchsorted(F1,r1)-1]
                if t1 < t_lim:
                    properties.append({
                        'ens_id': ensemble_id,
                        'E': Energies[i],
                        'id': seqs_flat[i],
                        't': t1,
                        'epi': epi+1,
                        'm' : 0
                    })

    if infection > 1:
        if len(memory_clones.index) > 0:
            memory_clones = memory_clones.loc[memory_clones['ens_id'] == ensemble_id]
            # sampled_memory = memory_clones.sample(n=50000, weights='N', replace=True)
            for index, row in memory_clones.iterrows():
                E = row['E'] + DDE
                N0 = row['N']
                epi = row['epi']
                id_clone = row['id']
                if E < E_lim:
                    F1 = 1-np.exp(-b0_scaled /(1+ (np.exp(E)/K_step))**pmem * (np.exp(lamA*times)-1)) # Here change p for 1 
                    # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                    r1 = np.random.random(int(N0))
                    #t1 = times[F1<r1][-1]
                    t1 = times[np.searchsorted(F1,r1)-1]
                    mask = t1 < t_lim
                    properties.extend([{
                                    'ens_id': ensemble_id,
                                    'E': E,
                                    'id': id_clone,
                                    't': t1_i,
                                    'epi': epi,
                                    'm': 1} for t1_i in t1[mask]])

    return properties

def generate_repertoire_Me_seqs(Alphabet, motif, cum_Omega_0, Es_avg, time_array, ensemble_id, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, infection, chunk_size, input_memory_file, N_epi):
    lamA = lamA
    k_on = 1e6*24*3600; #(M*days)^-1
    b0 = 1e5
    N0=1
    b0_scaled = (b0 * 10 * k_on) / (lamA * N_A)
    K_step = k_on / k_step
    times = time_array
    exp_lamA_times = np.exp(lamA * times)
    properties = []
    R = np.tile(np.arange(20), (int(chunk_size)*l, 1)).T
    for j in range(L0 // chunk_size):
        # print(j)
        seqs_flat = np.random.randint(0, 20, size=(int(chunk_size) * l))
        for epi in range(N_epi):
            Es = calculate_Es(motif[:, epi*l:(epi+1)*l], seqs_flat, R, l, 20, chunk_size, E_ms[epi])
            Es_idx = np.arange(int(chunk_size))[Es < np.min(Es) + 4]
            Es = Es[Es < np.min(Es) + 4]
            seqs = seqs_flat.reshape(int(chunk_size), l)[Es_idx]
            for i, E in enumerate(Es):
                factor = b0_scaled / (1 + (K_step * np.exp(E))**p)
                F1 = 1 - np.exp(-factor * (exp_lamA_times - 1))
                r1 = np.random.random()
                t1 = times[np.searchsorted(F1,r1)-1]
                if t1 < t_lim:
                    properties.append({
                        'ens_id': ensemble_id,
                        'E': E,
                        't': t1, 
                        # 'seq': from_i_to_aa_Alphabet(Alphabet, proto_E),
                        'seq': list(seqs[i]),
                        'epi': epi+1,
                        'm' : 0
                    })

    # for _ in range(L0 // chunk_size):
    #     proto_Es = np.random.randint(0, 20, size=(chunk_size, l))
    #     for proto_E in proto_Es:
    #         for epi in range(N_epi):
    #             E = calculate_energy(motif[:, epi*l:(epi+1)*l], proto_E) + E_ms[epi]
    #             # print(E)
    #             if E < E_lim:
    #                 F1 = 1-np.exp(-(b0*10*k_on)/(lamA*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lamA*times)-1))
    #                 r1 = np.random.random()
    #                 t1 = times[np.searchsorted(F1,r1)-1]
    #                 if t1 < t_lim:
    #                     properties.append({
    #                         'ens_id': ensemble_id,
    #                         'E': E,
    #                         't': t1, 
    #                         # 'seq': from_i_to_aa_Alphabet(Alphabet, proto_E),
    #                         'seq': list(proto_E),
    #                         'epi': epi+1,
    #                         'm' : 0
    #                     })

    if infection > 1:
        memory_clones = pd.read_csv(input_memory_file, converters={"seq": literal_eval})
        memory_clones = memory_clones.loc[memory_clones['ens_id'] == ensemble_id]
        if len(memory_clones.index) > 0:
            pre_memory_idx = []
            for index, row in memory_clones.iterrows():
                # seqs_memory += [row['seq']]*row['N']
                pre_memory_idx += [index]*int(row['N'])
            memory_idx = np.random.choice(pre_memory_idx, 1000, replace = False)
            memory = np.array(memory_clones['seq'].iloc[memory_idx])
            #memory = memory_clones['seq'].sample(n = 1000, replace = True, weights = np.array(memory_clones['n']))
            # print(len(seqs_memory))
            
            for seq in memory:
                # proto_E = from_aa_to_i_Alphabet(Alphabet, seq_aa)
                proto_E = seq
                for epi in range(N_epi):
                    E = calculate_energy(motif[:, epi*l:(epi+1)*l], proto_E) + E_ms[epi]
                    #for cell in range(int(clone['n'])):
                    if E < E_lim:
                        F1 = 1-np.exp(-(b0*10*k_on)/(lamA*N_A*(1+ (k_on*np.exp(E))/k_step)**p)*(np.exp(lamA*times)-1))
                        # for n0 in range(int(N0)): # to account for more than 1 cell per lineage
                        r1 = np.random.random()
                        #t1 = times[F1<r1][-1]
                        t1 = times[np.searchsorted(F1,r1)-1]
                        if t1 < t_lim:
                            properties.append({
                                'ens_id': ensemble_id,
                                'E': E,
                                't': t1,
                                'seq': list(seq),
                                'epi': epi+1,
                                'm' : 1
                            })


    return properties

# Define the merged and unified version of generate_repertoire_Me

def generate_repertoire_Me(Alphabet, motif, Q0s=None, Ess=None, dEs=None, time_array=None, ensemble_id=0, L0=1000, l=10, t_lim=5.0, E_lim=20.0, Es_ms=None, p=2, pmem=2, k_step=1.0, lamA=0.1, infection=1, chunk_size=100, memory_clones=None, N_epi=1, DDE=0.0, use_seqs=False):

    """
    Simulates the activation of immune cell repertoires under naive or memory conditions.

    This function supports two modes:
    - Energy-based: Directly samples energies from precomputed distributions (default).
    - Sequence-based: Samples sequences and computes energies from a motif matrix.

    Parameters
    ----------
    Alphabet : list
        List of amino acids or sequence characters (used in sequence reconstruction).
    motif : np.ndarray
        Motif matrix representing binding preferences per epitope.
    Q0s : list of np.ndarray, optional
        Base distributions used for energy sampling (used if `use_seqs=False`).
    Ess : list of np.ndarray, optional
        Energy values corresponding to Q0s (used if `use_seqs=False`).
    dEs : list of np.ndarray, optional
        Energy bin widths used with Q0s (used if `use_seqs=False`).
    time_array : np.ndarray
        Time grid over which activation is evaluated.
    ensemble_id : int, default=0
        ID of the simulated repertoire ensemble.
    L0 : int, default=1000
        Total number of cells to simulate.
    l : int, default=10
        Length of each sequence.
    t_lim : float, default=5.0
        Time threshold for successful activation.
    E_lim : float, default=20.0
        Maximum allowed energy for activation.
    Es_ms : list or np.ndarray
        Energy offsets to apply for each epitope.
    p : float, default=2
        Hill coefficient for naive clone activation.
    pmem : float, default=2
        Hill coefficient for memory clone activation.
    k_step : float, default=1.0
        Step affinity constant (controls nonlinearity).
    lamA : float, default=0.1
        Activation rate constant.
    infection : int, default=1
        Whether to include memory clones (>1 means include).
    chunk_size : int, default=100
        Number of clones processed per simulation chunk.
    memory_clones : pd.DataFrame, optional
        Preloaded memory clone data (required if `use_seqs=False` and `infection > 1`).
    N_epi : int, default=1
        Number of epitopes in the simulation.
    DDE : float, default=0.0
        Energy correction to apply to memory clones.
    use_seqs : bool, default=False
        Whether to use sequence-based simulation (True) or energy-based (False).

    Returns
    -------
    properties : list of dict
        A list of dictionaries, each representing an activated clone with:
        - 'ens_id': int, ensemble ID
        - 'E': float, energy
        - 't': float, activation time
        - 'epi': int, epitope index
        - 'm': int, memory flag (0 = naive, 1 = memory)
        - Optionally 'seq' or 'id' depending on mode
    """

    k_on = 1e6 * 24 * 3600
    b0 = 1e5
    N_A = 6.022e23
    b0_scaled = (b0 * 10 * k_on) / (lamA * N_A)
    K_step = k_step / k_on 
    times = time_array
    exp_lamA_times = np.exp(lamA * times)
    properties = []

    for j in range(L0 // chunk_size):
        if use_seqs:
            seqs_flat = np.random.randint(0, 20, size=(int(chunk_size) * l))
            R = np.tile(np.arange(20), (int(chunk_size)*l, 1)).T
        else:
            seqs_flat = np.random.rand(chunk_size)

        for epi in range(N_epi):
            if use_seqs:
                motif_epi = motif[:, epi*l:(epi+1)*l]
                Energies = calculate_Es(motif_epi, seqs_flat, R, l, 20, chunk_size, Es_ms[epi])
                Es_idx = np.arange(int(chunk_size))[Energies < np.min(Energies) + 4]
                Energies = Energies[Energies < np.min(Energies) + 4]
                factors = b0_scaled / (1 + (K_step * np.exp(Energies))**p)
                seqs = seqs_flat.reshape(int(chunk_size), l)[Es_idx]
                for i, factor in enumerate(factors):
                    F1 = 1 - np.exp(-factor * (exp_lamA_times - 1))
                    r1 = np.random.random()
                    t1 = times[np.searchsorted(F1, r1) - 1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': Energies[i],
                            'seq': list(seqs[i]),
                            't': t1,
                            'epi': epi + 1,
                            'm': 0
                        })
            else:
                cum_Omega = np.cumsum(Q0s[epi] * dEs[epi])[::]
                Es_idx = np.searchsorted(cum_Omega, seqs_flat) - 1
                Energies = Ess[epi][Es_idx]
                Energies = Energies[Energies < np.min(Energies) + 4]
                factors = b0_scaled / (1 + (np.exp(Energies) / K_step))**p
                for i, factor in enumerate(factors):
                    F1 = 1 - np.exp(-factor * (exp_lamA_times - 1))
                    r1 = np.random.random()
                    t1 = times[np.searchsorted(F1, r1) - 1]
                    if t1 < t_lim:
                        properties.append({
                            'ens_id': ensemble_id,
                            'E': Energies[i],
                            'seq': seqs_flat[i],
                            't': t1,
                            'epi': epi + 1,
                            'm': 0
                        })

    # Memory clones
    if infection > 1:
        if memory_clones is not None and len(memory_clones.index) > 0:
            memory_clones = memory_clones.loc[memory_clones['ens_id'] == ensemble_id]

            if use_seqs:
                pre_memory_idx = []
                for index, row in memory_clones.iterrows():
                    pre_memory_idx += [index] * int(row['N'])
                memory_idx = np.random.choice(pre_memory_idx, 1000, replace=False)
                memory = np.array(memory_clones['seq'].iloc[memory_idx])

                for seq in memory:
                    proto_E = seq
                    for epi in range(N_epi):
                        E = calculate_energy(motif[:, epi*l:(epi+1)*l], proto_E) + Es_ms[epi]
                        if E < E_lim:
                            F1 = 1 - np.exp(-(b0 * 10 * k_on) / (lamA * N_A * (1 + (k_on * np.exp(E)) / k_step)**p) * (np.exp(lamA * times) - 1))
                            r1 = np.random.random()
                            t1 = times[np.searchsorted(F1, r1) - 1]
                            if t1 < t_lim:
                                properties.append({
                                    'ens_id': ensemble_id,
                                    'E': E,
                                    't': t1,
                                    'seq': list(seq),
                                    'epi': epi + 1,
                                    'm': 1
                                })
            else:
                for index, row in memory_clones.iterrows():
                    E = row['E'] + DDE
                    N0 = row['N']
                    epi = row['epi']
                    id_clone = row['id']
                    if E < E_lim:
                        F1 = 1 - np.exp(-b0_scaled / (1 + (np.exp(E) / K_step))**pmem * (np.exp(lamA * times) - 1))
                        r1 = np.random.random(int(N0))
                        t1 = times[np.searchsorted(F1, r1) - 1]
                        mask = t1 < t_lim
                        properties.extend([{
                            'ens_id': ensemble_id,
                            'E': E,
                            'id': id_clone,
                            't': t1_i,
                            'epi': epi,
                            'm': 1
                        } for t1_i in t1[mask]])

    return properties

# -------

def expansions(data, time_array, p, lamB, C, dT):
    data_active = data.copy
    data_active = data.loc[data['t']<= np.min(data['t']) + 1/lamB*np.log(C/100)]
    clone_size_total = []
    clone_size_total_time = []
    ids_C_total = []
    lim_size = 2
    
    t_act_data = np.min(data_active['t'])
    activation_times = data_active['t'].values
    #---------------------------- B cell linages ----------------------
    clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lamB, C, dT)
    clone_sizes[clone_sizes == 1] = 0
    clone_size_total.extend(map(int, clone_sizes[:,-1]))
    clone_size_total_time.extend([list(clone_sizes[i, ::100]) for i in range(len(data_active))])

    # print(sum(clone_size_total))
    data_active = data_active.assign(N=clone_size_total)
    # data_active = data_active.assign(N_t=clone_size_total_time)
    data_active = data_active.loc[data_active['N']>=lim_size]
    return data_active

def expansionsN0(data, time_array, p, lamB, C, dT):
    data_active = data.copy
    data_active = data.loc[data['t']<= np.min(data['t']) + 1/lamB*np.log(C/10)]
    clone_size_total = []
    clone_size_total_time = []
    ids_C_total = []
    lim_size = 2
    
    t_act_data = np.min(data_active['t'])
    activation_times = data_active['t'].values
    N0s = data_active['N0'].values
    #---------------------------- B cell linages ----------------------
    clone_sizes = get_clones_sizes_C_new(len(activation_times), time_array, activation_times, N0s, lamB, C, dT)
    clone_sizes[clone_sizes == 1] = 0
    clone_size_total.extend(map(int, clone_sizes[:,-1]))
    clone_size_total_time.extend([list(clone_sizes[i, ::100]) for i in range(len(data_active))])

    # print(sum(clone_size_total))
    data_active = data_active.assign(N=clone_size_total)
    # data_active = data_active.assign(N_t=clone_size_total_time)
    data_active = data_active.loc[data_active['N']>=lim_size]
    return data_active

def ensemble_of_expansions_feedback(data, time_array, N_ens, p, lamB, C, dT):
    data = data.loc[data['t']<= np.min(data['t']) + 1/lamB*np.log(C/10)]
    clone_size_total = []
    clone_size_total_time = []
    ids_C_total = []
    lim_size = 2
    for i_ens in np.arange(N_ens):
        data_active = data.loc[data['ens_id']==i_ens]
        t_act_data = np.min(data_active['t'])
        activation_times = data_active['t'].values
        #---------------------------- B cell linages ----------------------
        clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lamB, C, dT)
        clone_sizes[clone_sizes == 1] = 0
        clone_size_total.extend(map(int, clone_sizes[:,-1]))
        clone_size_total_time.extend([list(clone_sizes[i, ::100]) for i in range(len(data_active))])

    print(sum(clone_size_total))
    data['N'] = clone_size_total
    data['N_t'] = clone_size_total_time
    data = data.loc[data['N']>=lim_size]
    return data

#used to get clone size time trajectories
def ensemble_of_expansions_time(data, N_ens, p, time_array, lambda_B, C, dT):

    data_active = data.loc[data['t']<= np.min(data['t']) + 1/lambda_B*np.log(C/10)]
    activation_times = np.array(data_active['t'])
    energies  = np.array(data_active['E'])
    ids = data_active.index.tolist()
    # N0s = data_active['N0'].values
    #---------------------------- B cell linages ----------------------
    clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lambda_B, C, dT)
    # print(clone_sizes)
    #--------------------------t_C filter-------------------------
    lim_size = 2
    clone_sizes_C, activation_times_C, energies_C, filter_C, n_C = apply_filter_C(clone_sizes, activation_times, energies, lim_size)
    ids_C = np.array(ids)[filter_C]
    clone_sizes_in_time = [clone_sizes_C[i] for i in range(len(ids_C))]
    data_active = data_active.iloc[ids_C]
    data_active.loc[:,'N_t'] = clone_sizes_in_time
    return data_active

def group_by_clones(df_expansion):
    df_grouped = pd.DataFrame([], columns = df_expansion.columns)
    for ens_id in range(np.max(df_expansion['ens_id'])):
        df_ens = df_expansion[df_expansion['ens_id']==ens_id]
        for epi in range(np.max(df_ens['epi'])):
            df_epi = df_ens[df_ens['ens_id']==epi]
            for seq in set(df_epi['seq']):
                df_seq = df_epi[df_epi['seq']==seq]
                avg_time = np.mean(df_seq['t'])
                clone_size = np.sum(df_seq['n'])
                memory = int(np.mean(df_seq['m']))
                E = np.mean(df_seq['E'])
                line = [ens_id, E, avg_time, seq, epi, memory, clone_size]
                line = pd.DataFrame(line, columns = df_expansion.columns)
                df_grouped = np.concatenate([df_grouped, line])
    return df_grouped

# -------

def response(Alphabet, motif, Q0s, Ess, dEs, time_array, dT, ensemble_id, L0, l, t_lim, E_lim, Es_ms, p, pmem, k_step, lamA, lamB, C, infection, chunk_size, memory_clones, N_epi, DDE, use_seqs=False):
    # props_activation = generate_repertoire_Me(Alphabet, motif, Q0s, Ess, dEs, time_array, ensemble_id, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, infection, chunk_size, memory_clones, N_epi, seqs, DDE)
    props_activation = generate_repertoire_Me(Alphabet, motif, Q0s=Q0s, Ess=Ess, dEs=dEs, time_array=time_array, ensemble_id=ensemble_id, L0=L0, l=l, t_lim=t_lim, E_lim=E_lim, Es_ms=Es_ms, p=p, pmem=pmem, k_step=k_step, lamA=lamA, infection=infection, chunk_size=chunk_size, memory_clones=memory_clones, N_epi=1, DDE=0.0, use_seqs=use_seqs) 
    df_props_activation = pd.DataFrame(props_activation)
    df_props_expansion = expansions(df_props_activation, time_array, p, lamB, C, dT)
    return df_props_expansion

def responseN0(Alphabet, motif, Q0s, Ess, dEs, time_array, dT, ensemble_id, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, lamB, C, infection, chunk_size, memory_clones, N_epi, DDE):
    props_activation = generate_repertoire_Me(Alphabet, motif, Q0s, Ess, dEs, time_array, ensemble_id, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, infection, chunk_size, memory_clones, N_epi, DDE)
    df_props_activation = pd.DataFrame(props_activation)
    df_props_expansion = expansionsN0(df_props_activation, time_array, p, lamB, C, dT)
    return df_props_expansion

def ensemble_of_responses(Alphabet, motif, Q0s, Ess, dEs, time_array, dT, N_ens, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, lamB, C, infection, chunk_size, input_memory_file, N_epi, DDE, use_seqs=False, n_jobs=-1):

    if input_memory_file != '':
        if use_seqs:
            memory_clones = pd.read_csv(input_memory_file, converters={"seq": literal_eval})
        else:
            memory_clones = pd.read_csv(input_memory_file)
    else:
        memory_clones = ''

    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        delayed(response)(Alphabet, motif, Q0s, Ess, dEs, time_array, dT, i, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, lamB, C, infection, chunk_size, memory_clones, N_epi, DDE, use_seqs=use_seqs) for i in range(N_ens)
    )
    df = pd.concat(results, ignore_index=True)
    return df

def ensemble_of_responsesN0(Alphabet, motif, Q0s, Ess, dEs, time_array, dT, N_ens, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, lamB, C, infection, chunk_size, input_memory_file, N_epi, DDE, n_jobs=-1):

    if input_memory_file != '':
        memory_clones = pd.read_csv(input_memory_file)
    else:
        memory_clones = ''

    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        delayed(responseN0)(Alphabet, motif, Q0s, Ess, dEs, time_array, dT, i, L0, l, t_lim, E_lim, E_ms, p, pmem, k_step, lamA, lamB, C, infection, chunk_size, memory_clones, N_epi, DDE) for i in range(N_ens)
    )
    df = pd.concat(results, ignore_index=True)
    return df

def response_seqs(Alphabet, motif, cum_Omega_0, Es_avg, time_array, dT, ensemble_id, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, lamB, C, infection, chunk_size, input_memory_file, N_epi):
    props_activation = generate_repertoire_Me_seqs(Alphabet, motif, cum_Omega_0, Es_avg, time_array, ensemble_id, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, infection, chunk_size, input_memory_file, N_epi)
    df_props_activation = pd.DataFrame(props_activation)
    df_props_expansion = expansions(df_props_activation, time_array, p, lamB, C, dT)
    return df_props_expansion

def ensemble_of_responses_seqs(Alphabet, motif, cum_Omega_0, Es_avg, time_array, dT, N_ens, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, lamB, C, infection, chunk_size, input_memory_file, N_epi, n_jobs=-1):

    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        delayed(response_seqs)(Alphabet, motif, cum_Omega_0, Es_avg, time_array, dT, i, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, lamB, C, infection, chunk_size, input_memory_file, N_epi) for i in range(N_ens)
    )
    df = pd.concat(results, ignore_index=True)
    return df

