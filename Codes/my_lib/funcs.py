
from funcs_mini import*

# ------- functions -------

# Function to generate random sequences and compute properties

def generate_repertoire_Me(
    Alphabet,
    motif,
    Q0s=None, Ess=None, dEs=None, Es_ms=None,
    time_array=None,
    ensemble_id=0, L0=1000, l=10,
    t_lim=5.0, E_lim=20.0,
    use_seqs=False,
    fixed_repertoire=None,
    **kwargs
):

    # Default values via kwargs (can be omitted in old code)
    p = kwargs.get('p', 4)
    pmem = kwargs.get('pmem', 2)
    k_step = kwargs.get('k_step', 1.0)
    lamA = kwargs.get('lamA', 0.1)
    infection = kwargs.get('infection', 1)
    chunk_size = kwargs.get('chunk_size', 100)
    memory_clones = kwargs.get('memory_clones', None)
    N_epi = kwargs.get('N_epi', 1)
    DDE = kwargs.get('DDE', 0.0)

    # Constants
    k_on = 1e6 * 24 * 3600
    b0 = 1e5
    N_A = 6.022e23
    b0_scaled = (b0 * 10 * k_on) / (lamA * N_A)
    K_step = k_step / k_on 
    times = time_array
    exp_lamA_times = np.exp(lamA * times)
    properties = []
    R = np.tile(np.arange(20), (int(chunk_size)*l, 1)).T

    for j in range(L0 // chunk_size):
        if fixed_repertoire is not None:
            seqs_flat = fixed_repertoire[j]
        else:
            if use_seqs:
                seqs_flat = np.random.randint(0, 20, size=(int(chunk_size) * l)) # This is the line where the repertoire is created
            else:
                seqs_flat = np.random.rand(chunk_size) # This is the line where the repertoire is created

        for epi in range(N_epi):
            if use_seqs:
                motif_epi = motif[:, epi*l:(epi+1)*l]
                Energies = calculate_Es(motif_epi, seqs_flat, R, l, 20, chunk_size, Es_ms[epi])
                min_energy = np.min(Energies)
                Es_idx = np.arange(int(chunk_size))[Energies < min_energy + 4]
                Energies = Energies[Energies < min_energy + 4]
                factors = b0_scaled / (1 + (np.exp(Energies) / K_step))**p
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
                min_energy = np.min(Energies)
                Energies = Energies[Energies < min_energy + 4]
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
                memory_idx = np.random.choice(pre_memory_idx, 100, replace=False)
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
                                    'seq': list(seq),
                                    't': t1,
                                    'epi': epi + 1,
                                    'm': 1
                                })
            else:
                weights = memory_clones['N'].values / memory_clones['N'].sum()
                chosen_indices = np.random.choice(memory_clones.index, size=200, replace=True, p=weights)
                selected_clones = memory_clones.loc[chosen_indices]
                for _, row in selected_clones.iterrows():
                    E = row['E'] + DDE
                    # N0 = row['N']
                    epi = row['epi']
                    seq = row['seq']
                    if E < E_lim:
                        F1 = 1 - np.exp(-b0_scaled / (1 + (np.exp(E) / K_step))**pmem * (np.exp(lamA * times) - 1))
                        r1 = np.random.random()
                        t1 = times[np.searchsorted(F1, r1) - 1]
                        if t1 < t_lim:
                            properties.append({
                            'ens_id': ensemble_id,
                            'E': E,
                            'seq': seq,
                            't': t1,
                            'epi': epi,
                            'm': 1
                        })

    return properties

# -------

def expansions(data, time_array, dT, **kwargs):
    """
    Simulates B-cell clonal expansions following activation.

    Parameters
    ----------
    data : pd.DataFrame
        Activation data containing clone activation times (`t`).
    time_array : np.ndarray
        Time grid used in the simulation.
    dT : float
        Time resolution for clone size evolution.
    kwargs : dict
        Additional parameters:
            - p : Hill coefficient (currently unused here)
            - lamB : Growth rate of clones
            - C : Carrying capacity

    Returns
    -------
    data_active : pd.DataFrame
        Activation data with estimated clone sizes (column 'N').
    """

    lamB = kwargs.get('lamB', 2)
    C = kwargs.get('C', 2e4)
    p = kwargs.get('p', 3)  # Currently unused but kept for compatibility/future use

    lim_size = 2
    data_active = data.copy()

    # Filter based on expansion time window
    t_cutoff = np.min(data['t']) + (1 / lamB) * np.log(C / 100)
    data_active = data_active.loc[data_active['t'] <= t_cutoff]
    min_energy = np.min(data_active.loc[data_active['epi']==2]['E'])
    # print(min_energy, '*')
    # print(data_active.loc[data_active['E']==min_energy]['t'], '*')
    # first_times = sorted(data_active['t'].unique())[:200]
    # print('max_time=', t_cutoff)
    # data_active = data_active.loc[data_active['t'].isin(first_times)]

    # print(data_active.loc[data_active['E'].isin([min_energy])])

    activation_times = data_active['t'].values

    # Simulate clone growth
    clone_sizes = get_clones_sizes_C(len(activation_times), time_array, activation_times, lamB, C, dT)
    clone_sizes[clone_sizes == 1] = 0

    clone_size_total = list(map(int, clone_sizes[:, -1]))
    # clone_size_total_time = [list(clone_sizes[i, ::100]) for i in range(len(data_active))]

    data_active = data_active.assign(N=clone_size_total)

    # print(data_active.loc[data_active['E']==min_energy]['t'])

    # data_active = data_active.assign(N_t=clone_size_total_time)  # Optionally store time-resolved sizes

    data_active = data_active.loc[data_active['N'] >= lim_size]
    # print(data_active.loc[data_active['E'].isin([min_energy])])
    return data_active

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

# -------

def response(
    Alphabet,
    motif,
    Q0s,
    Ess,
    dEs,
    time_array,
    dT,
    ensemble_id,
    L0,
    l,
    t_lim,
    E_lim,
    Es_ms,
    lamB,
    C,
    use_seqs=False,
    fixed_repertoire=None,
    **kwargs
):
    """
    Wrapper function to simulate immune response.

    Parameters
    ----------
    Alphabet, motif, Q0s, Ess, dEs, time_array, dT, etc.: core simulation inputs
    kwargs : dict
        Additional parameters passed to generate_repertoire_Me (e.g., p, pmem, k_step, lamA, infection, chunk_size, memory_clones, N_epi, DDE)
    """
    props_activation = generate_repertoire_Me(
        Alphabet,
        motif,
        Q0s=Q0s,
        Ess=Ess,
        dEs=dEs,
        time_array=time_array,
        ensemble_id=ensemble_id,
        L0=L0,
        l=l,
        t_lim=t_lim,
        E_lim=E_lim,
        Es_ms=Es_ms,
        use_seqs=use_seqs,
        fixed_repertoire=fixed_repertoire,
        **kwargs
    )

    df_props_activation = pd.DataFrame(props_activation)

    # Expansion still receives p, lamB, C as individual args
    # (you can also use kwargs.get() here if needed)
    df_props_expansion = expansions(df_props_activation, time_array, dT, **kwargs)

    return df_props_expansion

def ensemble_of_responses(
    Alphabet,
    motif,
    Q0s,
    Ess,
    dEs,
    time_array,
    dT,
    N_ens,
    L0,
    l,
    t_lim,
    E_lim,
    Es_ms,
    lamB,
    C,
    input_memory_file,
    use_seqs=False,
    n_jobs=-1,
    **kwargs
):
    """
    Simulates multiple immune repertoires in parallel, using `response`.

    Parameters
    ----------
    Alphabet, motif, Q0s, ... : core simulation inputs
    input_memory_file : str
        File containing memory clone data.
    use_seqs : bool
        Whether to run sequence-based simulation.
    n_jobs : int
        Number of parallel jobs.
    kwargs : dict
        Additional options passed down to `response` and `generate_repertoire_Me`.

    Returns
    -------
    df : pd.DataFrame
        Combined simulation results from all ensembles.
    """

    # Load memory clones
    if input_memory_file != '':
        if use_seqs:
            memory_clones = pd.read_csv(input_memory_file, converters={"seq": literal_eval})
        else:
            memory_clones = pd.read_csv(input_memory_file)
    else:
        memory_clones = None

    # Create fixed reperotire

    if kwargs.get('reuse_repertoire', False):
        chunk_size = kwargs.get('chunk_size', 100)
        fixed_repertoire = []
        for j in range(L0 // chunk_size):
            if use_seqs:
                fixed_repertoire.append(np.random.randint(0, 20, size=(chunk_size * l)))
            else:
                fixed_repertoire.append(np.random.rand(chunk_size))
        kwargs['fixed_repertoire'] = fixed_repertoire

    # Parallel execution of ensembles
    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        delayed(response)(
            Alphabet, motif, Q0s, Ess, dEs,
            time_array, dT, i, L0, l, t_lim, E_lim, Es_ms,
            lamB, C,
            use_seqs=use_seqs,
            memory_clones=memory_clones,
            **kwargs
        ) for i in range(N_ens)
    )

    df = pd.concat(results, ignore_index=True)
    return df
