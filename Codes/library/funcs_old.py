
def ensemble_of_activations_Me(Alphabet, motif, cum_Omega_0, Es_avg, time_array, N_ens, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, infection, chunk_size, input_memory_file, N_epi, n_jobs=-1, seqs = False):

    results = Parallel(n_jobs=n_jobs, backend='loky', verbose=0)(
        delayed(generate_repertoire_Me_seqs)(Alphabet, motif, cum_Omega_0, Es_avg, time_array, i, L0, l, t_lim, E_lim, E_ms, p, k_step, lamA, infection, chunk_size, input_memory_file, N_epi) for i in range(N_ens)
    )
    all_properties = [prop for sublist in results for prop in sublist]
    df = pd.DataFrame(all_properties)
    return df
