from funcs_mini import*

# ------- functions analysis -------

def run_essay_time(data, antigen_past, past_id, antigen_current, current_id, antigen_test, test_id, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, beta, time_array, sera):
    # Es_update = []
    # antigen_seq_test = from_aa_to_i(antigen_test, energy_model, '../../')
    # motif = get_motif(antigen_seq_test, energy_model, '../../')*1.2
    # E_ms = np.ones(N_epi)
    # for epi in range(N_epi):
    #     E_m = -3
    #     # Normalize motif
    #     for i in range(l):
    #         E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
    #         motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
    #     E_ms[epi] = E_m
    #     for seq in data['seq']:
    #         E = calculate_energy(motif, from_aa_to_i(seq, energy_model, '../../')) + E_m
    #         Es_update.append(E)
    # data['E'] = Es_update

    # Find an input directory 
    Z = np.zeros_like(time_array)
    counter = 0
    exp_Z = 0
    for i in range(N_ens): # Use N_ens normally. For the moment let's use 1
        data_i = data.loc[data['ens_id']==i]
        data_i.reset_index(drop=True, inplace=True)
        data_i = ensemble_of_expansions_time(data_i, N_ens, p, time_array, lamB, C, dT) 
        Z_i = np.zeros_like(time_array)

        for epi in range(N_epi):
            data_epi = data_i.loc[data_i['epi']==epi+1]
            data_epi.reset_index(drop=True, inplace=True)        
            
            # print(data_epi['n_t'][data_epi['n_t']==1])
            if(len(data_epi)>0):
                if past_id == 0:
                    min_time = data_epi.min(numeric_only=True)['t']
                    min_time_idx = data_epi.idxmin(numeric_only=True)['t']
                    sizes = np.stack(data_epi['n_t'])
                    sizes[sizes == 1] = 0
                    Z_epi = (sizes/np.exp(np.array(data_epi['E']))[:, np.newaxis]).sum(axis = 0)
                else:
                    data_epi = data_epi.loc[data_epi['m']==1]
                    data_epi.reset_index(drop=True, inplace=True)  
                    min_time = data_epi.min(numeric_only=True)['t']
                    min_time_idx = data_epi.idxmin(numeric_only=True)['t']
                    sizes = np.stack(data_epi['n_t'])
                    sizes[sizes == 1] = 0
                    Z_epi = (sizes/np.exp(np.array(data_epi['E']))[:, np.newaxis]).sum(axis = 0)
                Z_i += Z_epi
                Z_i_min = Z_i[Z_i>0][0]
                Z_i_max = Z_i[-1]

        try: # this is the fitting setup that properly fits the exponent for the 1st infection. It overestimates the exponent for the 2nd infection
            if past_id == 0:
                popt, pcov = curve_fit(my_linear_func, time_array[(time_array>min_time*1.01) & (Z_i<0.1*Z_i_max)], np.log(Z_i[(time_array>min_time*1.01) & (Z_i<0.1*Z_i_max)]))
            else:
                popt, pcov = curve_fit(my_linear_func, time_array[(time_array>min_time*1.01) & (Z_i<0.1*Z_i_max)], np.log(Z_i[(time_array>min_time*1.01) & (Z_i<0.1*Z_i_max)]))
            exp_Z_i = popt[1]
            exp_Z += exp_Z_i
            Z += Z_i
            counter+=1
        except:
            print('Exponent for Z could not be fitted' )

    Z/=counter
    exp_Z/=counter

    sera['past'].append(antigen_past)
    sera['past id'].append(past_id)
    sera['current'].append(antigen_current)
    sera['current id'].append(current_id)
    sera['test'].append(antigen_test)
    sera['test id'].append(test_id)
    sera['beta'].append(beta)
    sera['exp_Z'].append(exp_Z)
    sera['Z'].append(list(Z))

    return sera

def run_essay_recurrent_time(data, inf, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, beta, time_array, sera):

    Z = []
    lamZ = []
    Lact = []

    if inf > 0:
        data = data.loc[data['m']==1]
    for i in range(N_ens): # Use N_ens normally. For the moment let's use 1
        data_i = data.loc[data['ens_id']==i]
        data_i.reset_index(drop=True, inplace=True)
        data_i = ensemble_of_expansions_time(data_i, N_ens, p, time_array, lamB, C, dT) 
        logZ_i = np.zeros_like(time_array)
        lamZ_i = np.zeros_like(time_array[::80][:-1])

        for epi in range(N_epi):
            data_epi = data_i.loc[data_i['epi']==epi+1]
            data_epi.reset_index(drop=True, inplace=True)        
            
            if(len(data_epi)>0):
                # min_time = data_epi.min(numeric_only=True)['t']
                max_time = data_epi.max(numeric_only=True)['t'] 
                min_time_idx = data_epi.idxmin(numeric_only=True)['t']
                sizes = np.stack(data_epi['n_t'])
                sizes[sizes == 1] = 0
                activated = sizes.copy()
                activated[activated > 0] = 1
                # sizes = sizes[min_time_idx, :]
                Ks = np.exp(np.array(data_epi['E']))[:, np.newaxis]
                # Ks = Ks[min_time_idx]
                Z_epi = (sizes/Ks).sum(axis = 0)
                Lact_epi = (activated).sum(axis = 0)
                Z_epi[Z_epi==0] = 1/Ks[min_time_idx]
                # Z_epi[Z_epi==0] = 1/Ks
                logZ_i += np.log(Z_epi)
                Z_i_min = logZ_i[logZ_i>0][0]
                Z_i_max = logZ_i[-1]
                lamZ_i += np.diff(logZ_i[::80])/np.diff(time_array[::80])

        Z.append(list(Z_epi))
        lamZ.append(list(lamZ_i))
        Lact.append(list(Lact_epi))

    sera['inf'].append(inf)
    sera['beta'].append(beta)
    sera['Z'].append(list(Z))
    sera['lamZ'].append(list(lamZ))
    sera['Lact'].append(list(Lact))

    return sera

def run_essay_Ab_time(data, antigen_past, past_id, antigen_current, current_id, antigen_test, test_id, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, beta, time_array, sera):
    # Es_update = []
    # antigen_seq_test = from_aa_to_i(antigen_test, energy_model, '../../')
    # motif = get_motif(antigen_seq_test, energy_model, '../../')*1.2
    # E_ms = np.ones(N_epi)
    # for epi in range(N_epi):
    #     E_m = -3
    #     # Normalize motif
    #     for i in range(l):
    #         E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
    #         motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
    #     E_ms[epi] = E_m
    #     for seq in data['seq']:
    #         E = calculate_energy(motif, from_aa_to_i(seq, energy_model, '../../')) + E_m
    #         Es_update.append(E)
    # data['E'] = Es_update

    # Find an input directory 
    T_1 = np.zeros_like(time_array)
    T_2 = np.zeros_like(time_array)
    counter = 0
    exp_Z = 0
    for i in range(N_ens): # Use N_ens normally. For the moment let's use 1
        data_i = data.loc[data['ens_id']==i]
        data_i.reset_index(drop=True, inplace=True)
        data_i = ensemble_of_expansions_time(data_i, N_ens, p, time_array, lamB, C, dT) 
        T_1_i = np.zeros_like(time_array)
        T_2_i = np.zeros_like(time_array)

        for epi in range(N_epi):
            data_epi = data_i.loc[data_i['epi']==epi+1]
            data_epi.reset_index(drop=True, inplace=True)        
            
            # print(data_epi['n_t'][data_epi['n_t']==1])
            if(len(data_epi)>0):
                if past_id == 0:
                    min_time = data_epi.min(numeric_only=True)['t']
                    min_time_idx = data_epi.idxmin(numeric_only=True)['t']
                    sizes = np.stack(data_epi['n_t'])
                    sizes[sizes == 1] = 0
                    Z_epi = (sizes/np.exp(np.array(data_epi['E']))[:, np.newaxis]).sum(axis = 0)
                    t12 = time_array[Z_epi>Z_epi[-1]/2][0]
                    Ab_1 = np.cumsum(2000*60*60*24*sizes*dT, axis = 1)
                    Ab_2 = np.cumsum(2000*60*60*24*sizes*(1/(1+np.exp(-4.*(time_array - t12))))*dT, axis = 1)
                    T_1_epi = (Ab_1/np.exp(np.array(data_epi['E']))[:, np.newaxis]).sum(axis = 0)
                    T_2_epi = (Ab_2/np.exp(np.array(data_epi['E']))[:, np.newaxis]).sum(axis = 0)

                else:
                    data_epi = data_epi.loc[data_epi['m']==1]
                    data_epi.reset_index(drop=True, inplace=True)  
                    min_time = data_epi.min(numeric_only=True)['t']
                    min_time_idx = data_epi.idxmin(numeric_only=True)['t']
                    sizes = np.stack(data_epi['n_t'])
                    sizes[sizes == 1] = 0
                    Z_epi = (sizes/np.exp(np.array(data_epi['E']))[:, np.newaxis]).sum(axis = 0)
                    t12 = time_array[Z_epi>Z_epi[-1]/2][0]
                    Ab_1 = np.cumsum(2000*60*60*24*sizes*dT, axis = 1)
                    Ab_2 = np.cumsum(2000*60*60*24*sizes*(1/(1+np.exp(-4.*(time_array - t12))))*dT, axis = 1)
                    T_1_epi = (Ab_1/np.exp(np.array(data_epi['E']))[:, np.newaxis]).sum(axis = 0)
                    T_2_epi = (Ab_2/np.exp(np.array(data_epi['E']))[:, np.newaxis]).sum(axis = 0)

                T_1_i += T_1_epi
                T_1_i_min = T_1_i[T_1_i>0][0]
                T_1_i_max = T_1_i[-1]

                T_2_i += T_2_epi
                T_2_i_min = T_2_i[T_2_i>0][0]
                T_2_i_max = T_2_i[-1]

        T_1 += T_1_i
        T_2 += T_2_i
        counter+=1

                # try: # this is the fitting setup that properly fits the exponent for the 1st infection. It overestimates the exponent for the 2nd infection
                #     if past_id == 0:
                #         popt, pcov = curve_fit(my_linear_func, time_array[(time_array>min_time*1.01) & (Z_i<0.1*Z_i_max)], np.log(Z_i[(time_array>min_time*1.01) & (Z_i<0.1*Z_i_max)]))
                #     else:
                #         popt, pcov = curve_fit(my_linear_func, time_array[(time_array>min_time*1.01) & (Z_i<0.1*Z_i_max)], np.log(Z_i[(time_array>min_time*1.01) & (Z_i<0.1*Z_i_max)]))
                #     exp_Z_i = popt[1]
                #     exp_Z += exp_Z_i
                #     Z += Z_i
                #     counter+=1
                # except:
                #     print('Exponent for Z could not be fitted' )

    T_1/=counter
    T_2/=counter
    # exp_Z/=counter

    sera['past'].append(antigen_past)
    sera['past id'].append(past_id)
    sera['current'].append(antigen_current)
    sera['current id'].append(current_id)
    sera['test'].append(antigen_test)
    sera['test id'].append(test_id)
    sera['beta'].append(beta)
    # sera['exp_Z'].append(exp_Z)
    sera['T_1'].append(list(T_1))
    sera['T_2'].append(list(T_2))

    return sera

def run_essay(data, antigen_past, past_id, antigen_current, current_id, antigen_test, test_id, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, time_array, sera):
    Es_update = []
    antigen_seq_test = from_aa_to_i(antigen_test, energy_model, '../../')
    motif = get_motif(antigen_seq_test, energy_model, '../../')*1.2
    E_ms = np.ones(N_epi)
    for epi in range(N_epi):
        E_m = -3
        # Normalize motif
        for i in range(l):
            E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
            motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
        E_ms[epi] = E_m
        for seq in data['seq']:
            E = calculate_energy(motif, from_aa_to_i(seq, energy_model, '../../')) + E_m
            Es_update.append(E)
    data['E'] = Es_update

    # Find an input directory 
    # Z = np.zeros_like(time_array)
    log_Z = 0
    counter = 0
    # Z_naive = np.zeros_like(time_array)
    # Z_memory = np.zeros_like(time_array)
    log_Z_naive = 0
    log_Z_memory = 0
    counter_0 = 0
    counter_1 = 0
    for i in range(N_ens):
        data_i = data.loc[data['ens_id']==i]
        data_i.reset_index(drop=True, inplace=True)
        # data_i = ensemble_of_expansions_time(data_i, N_ens, p, time_array, lamB, C, dT) 
        Z_i = 0
        Z_i_0 = 0
        Z_i_1 = 0
        for epi in range(N_epi):
            data_epi = data_i.loc[data_i['epi']==epi+1]
            data_epi.reset_index(drop=True, inplace=True)

            data_epi_0 = data_epi.loc[data_epi['m']==0]
            data_epi_1 = data_epi.loc[data_epi['m']==1]

            if(len(data_epi)>0):
                Z_epi = (data_epi['n']/np.exp(data_epi['E']))[data_epi.index].sum(axis = 0)
                Z_i += Z_epi
                log_Z += np.log(Z_i)
                counter+=1
                if len(data_epi_0.index)>0:
                    counter_0+=1
                    Z_epi_0 = (data_epi['n']/np.exp(data_epi['E']))[data_epi_0.index].sum(axis = 0)
                    Z_i_0 = Z_epi_0
                    log_Z_naive += np.log(Z_i_0)
                if len(data_epi_1.index)>0:
                    counter_1+=1
                    Z_epi_1 = (data_epi['n']/np.exp(data_epi['E']))[data_epi_1.index].sum(axis = 0)
                    Z_i_1 = Z_epi_1
                    log_Z_memory += np.log(Z_i_1)
        

    Z = np.exp(log_Z/counter)
    Z_naive = np.exp(log_Z_naive/counter_0)
    Z_memory= np.exp(log_Z_memory/np.max([1,counter_1]))

    sera['past'].append(antigen_past)
    sera['past id'].append(past_id)
    sera['current'].append(antigen_current)
    sera['current id'].append(current_id)
    sera['test'].append(antigen_test)
    sera['test id'].append(test_id)
    sera['Z'].append(Z)
    sera['Z naive'].append(Z_naive)
    sera['Z memory'].append(Z_memory)

    return sera

def run_essay_homo(data, antigen_past, past_id, antigen_current, current_id, antigen_test, test_id, energy_model, N_epi, N_ens, l, p, lamA, lamB, C, dT, time_array, sera):
    Es_update = []
    antigen_seq_test = from_aa_to_i(antigen_test, energy_model, '../../')
    motif = get_motif(antigen_seq_test, energy_model, '../../')*1.2
    E_ms = np.ones(N_epi)
    for epi in range(N_epi):
        E_m = -3
        # Normalize motif
        for i in range(l):
            E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
            motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
        E_ms[epi] = E_m
        for seq in data['seq']:
            E = calculate_energy(motif, from_aa_to_i(seq, energy_model, '../../')) + E_m
            Es_update.append(E)
    data['E'] = Es_update

    log_Z = np.zeros_like(time_array)
    counter = 0

    for i in range(N_ens):
        data_i = data.loc[data['ens_id']==i]
        data_i.reset_index(drop=True, inplace=True)
        data_i = ensemble_of_expansions_time(data_i, N_ens, p, time_array, lamB, C, dT) 
        Z_i = 1

        for epi in range(N_epi):
            data_epi = data_i.loc[data_i['epi']==epi+1]
            data_epi.reset_index(drop=True, inplace=True)

            if(len(data_epi)>0):
                counter+=1
                Z_epi = (data_epi['n_t']/np.exp(data_epi['E']))[data_epi.index].sum(axis = 0)
                Z_i += Z_epi
                log_Z += np.log(Z_i)

    Z = np.exp(log_Z/counter)
    Z12 = Z[Z<Z[-1]/2][-1]
    t12 = time_array[Z<Z12][-1]
    sera['current'].append(antigen_current)
    sera['current id'].append(current_id)
    sera['test'].append(antigen_test)
    sera['test id'].append(test_id)
    sera['Z12'].append(Z12)
    sera['t12'].append(t12)

    return sera

def DDG_distributions(data, DDF_df, antigen_past, past_id, antigen_current, current_id, antigen_test, test_id, energy_model, N_epi, l):
    DEs_update = []
    # antigen_seq_test = from_aa_to_i(antigen_test, energy_model, '../../')
    motif = get_motif(antigen_test, energy_model, '../../')*1.2
    E_ms = np.ones(N_epi)
    for epi in range(N_epi):
        E_m = -3
        # Normalize motif
        for i in range(l):
            E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
            motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
        E_ms[epi] = E_m
        for j, seq in enumerate(DDF_df['seq']):
            E = calculate_energy(motif, from_aa_to_i(seq, energy_model, '../../')) + E_m
            DEs_update.append(E - data[j])

    DDF_df[test_id] = DEs_update
    
    return DDF_df

def DDG(data, antigen_past, past_id, antigen_current, current_id, antigen_test, test_id, energy_model, N_epi, l, sera):
    Es_update = []
    antigen_seq_test = from_aa_to_i(antigen_test, energy_model, '../../')
    motif = get_motif(antigen_seq_test, energy_model, '../../')*1.2
    E_ms = np.ones(N_epi)
    for epi in range(N_epi):
        E_m = -3
        # Normalize motif
        for i in range(l):
            E_m+=np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
            motif[:, epi*l:(epi+1)*l][:, i] -= np.min(motif[:, epi*l:(epi+1)*l][:, i], axis=0)
        E_ms[epi] = E_m
        for seq in data['seq']:
            E = calculate_energy(motif, from_aa_to_i(seq, energy_model, '../../')) + E_m
            Es_update.append(E)
    # data['E'] = Es_update
    DDGs = [Es_update[j] - data['E'][j] for j in range(len(data.index))]
    
    sera['past'].append(antigen_past)
    sera['past id'].append(past_id)
    sera['current'].append(antigen_current)
    sera['current id'].append(current_id)
    sera['test'].append(antigen_test)
    sera['test id'].append(test_id)
    sera['DDG'].append(np.average(DDGs, weights = data['n']))
    sera['var_DDG'].append(np.var(DDGs, weights = data['n']))


    return sera




