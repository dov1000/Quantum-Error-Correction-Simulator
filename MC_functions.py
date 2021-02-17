import sys
import os
import math
import json
import error
import cross
import multiprocessing as mp
import chper_wrapper as wrapper
import circuit as c
import steane
import correction as cor
from visualizer import browser_vis as brow


gates = ['PrepareXPlus', 'PrepareXMinus', 'PrepareYPlus', 'PrepareYMinus', 
         'PrepareZPlus', 'PrepareZMinus', 'PrepareX', 'PrepareY', 'PrepareZ', 
         'I', 'X', 'Y', 'Z', 'H', 'S', 'T', 'CX', 'CZ', 'CS', 'CY']



timing_dict = {'PrepareXPlus':     50,
               'PrepareXMinus':    50,
               'PrepareYPlus':     50,
               'PrepareYMinus':    50,
               'PrepareZPlus':     50,
               'PrepareZMinus':    50,
               'PrepareX':         50,
               'PrepareY':         50,
               'PrepareZ':         50,
               'I':                10,
               'X':                10,
               'Y':                10,
               'Z':                10,
               'H':                10,
               'S':                10,
               'T':                10,
               'CX':               80,
               'CZ':               80,
               'CS':               80,
               'MeasureX':         500,
               'MeasureY':         500,
               'MeasureZ':         500,
               'MeasureXDestroy':  500,
               'MeasureYDestroy':  500,
               'MeasureZDestroy':  500
              }



def read_error_info(error_info_or_file):
    '''
    Reads in the error information for the Monte Carlo simulation.
    The input can be:
        (a) a float, in which case the noise is assumed to be depolarizing
            (symmetric Pauli).
        (b) a string, in which case it is assumed to be the name of the json
            file containing the error information.
        (c) a dictionary, in which case it is assumed to be the already the
            the dictionary with the error information.
    '''

    try:
        p = float(error_info_or_file)
        error_kind = 1     # error after every gate
        error_dict = wrapper.create_dep_noise_error_dict(gates, p, error_kind)

    except TypeError:

        # if error_info_or_file is a filename
        if type(error_info_or_file) == type(''):
            error_info = error.read_error(error_info_or_file, alternative)
           
        # if error_info_or_file is a dictionary
        elif type(error_info_or_file) == type({}):
            error_dict = error_info_or_file

        else:
            error_message = 'The error info needs to be either a float, '
            error_message += 'a filename pointing to a json file, or a '
            error_message += 'dictionary.'
            raise ValueError(error_message)

    if 'error_dict' in locals():
        if 'error_rates' in error_dict:
            rates = error_dict['error_rates']
            one_qubit = [(str(x),y) for x,y in error_dict['one_qubit_ratio'].items()]
            two_qubit = [(str(x),y) for x,y in error_dict['two_qubit_ratio'].items()]
            error_info = error.ErrorRates(rates, one_qubit, two_qubit, error_kind)
        else:
            error_info = error.ErrorRatesAlternative(error_dict)

    return error_info


def pre_run_MC(error_subsets, operation, error_info, n_proc, output_folder, n_runs_total,
               n_json_files, init_state, initial_I, initial_trans, sampling, 
               error_subsets_tol, Is_after_two_qubit):

    #print 'error subsets =', error_subsets
    #print 'operation =', operation
    #print 'error info =', error_info
    #print 'n_proc =', n_proc
    #print 'output folder =', output_folder
    #print 'n_runs_total =', n_runs_total
    #print 'n_json_files =', n_json_files
    #print 'sampling =', sampling
    #sys.exit(0)

    if sampling == 'old':
        run_MC(operation, error_info, n_proc, output_folder, n_runs_total,
               n_json_files, init_state, initial_trans)

    elif sampling == 'Muyalon':

        # We need to calculate n_runs_total as a function of the number of elements
        # in the set with the maximum number of elements or something like that.
        # Mauricio 2/2/2016.

        faulty_gates = error_info.dic.keys()
        gate_indices = wrapper.gates_list_for_operation(operation,
                                                        faulty_gates,
                                                        Is_after_two_qubit)

        # error_subsets = wrapper.find_subsets(p,p,len(gate_indices[0]), len(gate_indices[1]), 
        #                                       error_subsets_tol)
        print 'subset is', error_subsets

        # total_prob = wrapper.prob_for_subset(p,p,len(gate_indices[0]),len(gate_indices[1]),0, 0)
        # # total_prob = 0
        # for n_errors in error_subsets:
        #     total_prob += wrapper.prob_for_subset(p,p,len(gate_indices[0]),
        #                                            len(gate_indices[1]),n_errors[0], n_errors[1])

        # for n_errors in error_subsets:
            # runs = int(round(n_runs_total * (wrapper.prob_for_subset(p,p,len(gate_indices[0]), 
            #                                   len(gate_indices[1]), n_errors[0], n_errors[1])/total_prob)))
            # if runs >0:
            # print 'n errors =', n_errors
        specific_output_folder = output_folder + str(error_subsets) + '/'
        run_MC(operation, error_info, n_proc, specific_output_folder,
               n_runs_total, n_json_files, init_state, initial_I, initial_trans,
               sampling, error_subsets, gate_indices, Is_after_two_qubit)


def run_MC(operation, error_info, n_proc, output_folder, n_runs_total,
           n_json_files, init_state, initial_I, initial_trans,
           sampling, n_errors, gate_indices, Is_after_two_qubit):
    '''
    initial_trans refers to the logical transversal gate that the user
    would like to add before the EC circuit.  If False, no gate is added
    (just the Identity).
    '''   

    n_runs_partial = int(n_runs_total/n_json_files)
    n_runs_parallel = int(float(n_runs_partial)/n_proc)

    if output_folder[-1] != '/':  output_folder += '/'
    chp_loc = './chp_extended'
    CHP_IO_files = False
    add_error = True
    if operation == 'Cross' or operation == 'CrossEC':
        prep_stabs = wrapper.prepare_stabs_Cross
    elif operation == 'ShorEC_fivequbit':
        prep_stabs = wrapper.prepare_stabs_5qubit
    else:
        prep_stabs = wrapper.prepare_stabs_Steane
        
    init_state_stab, init_state_destab = prep_stabs(init_state,
                                                    chp_loc,
                                                    CHP_IO_files)

    print 'Total number of runs = %i' %int(n_runs_total)
    # print 'Running the Monte Carlo simulations ...'

    for json_index in range(n_json_files):

        # print 'Progress = %i/%i' %(json_index, n_json_files)

        json_filename = str(json_index) + '.json'
        abs_filename = output_folder + json_filename
        sim_func = wrapper.run_simulation_EC

        if n_proc == 1:
            final_dict = sim_func(n_runs_parallel, init_state_stab,
                                  init_state_destab, chp_loc,
                                  add_error, error_info, n_errors, gate_indices,
                                  CHP_IO_files, 0, operation, initial_I, initial_trans, 
                                  init_state, sampling, Is_after_two_qubit)

        else:

            ########## Parallel part ##########

            pool = mp.Pool(n_proc)
            results = [pool.apply_async(sim_func, (n_runs_parallel, init_state_stab,
                                  init_state_destab, chp_loc,
                                  add_error, error_info, n_errors, gate_indices, CHP_IO_files,
                                  i*n_runs_parallel, operation, initial_I, initial_trans, 
                                  init_state, sampling, Is_after_two_qubit))
                                  for i in range(n_proc)]
            pool.close()
            pool.join()

            dicts = [r.get() for r in results]

            l = []
            # fix the problem that parallel runs' results will overwrite each other
            #for i in range(len(dicts)):
            #    d = dicts[i]
            #    for j in range(len(l)):
            #        for key in d.keys():
            #            if key == l[j][0]:
            #                d[key + i*n_runs_parallel] = d.pop(key)
            #                break

            #    l += d.items()
            #final_dict = dict(l)

            for d in dicts:
                l += d.items()
            final_dict = dict(l)

            ###################################  

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        json_file = open(abs_filename, 'w')
        json.dump(final_dict, json_file, indent=4, separators=(',', ':'),
                  sort_keys=True)
        json_file.close()


    return None



def run_MC_until_failure(n_proc, chp_loc, decoder, Is_after2q, Is_after1q,
                         MS_heating, Stark, qubit_assignment, time_MS,
                         time_1q, output_folder, error_info, n_sim_total, init_state='+Z',
                         initial_I=False, initial_trans=False, ion_trap=False,
                         out_name='', corrtable=None, vec_path=None, init_vec=None,
                         gamma=None, delta=None, save_error_dict=False):

    if decoder[:3] == 'max':   
    
        #always forget whether to put the final slash for the path
        if vec_path[-1] != '/':
            vec_path = vec_path+'/'

        errorflag = 0
        #Test the inputs for the MLD
        if decoder[0:9] == 'max_like_' and (corrtable or vec_path or init_vec or gamma ) == None:
            print 'MLD chosen but correction table, vector path, initial vector, or gamma vector not input'
            errorflag = 1
        if decoder == 'max_like_onestep' and delta == None:
            print 'One step MLD chosen but delta vector not input.'
        #Further test if the inputs are the correct format
        if corrtable != None:
            try:
                with open(corrtable) as f:
                    json.load(f)
            except (ValueError, IOError):
                print 'Correction table file either does not exist or is not a *.json file.'
                errorflag = 1
        if vec_path != None and os.path.isdir(vec_path) == False:
            print 'MLD vector path is not a directory.'
            errorflag = 1
        if init_vec != None:
            try:
                np.load(vec_path+init_vec)
            except IOError:
                print 'Init vector file does not exist or is not a numpy vector.'
                errorflag = 1
        if gamma != None:
            try:
                np.load(vec_path+gamma)
            except IOError:
                print 'Gamma vector file does not exist or is not a numpy vector.'
                errorflag = 1
        if delta != None:
            try:
                np.load(vec_path+delta)
            except IOError:
                print 'Delta vector file does not exist or is not a numpy vector.'
                errorflag = 1
   

        #Attempt to load the correction dictionary if it exists
        if corrtable != None:
            try:
                with open(corrtable) as f:
                    synd2corr_dict = json.load(f)
            except IOError:
                pass

        #check consistency of gamma vector with correction table
        try:
            if gamma != None and init_vec != None:
                tgvec=np.load(vec_path+gamma)
                tistate=np.load(vec_path+init_vec)
                if len(synd2corr_dict.keys()[0]) != int(log(len(tgvec),2)) or len(synd2corr_dict.keys()[0]) != int(log(len(tistate),2)):
                    errorflag = 1
                    print 'The length of the gamma/initial vector is not consistent with syndrome lengths in the correction table.'
                del tgvec
                del tistate
        except (NameError, IOError):
            pass

        if errorflag == 1:
            sys.exit('Fatal error in decoder inputs!')

        #clean up inputs
        if init_vec == None:
            tivec = None
        else:
            tivec = vec_path+init_vec
        if gamma == None:
            tgam = None
        else:
            tgam = vec_path+gamma
        if delta == None:
            tdel = None
        else:
            tdel = vec_path+delta


    else:
        synd2corr_dict = None
        tivec = None
        tgam = None
        tdel = None
        
 
    '''
    initial_trans refers to the logical transversal gate that the user
    would like to add before the EC circuit.  If False, no gate is added
    (just the Identity).
    '''   

    if output_folder[-1] != '/':  output_folder += '/'
    n_json_files = 1
    n_sim = n_sim_total/n_proc
       
    print 'Running the Monte Carlo simulations ...'

    for json_index in range(n_json_files):

        # print 'Progress = %i/%i' %(json_index, n_json_files)
        sim_func = wrapper.run_several_sim_until_failure


        if n_proc == 1:

            if len(out_name) > 0:
                sim_func(chp_loc, decoder, Is_after2q, Is_after1q,
                         MS_heating, Stark, qubit_assignment,
                         time_MS, time_1q, output_folder,
                         out_name, error_info, n_sim, init_state,
                         initial_I, initial_trans, ion_trap,
                         synd2corr_dict, tivec, tgam, tdel,
                         save_error_dict)
                

            else:
                sim_func(chp_loc, decoder, Is_after2q, Is_after1q,
                         MS_heating, Stark, qubit_assignment,
                         time_MS, time_1q, output_folder,
                         0, error_info, n_sim, init_state,
                         initial_I, initial_trans, ion_trap,
                         synd2corr_dict, tivec, tgam, tdel,
                         save_error_dict)


        else:

            ########## Parallel part ##########

            pool = mp.Pool(n_proc)
            results = [pool.apply_async(sim_func, (chp_loc, decoder, Is_after2q,
                                                   Is_after1q, MS_heating, Stark,
                                                   qubit_assignment, times_MS,
                                                   time_1q, output_folder, i, error_info,
                                                   n_sim, init_state, initial_I,
                                                   initial_trans, ion_trap,
                                                   synd2corr_dict, tivec, tgam,
                                                   tdel, save_error_dict))
                                                for i in range(n_proc)]
            pool.close()
            pool.join()

            #dicts = [r.get() for r in results]

            #l = []
            # fix the problem that parallel runs' results will overwrite each other
            #for i in range(len(dicts)):
            #    d = dicts[i]
            #    for j in range(len(l)):
            #        for key in d.keys():
            #            if key == l[j][0]:
            #                d[key + i*n_runs_parallel] = d.pop(key)
            #                break

            #    l += d.items()
            #final_dict = dict(l)

            #for d in dicts:
            #    l += d.items()
            #final_dict = dict(l)

            ###################################  

        #if not os.path.exists(output_folder):
        #    os.makedirs(output_folder)
        #json_file = open(abs_filename, 'w')
        #json.dump(final_dict, json_file, indent=4, separators=(',', ':'),
        #          sort_keys=True)
        #json_file.close()


    return None
 

        

def run_surface17_prep(n_proc, n_times_total, chp_loc, Is_after2q,
                       error_info, decoder, output_folder):
    '''
    '''
    n_times = n_times_total/n_proc  # Make sure n_times_total is a multiple of n_proc
    sim_func = wrapper.repeat_surface17_prep
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if n_proc == 1:
        n_error = sim_func(n_times, chp_loc, Is_after2q,
                           error_info, decoder, output_folder,
                           0)
    else:

        ########## Parallel part ##########

        pool = mp.Pool(n_proc)
        results = [pool.apply_async(sim_func, (n_times, chp_loc, Is_after2q,
                                               error_info, decoder, output_folder,
                                               i))
                                            for i in range(n_proc)]
        pool.close()
        pool.join()
        
        n_error = [r.get() for r in results]

    return n_error



def parse_MC_results(operation, n_proc, output_folder, n_runs_total,
                     n_json_files=100, init_state='+Z', initial_trans=False):
    '''
    initial_trans will always be False for now.
    '''
    print 'Parsing the Monte Carlo results ...'
    n_per_json = int(float(n_runs_total)/n_json_files)
    n_serial = int(float(n_per_json)/n_proc)
    n_runs = n_serial*n_proc*n_json_files

    # n_hits_uncorr is the number of runs where there was an error
    # before applying perfet EC.  This means it counts single-qubit
    # errors that happened after the EC. 
    n_hits_uncorr = 0
    # n_hits_corr is the number of runs where there was an error only
    # after applying perfect EC.  This means that it only counts 
    # uncorrectable errors.
    n_hits_corr = 0
    n_tries = 0
    n_tries_X = 0
    n_tries_Z = 0
    
    if output_folder[-1] != '/':  output_folder += '/'

    for json_file_i in range(n_json_files):
        json_filename = str(json_file_i) + '.json'
        abs_filename = output_folder + json_filename
        json_file = open(abs_filename, 'r')
        results_dict = json.load(json_file)
        json_file.close()

        for run in results_dict:
            
            if operation[:4]=='Shor' or operation[:4] == 'Cros':
                stabs = results_dict[run][1]['stabs']
                corr_stabs = results_dict[run][1]['corr_stabs']
        
            else:
                # all the other cases: SteaneEC, KnillEC, PrepX, PrepZ
                stabs = results_dict[run]['stabs']
                corr_stabs = results_dict[run]['corr_stabs']
                num_tries = results_dict[run]['num_tries']
                
                if operation[:4]=='Prep':
                    n_tries += num_tries

                else:
                    n_tries_X += num_tries['X']
                    n_tries_Z += num_tries['Z']


            # we add 1 to n_hits_uncorr if there is at least 
            # one stabilizer with a negative eigenvalue.
            for stab in stabs:
                if stab[0] == '-':
                    n_hits_uncorr += 1
                    break

            # we add 1 to n_hits_corr only if there is one 
            # stabilizer with a negative eigenvalue after EC.
            for stab in corr_stabs:
                if stab[0] == '-':
                    n_hits_corr += 1
                    break
 
    final_dict = {}
    p_L_uncorr = float(n_hits_uncorr)/float(n_runs)
    stdev_uncorr = math.sqrt(p_L_uncorr - p_L_uncorr**2)
    p_L_corr = float(n_hits_corr)/float(n_runs)
    stdev_corr = math.sqrt(p_L_corr - p_L_corr**2)

    final_dict['p_L_uncorr'] = p_L_uncorr
    #final_dict['stdev_uncorr'] = stdev_uncorr
    final_dict['p_L_corr'] = p_L_corr
    #final_dict['stdev_corr'] = stdev_corr

    if operation[:4]=='Prep':
        final_dict['ave_n_tries'] = float(n_tries)/float(n_runs)
    elif operation[:3]=='Ste' or operation[:3]=='Kni':
        final_dict['ave_n_tries_X'] = float(n_tries_X)/float(n_runs)
        final_dict['ave_n_tries_Z'] = float(n_tries_Z)/float(n_runs)
        
    output_filename = 'MC_summary.json'
    output_filename = output_folder + output_filename
    out_file = open(output_filename, 'w')
    json.dump(final_dict, out_file, indent=4, separators=(',', ':'),
              sort_keys=True)
    out_file.close()

    return p_L_corr


def run_MC_and_parse(operation, error_info_or_file, n_proc, output_folder='./MC_results/',
                     n_runs_total=10, n_json_files=100, init_state='+Z', initial_trans=False,
                     sampling='old', error_subsets=[2,3,4], faulty_gates=['I']):
    '''
    '''
    # Interpret the error information and decide the number of runs to carry.
    error_info = read_error_info(error_info_or_file)

    try:
        n_runs_total = int(n_runs_total)
    except TypeError:
        # This occurs when the user doesn't specify the number of runs,
        # so it's decided based on the lower error rate.
        min_p = min(error_info.rates.values())
        n_runs_total = 1./(min_p**2)
        if min_p >= 1.0e-3:
            n_runs_total *= 10
        if min_p >= 1.0e-2:
            n_runs_total *= 10
        if min_p >= 0.1:
            n_runs_total = 500000

    # Run the Monte Carlo simulation
    pre_run_MC(operation, error_info, n_proc, output_folder, n_runs_total,
               n_json_files, init_state, initial_trans, sampling,
               error_subsets, faulty_gates)

    # Parse the output
    p_L_corr = parse_MC_results(operation, n_proc, output_folder, n_runs_total,
                                n_json_files, init_state, initial_trans)

    return p_L_corr

