import sys
import os
import copy
import json
import itertools as it
from circuit import *
import correction as cor
import steane as steane
import cross
import fivequbit
import surface17 as surf17
#import coder as coder
import chper_extended as chper
import error
# import visualizer as vis
from visualizer import browser_vis as brow
from scipy.stats import mode
import numpy as np
import scipy.misc
import random as rd
import time as t
#import MLD
import qcircuit_functions as qfun



def combine_stabilizers(data_stabs, data_destabs, anc_stabs, anc_destabs):
    '''
    combines data and ancillary stabilizers and destabilizers by
    adding Is at the end of the data and in front of the ancilla
    and them joining the two lists
    This function is copied exactly from chper_exteneded.py.
    The reason it's included here is to be able to call it without
    having to define a Chper object.
    '''
    n_data, n_anc = len(data_stabs), len(anc_stabs)
    new_data_stabs, new_data_destabs = [], []
    new_anc_stabs, new_anc_destabs = [] , []
    
    for i in range(n_data):
        new_data_stabs += [data_stabs[i][:] + ''.join(['I' for j in range(n_anc)])]
        new_data_destabs += [data_destabs[i][:] + ''.join(['I' for j in range(n_anc)])]
    
    extra_Is = ''.join(['I' for i in range(n_data)])
    
    for i in range(n_anc):
        new_anc_stabs += [anc_stabs[i][0] + extra_Is + anc_stabs[i][1:]]
        new_anc_destabs += [anc_destabs[i][0] + extra_Is + anc_destabs[i][1:]]

    comb_stabs = new_data_stabs[:] + new_anc_stabs[:]
    comb_destabs = new_data_destabs[:] + new_anc_destabs[:]

    return comb_stabs, comb_destabs



def create_dep_noise_error_dict(gates, error_rate, error_kind=1):
    '''
    creates a dictionary with the error information.
    Only applicable to depolarizing noise
    '''
    error_rates = {}
    for gate in gates:
        error_rates[gate] = error_rate

    one_qubit_ratio = {}
    for error in ['X', 'Y', 'Z']:
            one_qubit_ratio[error] = 1

    two_qubit_ratio = {}
    for error in itertools.product('IXYZ', 'IXYZ'):
            error_str = ''.join(error)
            two_qubit_ratio[error_str] = 1


    error_dict = {
                    'error_rates': error_rates,
                    'one_qubit_ratio': one_qubit_ratio,
                    'two_qubit_ratio': two_qubit_ratio,
                    'error_kind': error_kind
             }

    return error_dict



def split_circuit(circ, enc_circ):

    list_indexes = []
    for gate in circ.gates:
        if gate.gate_name[:2] == 'EC':
            list_indexes += [circ.gates.index(gate)]

    elem_remove = []
    for i in range(len(list_indexes[:-1])):
        if list_indexes[i+1] == list_indexes[i] + 1:
            elem_remove += [list_indexes[i]]
    for elem in elem_remove:
        list_indexes.remove(elem)

    elem_remove = []
    for index in list_indexes[:-1]:
        current_qubit = circ.gates[index].qubits[0].qubit_id
        next_qubits = [qubit.qubit_id for qubit in circ.gates[index+1].qubits]  
        if next_qubits.count(current_qubit) == 0:
            if current_qubit < max(circ.qubits()).qubit_id:
                elem_remove += [index]
    for elem in elem_remove:
        list_indexes.remove(elem)

    
    #print list_indexes

    list_subcircuits, list_subcircuits_enc = [], []
    first_index = 0
    for index in list_indexes:
        sub_circuit = Circuit(gates = circ.gates[first_index:index+1])
        sub_circuit_enc = Circuit(gates = enc_circ.gates[first_index:index+1])
        sub_circuit.update_map()
        sub_circuit_enc.update_map()
        list_subcircuits += [sub_circuit]
        list_subcircuits_enc += [sub_circuit_enc]
        first_index = index + 1
    sub_circuit = Circuit(gates = circ.gates[first_index:len(circ.gates)])
    sub_circuit_enc = Circuit(gates = enc_circ.gates[first_index:len(circ.gates)])
    sub_circuit.update_map()
    sub_circuit_enc.update_map()
    list_subcircuits += [sub_circuit]
    list_subcircuits_enc += [sub_circuit_enc]

    return [list_subcircuits, list_subcircuits_enc]


def outcome_to_syndrome_gate(EC_gate, dic, s=6, w=4):
    """
    Specific for the Steane code, but easliy generalizable.
    """
    EC = EC_gate.circuit_list[0].gates[0]
    stab_list = []
    for i in range(s):
        stab_list += [EC.circuit_list[0].gates[-w:]]
        if isinstance(EC, Encoded_Gate):
            EC = EC.circuit_list[0].gates[0]
    
    Z_outcome, X_outcome = [], []
    for stab in stab_list[:s/2]:
        Z_outcome += [sum([dic[gate][0] for gate in stab])%2]
    for stab in stab_list[s/2:]:
        X_outcome += [sum([dic[gate][0] for gate in stab])%2]
    
    X_error = steane.Code.stabilizer_syndrome_dict[tuple(Z_outcome)]
    Z_error = steane.Code.stabilizer_syndrome_dict[tuple(X_outcome)]

    try:
        X = X_error.index('E')
    except ValueError:
        X = 'n'
    try:
        Z = Z_error.index('E')
    except ValueError:
        Z = 'n'

    return [X, Z]



def parity_for_one_measurement(first_ancilla, dic, w=4):
    return sum([dic[first_ancilla + i][0] for i in range(w)])%2


def majority_vote_for_one_stabilizer(first_ancilla, dic, redun=3, w=4):
    majority_list = []
    for i in range(redun):
        majority_list += [parity_for_one_measurement(first_ancilla + i*w, dic, w)]
    
    return int(mode(majority_list)[0][0])


def outcome_to_syndrome(first_ancilla, dic, redun=3, s=6, w=4,
            diVincenzo=False):
    """
    Specific for the Steane code, but easliy generalizable.
    We're assuming that the Z stabilizers come first.
    """
    Z_outcome, X_outcome = [], []
    for i in range(s/2):
        Z_outcome += [majority_vote_for_one_stabilizer(first_ancilla + i*redun*w, dic, redun, w)]
    for i in range(s/2,s):
        X_outcome += [majority_vote_for_one_stabilizer(first_ancilla + i*redun*w, dic, redun, w)]
    
    X_error = steane.Code.stabilizer_syndrome_dict[tuple(Z_outcome)]
    Z_error = steane.Code.stabilizer_syndrome_dict[tuple(X_outcome)]

    #print Z_outcome
    #print X_outcome
    #print X_error 
    #print Z_error

    try:
        X = X_error.index('E')
    except ValueError:
        X = 'n'
    try:
        Z = Z_error.index('E')
    except ValueError:
        Z = 'n'

    return [X, Z]




def run_one_subcircuit(circ, enc_circ, num_d_qub, chp_location, stabs=[], destabs=[], add_error = False, error_file_or_info = 'error.json', redun=3, s=6, w=4):
    
    if add_error:
        error.add_error(enc_circ, error_file_or_info)

    enc_circ_chp = chper.Chper(circ=enc_circ,num_d_qub=num_d_qub,stabs=stabs,destabs=destabs)
    dic, final_stabs, final_destabs = enc_circ_chp.run(chp_location)
    
    EC = []
    for gate in circ.gates:
        if gate.gate_name[:2] == 'EC': EC += [gate.qubits[0].qubit_id]

    if len(EC) == 0:    # If there's no EC, we assume it's the last sub-circuit.
        return dic, final_stabs, final_destabs

    i = 0
    corrections = []
    
    for i in range(len(EC)):
        syndrome = outcome_to_syndrome(num_d_qub + i*redun*s*w, dic, redun, s, w)
        if syndrome[0] != 'n':
            correction = 'Z' + str( (s+1)*EC[i] + syndrome[0] )
            corrections += [correction]
        if syndrome[1] != 'n':
            correction = 'X' + str( (s+1)*EC[i] + syndrome[1] )
            corrections += [correction]
        

    #for gate in enc_circ.gates:
    #   if gate.gate_name[:2] == 'EC':
    #       correction = outcome_to_syndrome(gate, dic)
    #       log_qubit_num = circ.gates[i].qubits[0].qubit_id
    #       if correction[0] != 'n':
    #           corrections += [7*log_qubit_num + correction[0]]
    #       if correction[1] != 'n':
    #           corrections += [7*log_qubit_num + correction[1]]    
    #   i += 1

    return corrections, final_stabs, final_destabs


def run_whole_circuit(circ, enc_circ, chp_location, error_file_or_info):
    num_d_qub = len(enc_circ.data_qubits())
    sub_circuits = split_circuit(circ, enc_circ)
    num_sub_circuits = len(sub_circuits[0])
    stabs, destabs = [], []
    for i in range(num_sub_circuits):
        curr_circ = sub_circuits[0][i]
        curr_enc_circ = sub_circuits[1][i].unpack()

        #vis.Printer.print_circuit(sub_circuits[1][i], False, True, 'packed.txt')
        #vis.Printer.print_circuit(curr_enc_circ, False, True, 'unpacked.txt')


        #print "Subcircuit %i\n" %i
    
        if i == 0:  # temporary
            dic, stabs, destabs = run_one_subcircuit(curr_circ, curr_enc_circ, 
                                                     num_d_qub, chp_location, 
                                                     stabs, destabs)             

        else:
            dic, stabs, destabs = run_one_subcircuit(curr_circ, curr_enc_circ,
                                                     num_d_qub, chp_location, 
                                                     stabs, destabs, True, 
                                                     error_file_or_info)
        
        if i < num_sub_circuits-1:
            for corr in dic:
                sub_circuits[1][i+1].insert_gate('', [int(corr[1])], 'data', corr[0], True)
        else:
            #vis.Printer.print_circuit(curr_enc_circ)
            #print dic
            return measurement_interpretation_steane(dic)



def run_ShorEC(chp_location, initi_stabs, initi_destabs, Is_after2q,
               add_error=False, error_info=None, CHP_IO_files=False,
               state='+Z', carry_run=False, sampling='old', n_errors=2,
               gate_indices=[(0,0)], initial_I=True, initial_trans=False):
    '''
    sampling:  'old' or 'Muyalon', in honor of Muyuan and Alonzo.
    n_errors:  if sampling=='Muyalon', it refers to the number of errors
               to insert.
    gate_indices:  indices for the gates after which we can add errors.
    '''

    # 1. Generate the subcircuits.
    subcircs = create_EC_subcircs('Shor',
                                  Is_after2q,
                                  initial_I,
                                  initial_trans)


    run_dict = {}

    #brow.from_circuit(subcircs[1], True)

    # this is used in the error analysis to distinguish
    # between correctable and non-correctable errors.
    # we run a perfect EC step and get the final results.
    if carry_run:
        stabs, destabs, n_X, n_Z = run_CatEC_3(subcircs,
                                               7, False,
                                               chp_location,
                                               None,
                                               initi_stabs,
                                               initi_destabs,
                                               state,
                                               CHP_IO_files)
        return stabs, destabs

    

    carry_run = False
    
    if sampling == 'old':
        
        for i in range(6):
            subcirc = Circuit(gates=[EC_circ.gates[i]])
            subcirc.update_map()
            subcirc = subcirc.unpack()
            ### this line is a quick fix only for Alonzo's scripts ###
            error.add_Is(subcirc)
            ##########################################################
            if add_error:
                error.add_error(subcirc, error_info)

            subcircs += [subcirc]

            subcirc_dict, local_carry_run = get_errors_dict_Shor(subcirc)
            # As long as at least one of the local_carry_run is True,
            # we'll carry the run.
            if i >= 4:  local_carry_run = False
            carry_run = carry_run or local_carry_run
            errors_dict[i] = subcirc_dict


    elif sampling == "Muyalon":
        errors_dict, carry_run = add_errors_fast_sampler(gate_indices,
                                                         n_errors,
                                                         subcircs,
                                                         error_info)

    if carry_run:

        stabs, destabs, n_X, n_Z = run_CatEC_3(subcircs,
                                               7, False,
                                               chp_location,
                                               None,
                                               initi_stabs,
                                               initi_destabs,
                                               state,
                                               CHP_IO_files)
        
        #if stabs == initi_stabs:
            #run_dict[state] = ('ne', n_X, n_Z) 
        #else:
            #run_dict[state] = {'stabs': stabs, 'destabs': destabs} 
            #run_dict[state] = (stabs, destabs, n_X, n_Z)
            #condition = True

        run_dict = {'stabs':  stabs,
                    'destabs':  destabs,
                    'n_X':  n_X,
                    'n_Z':  n_Z
                   }

        return (errors_dict, run_dict)
       
    
    else:
        
        return None    



def run_5qubit(chp_location, initi_stabs, initi_destabs, Is_after2q, add_error=False,
               error_info=None, CHP_IO_files=False, state='+Z', carry_run=False,
               sampling='old', n_errors=2, gate_indices=[(0,0)], initial_I=True,
               initial_trans=False):
    '''
    '''


    # 1. Generate the subcircuits.
    subcircs = create_EC_subcircs('ShorEC_fivequbit',
                                  Is_after2q,
                                  initial_I,
                                  initial_trans)

    run_dict = {}

    #brow.from_circuit(subcircs[1], True)

    # 2. If carry_run is true, then call run_fivequbit_3 after doing the 
    #    initial subcircuit decomposition
    if carry_run:
        stabs, destabs, n_times = run_5qubit_3(subcircs,
                                               5, False,
                                               chp_location,
                                               None,
                                               initi_stabs,
                                               initi_destabs,
                                               state,
                                               CHP_IO_files)

        return stabs, destabs

    
    
    if sampling == 'old':
        pass    


    elif sampling == "Muyalon":
        errors_dict, carry_run = add_errors_fast_sampler(gate_indices,
                                                         n_errors,
                                                         subcircs,
                                                         error_info)


    # 4. If carry_run is set to True after the errors have been added,
    #    we run the circuit.
    if carry_run:
        stabs, destabs, n_times = run_5qubit_3(subcircs,
                                               5, False,
                                               chp_location,
                                               None,
                                               initi_stabs,
                                               initi_destabs,
                                               state,
                                               CHP_IO_files)
    
        # 5. Set up the run dictionary
        run_dict = {'stabs':  stabs,
                    'destabs': destabs,
                    'n':  n_times
                   }

        return (errors_dict, run_dict)


    else:

        return None


        
def run_CrossEC(chp_location, initi_stabs, initi_destabs, Is_after2q, add_error=False,
                error_info=None, CHP_IO_files=False, state='+Z', carry_run=False,
                sampling='old', n_errors=2, gate_indices=[(0,0)], initial_I=True, 
                initial_trans=False): 
    '''
    '''

    # 1. Generate the subcircuits
    subcircs = create_EC_subcircs('Cross',
                                  Is_after2q,
                                  initial_I,
                                  initial_trans)

    run_dict = {}

    # this is used in the error analysis to distinguish
    # between correctable and non-correctable errors.
    # we run a perfect EC step and get the final results.
    if carry_run:
        stabs, destabs, n_times = run_Cross_3(subcircs,
                                               7, False,
                                               chp_location,
                                               None,
                                               initi_stabs,
                                               initi_destabs,
                                               state,
                                               CHP_IO_files)
        return stabs, destabs


    carry_run = False

    if sampling == 'old':
        pass
    
    elif sampling == "Muyalon":
        errors_dict, carry_run = add_errors_fast_sampler(gate_indices,
                                                         n_errors,
                                                         subcircs,
                                                         error_info)
        

    if carry_run:

        stabs, destabs, n_times = run_Cross_3(subcircs,
                                               7, False,
                                               chp_location,
                                               None,
                                               initi_stabs,
                                               initi_destabs,
                                               state,
                                               CHP_IO_files)
 
        run_dict = {'stabs':  stabs,
                    'destabs':  destabs,
                    'n':  n_times
                   }

        return (errors_dict, run_dict)
       

    else:
        
        return None 


 

def run_simulation_CatEC_3(n_runs, EC_circuit, initial_states, add_error, error_info,
                           chp_location, json_filename, CHP_IO_files=True,
                           first_run_index=0):
    '''
    Think of a better name later on
    CHP_IO_files:    True if we want CHP to deal with (txt) files for the I/O operations.
                     False if we want strings
    first_run_index: the index of the first run.  This was added later on to make it easier
                     to merge the final dictionaries when running this function in a parallel
                     way.
    Notice that json_filename is not used anymore as the final dictionary is returned
    by the function rather than dumped to a json file.
    '''
    output_dict = {}
    run_range = range(n_runs)

    for run in run_range:

        #if run % notice == 0:  print 'I\'m on run %i' %run
        #condition = False
        run_dict = {}
        subcircs = []
        errors_dict = {}
        carry_run = False

        for i in range(6):
            subcirc = Circuit(gates=[EC_circuit.gates[i]])
            subcirc.update_map()
            subcirc = subcirc.unpack()
            if add_error:
                error.add_error(subcirc, error_info)

            subcircs += [subcirc]

            data_errors, anc_errors = [], []
            for gate in subcirc.gates:
                if gate.is_error:
                    if i < 4:
                        carry_run = True
                    qubit = gate.qubits[0]
                    if qubit.qubit_type == 'data':
                        data_errors += [(qubit.qubit_id,
                                         gate.gate_name)]
                    elif qubit.qubit_type == 'ancilla':
                        anc_errors += [(qubit.qubit_id,
                                        gate.gate_name)]

            if len(data_errors) > 0 or len(anc_errors) > 0:
                errors_dict[i] = {'d': data_errors, 'a': anc_errors}

        
        if carry_run:

            for state in initial_states:

                #print 'state =', state
                #print initial_states[state]

                stabs, destabs, n_X, n_Z = run_CatEC_3(subcircs,
                                                        7, False,
                                                        chp_location,
                                                        None,
                                                        initial_states[state][0],
                                                        initial_states[state][1],
                                                        state,
                                                        CHP_IO_files)
        
                if stabs == initial_states[state][0]:
                    run_dict[state] = ('ne', n_X, n_Z)
                else:
                    #run_dict[state] = {'stabs': stabs, 'destabs': destabs} 
                    run_dict[state] = (stabs, n_X, n_Z)
                    #condition = True

            output_dict[run + first_run_index] = (errors_dict, run_dict)

    #json_file = open(json_filename, 'w')
    #json.dump(output_dict, json_file, indent=4, separators=(',', ':'),
    #          sort_keys=True)
    #json_file.close()

    return output_dict



def run_CatEC_3(circ, num_d_qub, add_error, chp_location, error_file_or_info,
                initial_stabs=[], initial_destabs=[], init_state='+Z', 
                CHP_IO_files=True, always_run_3=False):
    '''
    Specific to Steane code with Shor ancilla and 
    redundancy 3
    '''

    # number of times we measure each set of stabilizers
    n_X, n_Z = 0, 0
    
    # If circ is a list, we assume that the sub-circuits are 
    # already unpacked and with the errors inserted.
    if type(circ) == type([]):
        Xs = [circ[2*i] for i in range(3)]
        Zs = [circ[2*i+1] for i in range(3)]
    
    else:
        EC_circ = circ.gates[0].circuit_list[0]
        Xs = [Circuit(gates=[EC_circ.gates[2*i]]) for i in range(3)]
        Zs = [Circuit(gates=[EC_circ.gates[2*i+1]]) for i in range(3)]
        for i in range(3):
            Xs[i].update_map()
            Xs[i] = Xs[i].unpack()
            Zs[i].update_map()
            Zs[i] = Zs[i].unpack()
            if add_error:
                error.add_error(Xs[i], error_file_or_info)
                error.add_error(Zs[i], error_file_or_info)
        

    # create initial stabs and destabs if non
    
    if initial_stabs == [] and initial_destabs == []:
        stabs, destabs = prepare_stabs_Steane(init_state,
                            chp_location,
                            CHP_IO_files)
    
    else:
        stabs = copy.deepcopy(initial_stabs)
        destabs = copy.deepcopy(initial_destabs)    

    # make this generazible to more initial states by
    # copying methods from new_musiqc/steane.py to
    # new_musiqc/scriptsIARPA/steane.py

    #browser_vis.from_circuit(Xs[0],True)   
    #vis.Printer.print_circuit(X1)
    #qubits = Zs[1].qubits()    
    #print qubits[0]

    #g1 = Zs[1].gates[0]
    #g2 = Zs[1].gates[2]
    #g3 = Zs[1].gates[21]
    #g4 = Zs[1].gates[38]
    #g5 = Zs[1].gates[45]
    #new_g1 = Zs[1].insert_gate(g1, [qubits[4]], '', 'Z', False)
    #new_g2 = Zs[1].insert_gate(g2, [qubits[9]], '', 'X', False)
    #new_g3 = Zs[1].insert_gate(g3, [qubits[13]], '', 'X', False)
    #new_g4 = Zs[1].insert_gate(g4, [qubits[17]], '', 'X', False)
    #new_g5 = Zs[1].insert_gate(g5, [qubits[6]], '', 'X', False)
    #new_g1.is_error = True
    #new_g2.is_error = True 
    #new_g3.is_error = True
    #new_g4.is_error = True 
    #new_g5.is_error = True 

    #browser_vis.from_circuit(Zs[1], True)
    
    #######################################################

    Z_corrs, X_corrs = [], []
    for i in range(2):
        
        n_X += 1        
        n_Z += 1    
        stabs, destabs, Z_corr = run_stabs_Steane(Xs[i], 
                               num_d_qub,
                               stabs, destabs,
                               chp_location,
                               'X',
                               CHP_IO_files)
        Z_corrs += [Z_corr]

        stabs, destabs, X_corr = run_stabs_Steane(Zs[i], 
                            num_d_qub,
                            stabs, destabs,
                            chp_location,
                            'Z',
                            CHP_IO_files)
        X_corrs += [X_corr]             

    # Decide if we need to repeat the X stabilizers
    if (Z_corrs[0] != Z_corrs[1]) or always_run_3:

        n_X += 1
        stabs, destabs, Z_corr = run_stabs_Steane(Xs[2], 
                               num_d_qub,
                               stabs, destabs,
                               chp_location,
                               'X',
                               CHP_IO_files)

    # Decide if we need to repeat the Z stabilizers
    if (X_corrs[0] != X_corrs[1]) or always_run_3:
        
        n_Z += 1
        stabs, destabs, X_corr = run_stabs_Steane(Zs[2], 
                               num_d_qub,
                               stabs, destabs,
                               chp_location,
                               'Z',
                               CHP_IO_files)
    if 'Z' in Z_corr:
        stabs, destabs = update_stabs(stabs, destabs, Z_corr)
    if 'X' in X_corr:
        stabs, destabs = update_stabs(stabs, destabs, X_corr)

    return stabs, destabs, n_X, n_Z



def run_Cross_3(circ, num_d_qub, add_error, chp_location, error_file_or_info,
                initial_stabs=[], initial_destabs=[], init_state='+Z', 
                CHP_IO_files=True):
    '''
    Specific to Cross code with redundancy 3
    '''

    # number of times we measure the stabilizers
    n_times = 0
    
    # If circ is a list, we assume that the sub-circuits are 
    # already unpacked and with the errors inserted.
    if type(circ) == type([]):
        circ_list = circ[:]
        #print 'We have a list'
        
    else:
        EC_circ = circ.gates[0].circuit_list[0]
        circ_list = [Circuit(gates=[EC_circ.gates[i]]) for i in range(3)]
        for i in range(3):
            circ_list[i].update_map()
            circ_list[i] = circ_list[i].unpack()
    
    # Should be False if we are running the MC simulations, because
    # the errors should have been added at a previous step.

    if add_error:
        error.add_error(circ_list[i], error_file_or_info)

    # create initial stabs and destabs if non 
    if initial_stabs == [] and initial_destabs == []:
        stabs, destabs = prepare_stabs_Cross(init_state,
                            chp_location,
                            CHP_IO_files)
    
    else:
        stabs = copy.deepcopy(initial_stabs)
        destabs = copy.deepcopy(initial_destabs)    

    
    # run and correct
    corrs = []
    for i in range(2):
        
        n_times += 1    
        stabs, destabs, corr = run_stabs_Cross(circ_list[i], 
                               num_d_qub,
                               stabs, destabs,
                               chp_location,
                               CHP_IO_files)
        corrs += [corr]


    # Decide if we need to repeat the stabilizers
    if corrs[0] != corrs[1]:

        n_times += 1
        stabs, destabs, corr = run_stabs_Cross(circ_list[2], 
                               num_d_qub,
                               stabs, destabs,
                               chp_location,
                               CHP_IO_files)

    stabs, destabs = update_stabs(stabs, destabs, corr)

    return stabs, destabs, n_times



def run_5qubit_3(circ, num_d_qub, add_error, chp_location,
                 error_file_or_info, initial_stabs=[],
                 initial_destabs=[], init_state='+Z',
                 CHP_IO_files=True):
    '''
    '''

    #1 If statement for circuit decomposition 
    n_times = 0

    if type(circ) == type([]):
        circ_list = circ[:]
    else:
        EC_circ = circ.gates[0].circuit_list[0]
        circ_list = [Circuit(gates=[EC_circ.gates[i]]) for i in range(3)]
        for i in range(3):
            circ_list[i].update_map()
            circ_list[i] = circ_list[i].unpack()

    #2 IF statement for adding the error.  No errors should be added at this
    # stage as they were added before. 
    if add_error:
        error.add_error(circ_list[i], error_file_or_info)

    if initial_stabs == [] and initial_destabs == []:
        stabs, destabs = prepare_stabs_5qubit(init_state,
                                              chp_location,
                                              CHP_IO_files)
    else:
        stabs = copy.deepcopy(initial_stabs)
        destabs = copy.deepcopy(initial_destabs)

    
    corrs = []
    for i in range(2):
        n_times += 1
        stabs, destabs, corr = run_stabs_5qubit(circ_list[i],
                                                num_d_qub,
                                                stabs, destabs,
                                                chp_location,
                                                CHP_IO_files)
        corrs += [corr]

    #5 If statement to decide whether we need to repeat the stabilizers
    if corrs[0] != corrs[1]:
        n_times += 1
        stabs, destabs, corr = run_stabs_5qubit(circ_list[2],
                                                num_d_qub,
                                                stabs, destabs,
                                                chp_location,
                                                CHP_IO_files)

    #6 After updating the stabs, return the stabilizers,destabilizers, and n_times parameters
    stabs, destabs = update_stabs(stabs, destabs, corr)

    return stabs, destabs, n_times

    



def run_stabs_Steane(circ, n_d_q, stabs, destabs, chp_location,
             stab_kind='X', CHP_IO_files=True):
    '''
    '''
    n_a_q = len(circ.ancilla_qubits())
    circ_chp = chper.Chper(circ=circ, num_d_qub=n_d_q, num_a_qub=n_a_q,
                           stabs=stabs, destabs=destabs,
                           anc_stabs=[], anc_destabs=[],
                           input_output_files=CHP_IO_files)
    dic, final_stabs, final_destabs = circ_chp.run(chp_location)
    data_corr, prop_corr = stabs_Steane(dic, n_d_q, 3, 4, stab_kind)

    #print 'data correction =', data_corr
    #print 'prop correction =', prop_corr   

    #print 'Stabs before X corr =', final_stabs
    #print 'Destabs before X corr =', final_destabs
    
    if stab_kind in prop_corr:
        corr_stabs, corr_destabs = update_stabs(final_stabs,
                                  final_destabs,
                                  prop_corr)
    else:
        corr_stabs, corr_destabs = final_stabs, final_destabs
    
    #print 'Stabs after X corr =', corr_stabs
    #print 'Destabs after X corr =', corr_destabs

    return corr_stabs, corr_destabs, data_corr  



def run_stabs_5qubit(circ, n_d_q, stabs, destabs, chp_location,
                     CHP_IO_files=True):
    '''
    '''
    #1 Get the number of data and ancilla qubits 
    n_a_q=len(circ.ancilla_qubits())
    n_d_q=len(circ.data_qubits())
    #2 Create chper object called circ_chp.run(chp_location) to run the simulation 
    #print 'stabs=',stabs
    #print 'destabs=', destabs
    circ_chp=chper.Chper(circ=circ, num_d_qub=n_d_q, num_a_qub=n_a_q,
                         stabs=stabs, destabs=destabs, anc_stabs=[],
                         anc_destabs=[], input_output_files=CHP_IO_files)
    dic, final_stabs, final_destabs=circ_chp.run(chp_location)
    data_corr, prop_corr=stabs_5qubit(dic, n_d_q, 4, 4)
    corr_stabs, corr_destabs=update_stabs(final_stabs,
                                          final_destabs,
                                          prop_corr)

    return final_stabs, final_destabs, data_corr



def run_stabs_Cross(circ, n_d_q, stabs, destabs, chp_location, CHP_IO_files=True):
    '''
    '''
    n_a_q = len(circ.ancilla_qubits())
    
    #n_d_q = 21
    #stabs, destabs = [], []
    #for i in range(21):
    #    stab = ['Z' if i==j else 'I' for j in range(21)]
    #    stab.insert(0, '+')
    #    destab = ['X' if i==j else 'I' for j in range(21)]
    #    destab.insert(0, '+')
    #    stabs += [''.join(stab)]
    #    destabs += [''.join(destab)]

    #print n_d_q, stabs, destabs
    circ_chp = chper.Chper(circ=circ, num_d_qub=n_d_q, num_a_qub=n_a_q,
                           stabs=stabs, destabs=destabs,
                           anc_stabs=[], anc_destabs=[],
                           input_output_files=CHP_IO_files)
    dic, final_stabs, final_destabs = circ_chp.run(chp_location)
    data_corr = stabs_Cross(dic, n_d_q)

    return final_stabs, final_destabs, data_corr 
    
    
    
def run_stabs_surface17(circ, n_d_q, n_a_q, stabs, destabs, chp_location,
                        CHP_IO_files=True, stabilizers='all'):
    '''
    '''
    circ_chp = chper.Chper(circ=circ, num_d_qub=n_d_q, num_a_qub=n_a_q,
                           stabs=stabs[:], destabs=destabs[:],
                           anc_stabs=[], anc_destabs=[],
                           input_output_files=CHP_IO_files)
    dic, final_stabs, final_destabs = circ_chp.run(chp_location)

    #brow.from_circuit(circ, True)
    #print dic
    #sys.exit(0)

    if stabilizers == 'all':
        synd = ''.join(map(str,[dic[10][0], dic[12][0], dic[13][0], dic[15][0], 
                                dic[9][0], dic[11][0], dic[14][0], dic[16][0]]))     
    else:
        synd = ''.join(map(str,[dic[9][0], dic[10][0], dic[11][0], dic[12][0]]))

    return final_stabs, final_destabs, synd



def run_surface17EC(chp_location, initi_stabs, initi_destabs, Is_after2q, 
                    Is_after1q, MS_heating, Stark, qubit_assignment,
                    time_MS, time_1q, add_error, error_info, n_errors, 
                    gate_indices, sampling, CHP_IO_files=False, 
                    carry_run=False, initial_I=True, initial_trans=False, 
                    previous_errors=False, add_errors_debugging=False, 
                    stabilizers='all', ion_trap=False): 
    '''
    add_errors_debugging should only be True for debugging purposes.
    if ion_trap == True:  the circuit uses MS gates
    else:                 the circuit uses CXs and CZs.
    Is_after_2q:  whether or not to add I gates after 2-qubit gates.
    Is_after_1q:  whether or not to add I gates after 1-qubit gaets.
    This makes the addition of depolarizing Pauli noise easier.
    MS_heating:  whether or not we will include heating after MS gaets.
    Stark:  whether or not we will include Stark shifts.
    qubit_assignment:  a list of integers indicating the qubit to ion
                       assignment.  This is needed to calculate the 
                       duration of each MS gate.
    time_MS:  a function of the duration of the MS gate in terms of the
              distance between the two ions.  Muyuan's fit of Luming's 
              data.
    time_1q:  the duration of 1-qubit rotations in microseconds (float).
    '''

    # 1. Generate the circuit
    meas_errors = True
    if stabilizers == 'all':
        stabs = surf17.Code.stabilizers[:]
        ancillae_indexes = surf17.Code.ancillae_indexes[:]
    elif stabilizers == 'X':
        stabs = surf17.Code.stabilizers[:4]
        ancillae_indexes = surf17.Code.ancillae_indexes[:4]
    elif stabilizers == 'Z':
        stabs = surf17.Code.stabilizers[4:]
        ancillae_indexes = surf17.Code.ancillae_indexes[4:]

    ECcirc = cor.Surface_Code_ion_trap.generate_stabs_meas(
                                stabs,
                                ancillae_indexes,
                                initial_I,
                                meas_errors,
                                ion_trap,
                                add_errors_debugging,
                                Is_after2q, 
                                Is_after1q,
                                MS_heating,
                                Stark,
                                qubit_assignment,
                                time_MS,
                                time_1q)
    ECcirc = ECcirc.unpack()
    if stabilizers == 'all':
        ECcirc.to_ancilla(range(9,17))
    elif stabilizers == 'X':
        ECcirc.to_ancilla([9,11,14,16])
    elif stabilizers == 'Z':
        ECcirc.to_ancilla([10,12,13,15])
    
    #print 'n_gates =', len(ECcirc.gates)
    #brow.from_circuit(ECcirc, True)
    #sys.exit(0)

    #for gate in ECcirc.gates:
    #    print gate.gate_name
    #sys.exit(0)

    run_dict = {}

    # this is used in the error analysis to distinguish
    # between correctable and non-correctable errors.
    # we run a perfect EC step and get the final results.
    if carry_run:
        if stabilizers == 'all':  n_anc = 8
        else:                     n_anc = 4
        stabs, destabs, synd = run_stabs_surface17(ECcirc, 9, n_anc,
                                                   initi_stabs,
                                                   initi_destabs,
                                                   chp_location,
                                                   CHP_IO_files,
                                                   stabilizers)
        
        return stabs, destabs, synd


    carry_run = False

    if sampling == 'old':
        if add_error:
            error.add_error(ECcirc, error_info)
            errors_dict, carry_run = get_errors_dict(ECcirc)

 
    elif sampling == "Muyalon":
        if add_error:
            errors_dict, carry_run = add_errors_fast_sampler_surface(
                                                             gate_indices,
                                                             n_errors,
                                                             ECcirc,
                                                             error_info)
    
    if previous_errors: carry_run = True
    if carry_run:

        if stabilizers == 'all':  n_anc = 8
        else:                     n_anc = 4
        stabs, destabs, synd = run_stabs_surface17(ECcirc, 9, n_anc,
                                                   initi_stabs,
                                                   initi_destabs,
                                                   chp_location,
                                                   CHP_IO_files,
                                                   stabilizers)
        run_dict = {'stabs':  stabs,
                    'destabs':  destabs,
                    'syndrome': synd
                   }

        return (errors_dict, run_dict)
       

    else:
        
        return None 



def update_stabs(stabs, destabs, operation):
    '''
    This function is used exclusively when operation only has
    X or Z.  This means that the only thing we change is the 
    sign of the stabilizers and destabilizers.  Notice that
    changing the sign of the destabilizers is unnecessary,
    because either way D_i would anticommute with S_i and 
    commute with D_j and S_j:
    If {A,B} = 0, then {A,-B} = 0
    If [A,B] = 0, then [A,-B] = 0
    For this reason, and to save us some time, I have
    commented it out.
    '''
    if len(stabs) != len(destabs):
        raise IndexError('Check out the stabs and destabs.')

    n = len(stabs)
    for i in range(n):
        stabs[i] = update_one_stab(stabs[i], operation)
        #destabs[i] = update_one_stab(destabs[i], operation)

    return stabs, destabs
    


def update_one_stab(state, operator):
    '''
    Assumes state has sign, but operator doesn't, i.e.,
    it's always positive.
    For example: 
       - state:      '+XIXIXIX'  (can have X, Y, or Z)
       - operator:    'ZIIIIII'  (can only have X or Z)
    '''
    state = list(state)
    operator = list(operator)
    par = 0
    for i in range(len(operator)):
        if str(operator[i]) == 'X':
            if str(state[i+1]) == 'Y' or str(state[i+1]) == 'Z':
                par += 1
        elif str(operator[i]) == 'Z':
            if str(state[i+1]) == 'X' or str(state[i+1]) == 'Y':
                par += 1
        elif str(operator[i]) == 'Y':
            if str(state[i+1]) == 'X' or str(state[i+1]) == 'Z':
                par += 1

    if par%2 == 1:
        if   state[0] == '+':   state[0] = '-'
        elif state[0] == '-':   state[0] = '+'

    return ''.join(state)



def stabs_QEC_diVin(dic, n_first_anc, code, stab_kind=None):
    '''
    After one round of either X or Z stabs (for a CSS code)
    or the whole set of stabs (for a non-CSS code),
    with an unverified cat state, we first correct the hook 
    errors (errors that propagated from ancilla to data),
    and then return the correction for the data errors.
    '''
    
    # so far only Steane and 5qubit codes
    
    extra_s = 0  # little hack to get the stabilizer right

    if code == 'Steane':
        n_q, s, w = 7, 3, 4
        code_class = steane.Code
        if type(stab_kind) != type(''):
            raise TypeError('stab_kind either X or Z.')
        if stab_kind == 'Z':
            extra_s = 3
            error = 'X'
        elif stab_kind == 'X':
            error = 'Z'


    elif code == '5qubit':
        n_q, s, w = 5, 4, 4
        code_class = fivequbit.Code


    # 1 Raise exception if dictionary is not right length 
    if len(dic) != s*w:
        raise IndexError('Dictionary does not have right length.')

    # 2 Assign data_corr and prop_corr
    prop_corr = ['I' for i in range(n_q)]
    data_error = []

    for i in range(s):
        outcomes = [dic[n_first_anc + i*w + j][0]
                                for j in range(w)]
        data_error += [outcomes.pop(1)]
        if 0 not in outcomes:
            stab = code_class.stabilizer[i + extra_s]
            prop_corr = update_errors_anc(prop_corr,
                                          stab)
    
    data_corr = code_class.stabilizer_syndrome_dict[tuple(data_error)]
    if code == 'Steane':
        data_corr = [i if i == 'I' else error for i in data_corr]   
    
    return data_corr, prop_corr



def stabs_Steane(dic, n_first_anc, s=3, w=4, stab_kind='X'):
    '''
    one round of either X or Z
    '''
    if len(dic) != s*w: 
        raise IndexError('The dict should be of length 12.')
    
    if stab_kind == 'X':    error = 'Z'
    elif stab_kind == 'Z':  error = 'X'

    prop_corr = ['I' for i in range(7)]
    data_error = []
    
    for i in range(s):
        outcomes = [dic[n_first_anc + i*w + j][0]
                                for j in range(w)]
        
        data_error += [outcomes.pop(1)]
        if 0 not in outcomes:
            stab = steane.Code.stabilizer[i]
            prop_corr = update_errors_anc(prop_corr,
                                          stab,
                                          stab_kind)

    data_corr = steane.Code.stabilizer_syndrome_dict[tuple(data_error)]
    data_corr = [i if i == 'I' else error for i in data_corr]   

    return data_corr, prop_corr


def stabs_5qubit(dic, n_first_anc, s=4, w=4):
    '''
    '''

    #1 Raise exception if dictionary is not right length 
    if len(dic) != 16:
        raise IndexError('Dictionary should be of length 16.')

    #2 Assign data_corr and prop_corr
    prop_corr = ['I' for i in range(5)]
    data_error = []

    for i in range(s):
        outcomes = [dic[n_first_anc + i*w + j][0]
                            for j in range(w)]
        data_error += [outcomes.pop(1)]
        if 0 not in outcomes:
            stab = fivequbit.Code.stabilizer[i]
            prop_corr = update_errors_anc(prop_corr, stab)

    data_corr = fivequbit.Code.stabilizer_syndrome_dict[tuple(data_error)]
    
    return data_corr, prop_corr



def stabs_Cross(dic, n_first_anc):
    '''
    '''
    if len(dic) != 6:
        raise IndexError('The dict should be of length 6.')

    outcome = [dic[n_first_anc + i][0] for i in range(6)]
    

    ##########################################################
    # new makeshift part added to start generating results
    # We need to change it.  MGA 12/18/2015.
    
    # First, turn outcome from binary to decimal
    outcome_integer = 0
    for i in range(len(outcome)):
        outcome_integer += outcome[i]*2**(5-i)

    # Then look up the corresponding syndrome.
    # for now, if the error syndrome is not in the 
    # lookuptable, we don't correct it.
    # This means that there are some 2-qubit errors
    # and possibly higher weight errors that we 
    # are not correcting, but at least we are 
    # correcting all the errors needed to achieve FT.
    #try:
    #    corr = cross.Code.decoding_lookuptable[outcome_integer]
    #except KeyError:
    #    corr = ['I','I','I','I','I','I','I']
            
    #return corr

    ##########################################################
    
    ##########################################################
    # The complete lookup table was finally generated. 
    # It's a json file named 'complete_lookup_table.json',
    # generated with the script 'Cross_decoding.py'.
    # The json file is imported by cross.py.

    corr = cross.Code.complete_lookup_table[outcome_integer]
  
    return corr

    ##########################################################
    

    ##########################################################
    # Original naive decoder.

    #synd_dic = {0:[('Z',0)],
    #            1:[('Z',1)],
    #            2:[('Z',2)],
    #            3:[('Z',3),('X',4)],
    #            4:[('X',6)],
    #            5:[('X',4)]}

    #corr = ['I','I','I','I','I','I','I']
    #for i in range(len(outcome)):
    #    if outcome[i] == 1:
    #        for j in synd_dic[i]:
    #            if corr[j[1]] == 'I':
    #                corr[j[1]] = j[0]
    #            else:
    #                corr[j[1]] = corr[j[1]] + j[0]
    #if corr[4] == 'XX':
    #    corr[4] = 'I'
    
    #return corr 




    
#def update_errors_anc(current_errors, stab, error='X'):
#    '''
#    not error on ancilla, but errors that propagated
#    from ancilla to data
#    '''
#    l = len(stab)   # should be 7 for Steane
#    qubits_indexes = [i for i in range(l) if stab[i]!='I']
#    qubits_to_correct = qubits_indexes[:2]
#    for q in qubits_to_correct:
#        if current_errors[q] == 'I':
#            current_errors[q] = error
#        elif current_errors[q] == error:
#            current_errors[q] = 'I'
#    
#    return current_errors 



def logical_error_surface17(stabs, state='+Z'):
    '''
    Determines if there has been a logical error
    on the surface17 code.
    Currently only implemented to detect X errors
    '''
    Z_stabs = stabs[4:]
    s0 = (Z_stabs[0][1:]=='ZIIIZZIII')
    s1 = (Z_stabs[1][1:]=='IZZIZZIII')
    s2 = (Z_stabs[2][1:]=='IIIZZZIII')
    s3 = (Z_stabs[3][1:]=='IIIIIZZZI')
    s4 = (Z_stabs[4][1:]=='IIIIIIZZZ')
    s_total = (s0 and s1 and s2 and s3 and s4)
    if not(s_total):
        raise NameError('The stabilizers are not in standard form.')
    else:
        signs = [Z_stabs[i][0] for i in range(len(Z_stabs))]
        signs_bin = []
        for sign in signs:
            if sign=='+':  signs_bin += [0]
            elif sign=='-':  signs_bin += [1]
        sign0 = (signs_bin[0] + signs_bin[2])%2 
        sign1 = signs_bin[1]
        sign2 = (signs_bin[2] + signs_bin[3])%2 
        sign3 = (signs_bin[3] + signs_bin[4])%2 
        logsign = signs_bin[4]
        on_codespace = (sign0==0 and sign1==0 and sign2==0 and sign3==0)
        log_err = (logsign==1)
        if (on_codespace and log_err):
            return True
        else:
            return False



def update_errors_anc(current_errors, stab):
    '''
    not error on ancilla, but errors that propagated
    from ancilla to data (hook errors)
    In this new implementation we don't need the input
    error, which was present in the old one.
    This is not the most elegant solution, but it works.
    MGA 6/29/16.
    
    Notice that this only works for weight-4 stabilizers.
    '''
    #print 'stab =', stab
    l = len(stab)   # should be 7 for Steane and 5 for 5qubit
    qubits_indexes = [i for i in range(l) if stab[i]!='I']
    qubits_to_correct = qubits_indexes[:2]
    for q in qubits_to_correct:
        error = stab[q]
        if current_errors[q] == 'I':
            current_errors[q] = error

        elif current_errors[q] == 'X':
            if error == 'X':
                current_errors[q] = 'I'
            elif error == 'Y':
                current_errors[q] = 'Z'
            elif error == 'Z':
                current_errors[q] = 'Y'

        elif current_errors[q] == 'Y':
            if error == 'X':
                current_errors[q] = 'Z'
            elif error == 'Y':
                current_errors[q] = 'I'
            elif error == 'Z':
                current_errors[q] = 'X'

        elif current_errors[q] == 'Z':
            if error == 'X':
                current_errors[q] = 'Y'
            elif error == 'Y':
                current_errors[q] = 'X'
            elif error == 'Z':
                current_errors[q] = 'I'

    return current_errors



def prepare_stabs_Steane(init_state, chp_loc='./chp_extended',
                         CHP_IO_files=True):
    '''
    Prepares initial stabilizers and destabilizers
    specifically for Steane
    '''

    try:
        stabs = steane.Code.stabilizer_logical_CHP[init_state]
        destabs = steane.Code.destabilizer_logical_CHP[init_state]
    except KeyError:
        if init_state == '+Y':
            init_circ = cor.Steane_Correct.encoded_plus_i_Steane()
            stabs, destabs = create_stabs_and_destabs(init_circ,
                                                      chp_loc,
                                                      CHP_IO_files)
        else:
            print 'Sorry, we cannot prepare the state %s for the Steane code.' %init_state
    
    return stabs, destabs
    

    
def prepare_stabs_Cross(init_state, chp_loc='./chp_extended', CHP_IO_files=False):
    '''
    Prepares initial stabilizers and destabilizers
    specifically for Cross code.  
    '''
    stabilizers = cross.Code.stabilizer_CHP[:]
    destabilizers = cross.Code.destabilizer_CHP[:]
    try:
        stabilizers += [cross.Code.stabilizer_logical_CHP[init_state]]
        destabilizers += [cross.Code.destabilizer_logical_CHP[init_state]]
    except KeyError:
        print 'Sorry, we cannot prepare the state %s for the Cross code.' %init_state

    return stabilizers, destabilizers



def prepare_stabs_5qubit(init_state, chp_loc='./chp_extended', CHP_IO_files=False):
    '''
    Prepares initial stabilizers and destabilizers
    for the 5-qubit code.
    '''
    stabilizers = fivequbit.Code.stabilizer_CHP[:]
    destabilizers = fivequbit.Code.destabilizer_CHP[:]
    try:
        stabilizers += [fivequbit.Code.stabilizer_logical_CHP[init_state]]
        destabilizers += [fivequbit.Code.destabilizer_logical_CHP[init_state]]
    except KeyError:
        print 'Sorry, we cannot prepare the state %s for the fivequbit code.' %init_state

    return stabilizers, destabilizers



def prepare_stabs_surface17(init_state, chp_loc='./chp_extended', CHP_IO_files=False):
    '''
    Prepares initial stabilizers and destabilizers
    for the 5-qubit code.
    '''
    if init_state != '+Z':
        print 'Only +Z for now'
        return
    stabilizers = surf17.Code.stabilizer_CHP[:]
    destabilizers = surf17.Code.destabilizer_CHP[:]
    #try:
    #    stabilizers += [fivequbit.Code.stabilizer_logical_CHP[init_state]]
    #    destabilizers += [fivequbit.Code.destabilizer_logical_CHP[init_state]]
    #except KeyError:
    #    print 'Sorry, we cannot prepare the state %s for the fivequbit code.' %init_state

    return stabilizers, destabilizers



def create_stabs_and_destabs(circ, chp_location, CHP_IO_files=True):
    '''
    return the stabs and destabs of a particular circuit
    with no ancilla qubits or errors.
    '''
    n_data = len(circ.data_qubits())
    circ_chp = chper.Chper(circ=circ, num_d_qub=n_data, num_a_qub=0,
                   stabs=[], destabs=[], anc_stabs=[], anc_destabs=[],
                   states='None', input_output_files=CHP_IO_files)
    
    dic, final_stabs, final_destabs = circ_chp.run(chp_location)
    #print 'stabs at the end outside method =', circ_chp.stabs

    return final_stabs, final_destabs



#def stabilizer_syndrome_steane(bin_error):
#    '''
#    outputs the syndrome based on the bin_error
#    '''
#    H = np.array([[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]])
#    bin_error = np.array(bin_error)
#    syn = np.dot(bin_error,H)
#    return syn


def measurement_parity_steane(meas_outcomes):
    """meas_outcomes needs to be a list of seven bits."""
    parity = sum(meas_outcomes)%2
    H = np.array([[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]])
    meas = np.array(meas_outcomes)
    s = np.dot(meas,H)
    for bit in s:
        if bit%2 != 0:
            return 1-parity   #, 'error'
    return parity   #, 'no error'
            
    
def measurement_interpretation_steane(dic):
    num_logical_meas = len(dic)/7
    #print dic
    #print dic[0]
    results = []
    for i in range(num_logical_meas):
        meas_outcomes = [dic[7*i + j][0] for j in range(7)]
        result = measurement_parity_steane(meas_outcomes)
        results += [result]
    return results



def FT_state_Steane(state, chp_location, add_error=False, error_info=None,
                    CHP_IO_files=False):
    '''
    prepare a FT logical 0 or logical + in the Steane code.
    state:  '0' or '+'
    
    (1) creates a circuit to generate FT a logical |0> or |+> in the Steane
    code
    (2) if add_error==True, adds errors to the circuit according to error_info
    (3) runs CHP as many times as necessary to obtain an even parity in
    the ancilla
    (4) outputs the number of tries, the stabilizers, and the destabilizers
    of the final state  
    '''

    # To make it faster we will assume that add_error is always True.
        
    if state == 'Z':  
        FT_func = cor.Steane_Correct.FT_encoded_zero_Steane
    elif state == 'X':
        FT_func = cor.Steane_Correct.FT_encoded_plus_Steane
    else:
        print 'Only 0 or +.  Sorry!'
        sys.exit(1)
    
    n_d_q = 7
    n_a_q = 7
    n_tries = 0
    parity = 1
    while parity == 1:
        n_tries += 1
        circ, meas_gates = FT_func()
        if add_error:
            error.add_error(circ, error_info)
        #brow.from_circuit(circ, True)

        run_circ = False
        for gate in circ.gates:
            #print gate.gate_name
            #t.sleep(1)
            if gate.is_error:
                #print 'gate is error'
                run_circ = True
                break           
        
        #print 'I broke off'
        #t.sleep(5)     

        # if it turns out that we don't add any errors to the circuit
        # (as might happen when the error rate is low), we don't even
        # care to run the circuit.  We just return the correct list
        # of stabilizers and destabilizers
        if not run_circ:
            final_stabs, final_destabs = steane.Code.stabs_dict[state]
            
            return n_tries, final_stabs, final_destabs, run_circ
    

        circ_chp = chper.Chper(circ=circ, num_d_qub=n_d_q, num_a_qub=n_a_q,
                       stabs=[], destabs=[], anc_stabs=[],
                       anc_destabs=[], states='None',
                       input_output_files=CHP_IO_files) 
        #print 'Before running circuit ...'
        #t.sleep(5)
        #print 'stabs =', circ_chp.stabs
        #t.sleep(5)
        dic, final_stabs, final_destabs = circ_chp.run(chp_location)
        keys = dic.keys()
        keys.sort()
        meas_outcomes = [dic[key][0] for key in keys]
        #print 'dic =', dic
        parity = steane.Code.parity_check(meas_outcomes)
        #print 'parity =', parity

    
    return n_tries, final_stabs, final_destabs, run_circ



def correct_one_error_SteaneEC(error_type, chp_location, initial_stabs, 
                  initial_destabs, add_error=False, error_info=None, 
                  CHP_IO_files=False, initial_trans=False):
    '''
    error can be either 'X' or 'Z'.  
    (1) If error == 'X', then prepares a FT |+> state, and uses it to 
    detect and correct X errors on the data.
    (2) If error == 'Z', then prepares a FT |0> state, and uses it to 
    detect and correct Z errors on the data.
    ''' 
    n_d_q = 7
    
    # prepare FT state
    n_tries, anc_stabs, anc_destabs, run_circ1 = FT_state_Steane(error_type, 
                            chp_location, add_error, 
                            error_info, CHP_IO_files)
    #print 'n tries =', n_tries, 'run_circ1 =', run_circ1
    #return n_tries, anc_stabs, anc_destabs, run_circ1  
    #sys.exit(0)

    # prepare error detection circuit
    circ, meas_gates = cor.Steane_Correct.detect_errors(error_type, n_d_q,
                                                        initial_trans)
    #brow.from_circuit(circ, True)
    #sys.exit(0)

    if add_error:
        error.add_error(circ, error_info)
    
    # run error detection circuit
    circ_chp = chper.Chper(circ=circ, num_d_qub=n_d_q, num_a_qub=n_d_q,
                   stabs=initial_stabs, destabs=initial_destabs,
                   anc_stabs=anc_stabs, anc_destabs=anc_destabs,
                   input_output_files=CHP_IO_files) 

    dic, final_stabs, final_destabs = circ_chp.run(chp_location)
    keys = dic.keys()
    keys.sort()
    meas_outcomes = [dic[key][0] for key in keys]
    correction = steane.Code.decode_syndrome_Steane_EC(meas_outcomes)
    if 'E' in correction:
        correction = [error_type if x=='E' else x for x in correction]
        corr_states =  update_stabs(final_stabs, final_destabs, correction)
        final_stabs, final_destabs = corr_states
    
    return n_tries, final_stabs, final_destabs, run_circ1   


    
def run_SteaneEC(chp_location, initi_stabs, initi_destabs, add_error=False,
                 error_info=None, CHP_IO_files=False, initial_trans=False):
    '''
    runs one round of Steane EC
    (1) Prepare FT |+>.  Use it to detect and correct X errors.
    (2) Prepare FT |0>.  Use it to detect and correct Z errors.
    Return final stabilizers and destabilizers and number of times
    we had to repeat the preparation of the FT |+> and |0> logical states.
    '''
    num_tries = {}
    run_circ = {}
    stabs = initi_stabs[:]
    destabs = initi_destabs[:]
    stab_index = 0

    for error_type in ['X', 'Z']:
        
        # only add the initial transversal gate at the very beginning
        # of the circuit.
        if stab_index == 1:  initial_trans = False

        results = correct_one_error_SteaneEC(error_type, chp_location, 
                                stabs, destabs, add_error, error_info, 
                                CHP_IO_files, initial_trans)
        num_tries[error_type], stabs, destabs, run_circ[error_type] = results

        stab_index += 1

    return num_tries, stabs, destabs, run_circ



def run_KnillEC(chp_location, initial_stabs, initial_destabs, add_error=False,
                error_info=None, CHP_IO_files=False, initial_trans=False,
                order='reverse'):
    '''
    Steane code with Knill correction
    '''
    n_tries = {}
    run_circ = {}
    n_d_q = 7
    
    # prepare logical |+> state
    n_tries['X'], plus_stabs, plus_destabs, run_circ['X'] = FT_state_Steane(
                        'X', chp_location, add_error,
                        error_info, CHP_IO_files)

    #print '|+> stabs =', plus_stabs

    
    # prepare logical |0> state
    n_tries['Z'], zero_stabs, zero_destabs, run_circ['Z'] = FT_state_Steane(
                        'Z', chp_location, add_error,
                        error_info, CHP_IO_files)
    #print '|0> stabs =', zero_stabs    


    if order == 'normal':
        # combine data with ancilla stabilizers and destabilizers
        stabs, destabs = combine_stabilizers(initial_stabs, initial_destabs,
                             plus_stabs, plus_destabs)
        stabs, destabs = combine_stabilizers(stabs, destabs,
                                 zero_stabs, zero_destabs)  
    
        final_stab_init_i = 2*n_d_q


    elif order == 'reverse':
        # combine data with ancilla stabilizers and destabilizers
        stabs, destabs = combine_stabilizers(zero_stabs, zero_destabs,
                             plus_stabs, plus_destabs)
        stabs, destabs = combine_stabilizers(stabs, destabs,
                                 initial_stabs, initial_destabs)    

        #print 'len(first stab) =', len(stabs[0])
        #print 'len(last stab) =', len(stabs[-1])
        #print 'zero stabs =', zero_stabs
        #print 'plus stabs =', plus_stabs

        final_stab_init_i = 0

    # prepare error detection circuit
    circ, meas_gates = cor.Knill_Correct.detect_errors(n_d_q, order,
                                                       initial_trans)
    #brow.from_circuit(circ, True)
        
    if add_error:
        error.add_error(circ, error_info)
    
    # run error detection circuit
    circ_chp = chper.Chper(circ=circ, num_d_qub=n_d_q, num_a_qub=2*n_d_q,
                   stabs=stabs, destabs=destabs,
                   anc_stabs=[], anc_destabs=[],
                   input_output_files=CHP_IO_files,
                   final_stab_init_i=final_stab_init_i) 
    
    dic, final_stabs, final_destabs = circ_chp.run(chp_location)    
    keys = dic.keys()
    keys.sort()
        
    Z_meas_outcomes = [dic[key][0] for key in keys[:n_d_q]]
    Z_parity = steane.Code.parity_meas_Steane_EC(Z_meas_outcomes)
    
    # apply X correction (which is obtained from Z_parity)
    if Z_parity == 1:
        X_corr_states = update_stabs(final_stabs, final_destabs, 
                         ['X' for i in range(n_d_q)])
        final_stabs, final_destabs = X_corr_states

    X_meas_outcomes = [dic[key][0] for key in keys[n_d_q:]]
    X_parity = steane.Code.parity_meas_Steane_EC(X_meas_outcomes)

    # apply Z_correction (which is obtained from X_parity)
    if X_parity == 1:
        Z_corr_states = update_stabs(final_stabs, final_destabs,
                         ['Z' for i in range(n_d_q)])
        final_stabs, final_destabs = Z_corr_states

    return n_tries, final_stabs, final_destabs, run_circ





def surface17_prep(chp_location, Is_after2q, error_info, decoder='lookuptable'):
    '''
    Simulates the prep of a logical |0> on the surface-17 code.
    '''
    
    # 1. Prepare all |0>.  Should be function, instead of doing it manually
    # but I'm in a rush.
    surf_stabs = ['+ZIIIIIIII',
                  '+IZIIIIIII',
                  '+IIZIIIIII',
                  '+IIIZIIIII',
                  '+IIIIZIIII',
                  '+IIIIIZIII',
                  '+IIIIIIZII',
                  '+IIIIIIIZI',
                  '+IIIIIIIIZ']
    
    surf_destabs = ['+XIIIIIIII',
                    '+IXIIIIIII',
                    '+IIXIIIIII',
                    '+IIIXIIIII',
                    '+IIIIXIIII',
                    '+IIIIIXIII',
                    '+IIIIIIXII',
                    '+IIIIIIIXI',
                    '+IIIIIIIIX']

    # 2. Run the circuit with the measurement of the X stabilizers.
    initial_I = True
    out_dict = run_surface17EC(chp_location, surf_stabs, surf_destabs,
                               Is_after2q, True, error_info, 0, [(0,0)],
                               'old', False, False, initial_I,
                               False, True, False, 'X')
   
    err_dict, run_dict = out_dict
    synd_Xstabs0 = run_dict['syndrome'][:]
    surf_stabs = run_dict['stabs'][:]
    surf_destabs = run_dict['destabs'][:]

    #print synd_Xstabs0
    #print surf_stabs
    
    # 3. Decide to stop if the X stabilizers were all 0 (1/16 chance)
    if synd_Xstabs0 != '0000':
        # 4. If the X stabilizers were not 0000, measure all stabs 1 time.
        out_dict = run_surface17EC(chp_location, surf_stabs, surf_destabs,
                                   Is_after2q, True, error_info, 0, [(0,0)],
                                   'old', False, False, False,
                                   False, True, False, 'all')
        err_dict, run_dict = out_dict
        synd_Zstabs1 = run_dict['syndrome'][:4]
        synd_Xstabs1 = run_dict['syndrome'][4:]
        surf_stabs = run_dict['stabs'][:]
        surf_destabs = run_dict['destabs'][:]
        
        #print synd_Xstabs1
        #print synd_Zstabs1
        #print surf_stabs

        if ((synd_Xstabs1 == synd_Xstabs0) and (synd_Zstabs1 == '0000')):
            corrZ = ['I' for i in range(9)]
            for i in surf17.Code.lookuptable['Xstabs'][synd_Xstabs1]:
                corrZ[i] = 'Z'
            corr_state = update_stabs(surf_stabs, 
                                      surf_destabs,
                                      ''.join(corrZ))
            surf_stabs = corr_state[0][:]
            surf_destabs = corr_state[1][:] 
           
        else:
            out_dict = run_surface17EC(chp_location, surf_stabs, surf_destabs,
                                       Is_after2q, True, error_info, 0, [(0,0)],
                                       'old', False, '', False, False,
                                       False, decoder, '', True, False, 'all')
            err_dict, run_dict = out_dict
            synd_Zstabs2 = run_dict['syndrome'][:4]
            synd_Xstabs2 = run_dict['syndrome'][4:]
            surf_stabs = run_dict['stabs'][:]
            surf_destabs = run_dict['destabs'][:]
            
            # Correct Z errors
            corrZ = ['I' for i in range(9)]
            for i in surf17.Code.lookuptable['Xstabs'][synd_Xstabs2]:
                corrZ[i] = 'Z'
            corr_state = update_stabs(surf_stabs, 
                                      surf_destabs,
                                      ''.join(corrZ))
            surf_stabs = corr_state[0][:]
            surf_destabs = corr_state[1][:] 
             
            # Correct X errors
            corrX = ['I' for i in range(9)]
            for i in surf17.Code.lookuptable['Zstabs'][synd_Zstabs2]:
                corrX[i] = 'X'
            corr_state = update_stabs(surf_stabs, 
                                      surf_destabs,
                                      ''.join(corrX))
            surf_stabs = corr_state[0][:]
            surf_destabs = corr_state[1][:] 
        
    
    # 5. Do perfect EC and determine if logical error    
              
    #surf_stabs[4] = surf_stabs[4].replace('+', '-')
    #surf_stabs[6] = surf_stabs[6].replace('+', '-')
    #surf_stabs[7] = surf_stabs[7].replace('+', '-')
    #surf_stabs[8] = surf_stabs[8].replace('+', '-')
 
    #print 'Doing perfect EC'
    surf_stabs, surf_destabs, synd = run_surface17EC(chp_location, surf_stabs, 
                                                     surf_destabs, Is_after2q, 
                                                     False, '', 0, [(0,0)],
                                                     'old', False, '', True)
       
    synd_Zstabs3 = synd[:4]
    synd_Xstabs3 = synd[4:]
        
    # Correct Z errors
    corrZ = ['I' for i in range(9)]
    for i in surf17.Code.lookuptable['Xstabs'][synd_Xstabs3]:
        corrZ[i] = 'Z'
    corr_state = update_stabs(surf_stabs, 
                              surf_destabs,
                              ''.join(corrZ))
    surf_stabs = corr_state[0][:]
    surf_destabs = corr_state[1][:] 
             
    # Correct X errors
    corrX = ['I' for i in range(9)]
    for i in surf17.Code.lookuptable['Zstabs'][synd_Zstabs3]:
        corrX[i] = 'X'
    corr_state = update_stabs(surf_stabs, 
                              surf_destabs,
                              ''.join(corrX))
    surf_stabs = corr_state[0][:]
    surf_destabs = corr_state[1][:] 

    #print surf_stabs 
        
    # 6. Determine if a logical error happened
    # This ugly, quick way to do it only works
    # it the state we're trying to prepare is |0>.
    # We just check that we are not |1>.
    if surf_stabs[-1][1:] != 'IIIIIIZZZ':
        raise NameError('Problems with stabs')
    else:
        if surf_stabs[-1][0] == '+':
            return 0
        else:
            return 1



def repeat_surface17_prep(n_times, chp_location, Is_after2q, error_info, 
                          decoder='lookuptable', output_folder='',
                          stream_index=0):
    '''
    '''
    n_error = 0
    n_print = 25000
    for n_run in range(n_times):
        n_error += surface17_prep(chp_location, Is_after2q, error_info,
                                  decoder)
        
        if n_run%n_print == 0:
            print '(%i) n_run = %i ; n_error = %i' %(stream_index, n_run, n_error) 
    

    out_dict = {'n_runs': n_times, 'n_errors': n_error}    
    out_filename = output_folder + str(stream_index) + '.json'
    out_file = open(out_filename, 'w')
    json.dump(out_dict, out_file, indent=4, separators=(',', ':'),
                  sort_keys=True)
    out_file.close()

    return n_error    



def run_surface17_until_failure(chp_location, decoder, Is_after2q, Is_after1q,
                                MS_heating, Stark, qubit_assignment, time_MS,
                                time_1q, output_folder, stream_index, error_info, 
                                init_state='+Z', initial_I=False, initial_trans=False,
                                ion_trap=False, corrtable=None, init_vect=None,
                                gamma=None, delta=None):
    '''
    For now, always slow sampler
    stream_index:  for parallel purposes
    ion_trap:  False if running abstract circuit.  True if running ion trap.
    '''
    surf_stabs, surf_destabs = prepare_stabs_surface17(init_state)    

    total_error_dict = {}
    logical_error = False
    n_runs = 0

    
    if decoder == 'max_like_onestep':

        err_vec = np.load(init_vec)
        Gamma_mat = np.load(gamma)
        Del_mat = np.load(delta)

        while(not(logical_error)):

            out_dict = run_surface17EC(chp_location, surf_stabs, surf_destabs,
                                       Is_after2q, Is_after1q, MS_heating,
                                       Stark, qubit_assignment, time_MS, time_1q,
                                       True, error_info, 0, [(0,0)],
                                       'old', False, False, False,
                                       False, True, False, 'all', ion_trap)

            if out_dict != None:
                err_dict, run_dict = out_dict
                total_error_dict[n_runs] = err_dict
                synd_Zstabs = run_dict['syndrome'][:4]
                synd_Xstabs = run_dict['syndrome'][4:]
                synd_stabs = synd_Zstabs + synd_Xstabs
                surf_stabs = run_dict['stabs'][:]
                surf_destabs = run_dict['destabs'][:]
            else:
                synd_stabs = '00000000'

            post_synd, err_vec = MLD.MLDecoder.decode_onestep(err_vec,
                                                              Gamma_mat,
                                                              Del_mat,
                                                              synd_stabs)

            #corr_dict = MLD.MLDecoder.Tilted17SC_DPcorr_MLD(post_synd)
            corr_dict = MLD.MLDecoder.correct_from_dict(post_synd,corrtable)

            if corr_dict['X corr'] != '000000000':
                corr_X = [ int(i) for i in corr_dict['X corr']]
                corr_X = ['X' if x==1 else 'I' for x in corr_X]

                corr_state = update_stabs(surf_stabs,
                                          surf_destabs,
                                          ''.join(corr_X))
                surf_stabs = corr_state[0][:]
                surf_destabs = corr_state[1][:]

            if corr_dict['Z corr'] != '000000000':
                corr_Z = [ int(i) for i in corr_dict['Z corr']]
                corr_Z = ['Z' if z==1 else 'I' for z in corr_Z]

                corr_state = update_stabs(surf_stabs,
                                          surf_destabs,
                                          ''.join(corr_Z))
                surf_stabs = corr_state[0][:]
                surf_destabs = corr_state[1][:]


            logical_error = logical_error_surface17(surf_stabs,
                                                    init_state)

            n_runs += 1


    
    
    elif decoder == 'max_like_twostep':

        err_vec = np.load(init_vec)
        Gamma_mat = np.load(gamma)

        while(not(logical_error)):
            

            out_dict1 = run_surface17EC(chp_location, surf_stabs, surf_destabs,
                                       Is_after2q, Is_after1q, MS_heating,
                                       Stark, qubit_assignment, time_MS, time_1q,
                                       True, error_info, 0, [(0,0)],
                                       'old', False, False, False,
                                       False, True, False, 'all', ion_trap)

            if out_dict1 != None:
                err_dict, run_dict = out_dict1
                total_error_dict[n_runs] = err_dict 
                synd_stabs1 = run_dict['syndrome'][:]
                surf_stabs = run_dict['stabs'][:]
                surf_destabs = run_dict['destabs'][:]
            

            out_dict2 = run_surface17EC(chp_location, surf_stabs, surf_destabs,
                                       Is_after2q, Is_after1q, MS_heating,
                                       Stark, qubit_assignment, time_MS, time_1q,
                                       True, error_info, 0, [(0,0)],
                                       'old', False, False, False,
                                       False, True, False, 'all', ion_trap)
            

            if out_dict2 != None:
                err_dict, run_dict = out_dict2
                total_error_dict[n_runs+1] = err_dict
                synd_stabs2 = run_dict['syndrome'][:]
                surf_stabs = run_dict['stabs'][:]
                surf_destabs = run_dict['destabs'][:]

            #process two rounds of syndrome measurement, bitwise AND of syndrome strings
            if out_dict1 != None and out_dict2 != None:
                #print 'syndrome step 1: '+str(synd_stabs1)
                #print 'syndrome step 2: '+str(synd_stabs2)
                synd_twostep = np.array([int(i) for i in synd_stabs1]) & np.array([int(j) for j in synd_stabs2])
                synd_twostep = ''.join(str(i) for i in synd_twostep)
                synd_twostep = synd_twostep[:4] + synd_twostep[4:]
                #print 'syndrome both: '+str(synd_twostep)
            else:
                synd_twostep = '00000000' #still need syndrome to update likelihood vector!

            post_synd, err_vec = MLD.MLDecoder.decode_twostep(err_vec, Gamma_mat, synd_twostep)

            #corr_dict = MLD.MLDecoder.Tilted17SC_DPcorr_MLD(post_synd)
            corr_dict = MLD.MLDecoder.correct_from_dict(post_synd,corrtable)

            if corr_dict['X corr'] != '000000000':
                corr_X = [int(i) for i in corr_dict['X corr']]
                corr_X = ['X' if x==1 else 'I' for x in corr_X]

                corr_state = update_stabs(surf_stabs, surf_destabs,''.join(corr_X))
                surf_stabs = corr_state[0][:]
                surf_destabs = corr_state[1][:]

            if corr_dict['Z corr'] != '000000000':
                corr_Z = [int(i) for i in corr_dict['Z corr']]
                corr_Z = ['Z' if z==1 else 'I' for z in corr_Z]

                corr_state = update_stabs(surf_stabs, surf_destabs,''.join(corr_Z))
                surf_stabs = corr_state[0][:]
                surf_destabs[1][:]

            logical_error = logical_error_surface17(surf_stabs, init_state)

            n_runs+=2




    elif decoder == 'lookuptable':
        
        Xerror_det, Zerror_det = False, False

        while(not(logical_error)):
           
            if (Xerror_det or Zerror_det):  previous_errors = True
            else:  previous_errors = False

            ######### just for debugging
            #if n_runs == 0:  add_errors = True
            #else:            add_errors = False
            add_errors = False
            ############################

            

            out_dict = run_surface17EC(chp_location, surf_stabs, surf_destabs,
                                       Is_after2q, Is_after1q, MS_heating,
                                       Stark, qubit_assignment, time_MS, time_1q,
                                       True, error_info, 0, [(0,0)],
                                       'old', False, False, False,
                                       False, previous_errors,
                                       add_errors, 'all', ion_trap)
            
 
            if out_dict != None:
                err_dict, run_dict = out_dict
                total_error_dict[n_runs] = [err_dict, run_dict['syndrome']]
                synd_Zstabs = run_dict['syndrome'][:4]
                synd_Xstabs = run_dict['syndrome'][4:]
                surf_stabs = run_dict['stabs'][:]
                surf_destabs = run_dict['destabs'][:]
            
                print n_runs  
                print 'Zstabs =', synd_Zstabs
                print Xerror_det
                print 'Xstabs =', synd_Xstabs 
                print Zerror_det
                print surf_stabs
    
                # if no error is detected
                if (synd_Zstabs=='0000' and synd_Xstabs=='0000'):
                    Xerror_det, Zerror_det = False, False
                    #print n_runs
                    #print surf_stabs
                    logical_error = logical_error_surface17(surf_stabs,
                                                            init_state)

                # if only a Z error is detected
                elif (synd_Zstabs=='0000' and synd_Xstabs!='0000'):
                    Xerror_det = False
                    if Zerror_det:
                        # correct final state with corrZ
                        corrZ = ['I' for i in range(9)]
                        for i in surf17.Code.lookuptable['Xstabs'][synd_Xstabs]:
                            corrZ[i] = 'Z'
                        corr_state = update_stabs(surf_stabs, 
                                                  surf_destabs,
                                                  ''.join(corrZ))
                        surf_stabs = corr_state[0][:]
                        surf_destabs = corr_state[1][:] 
                        #print corrZ
                        #print n_runs
                        #sys.exit(0)
                    else:
                        Zerror_det = True
                                      
         
                # if only an X error is detected 
                elif (synd_Zstabs!='0000' and synd_Xstabs=='0000'):
                    Zerror_det = False
                    if Xerror_det:
                        # correct final state with corrX
                        corrX = ['I' for i in range(9)]
                        for i in surf17.Code.lookuptable['Zstabs'][synd_Zstabs]:
                            corrX[i] = 'X'
                        corr_state = update_stabs(surf_stabs, 
                                                  surf_destabs,
                                                  ''.join(corrX))
                        surf_stabs = corr_state[0][:]
                        surf_destabs = corr_state[1][:] 
                        #print corrX
                        #print n_runs
                        #sys.exit(0)
                    else:
                        Xerror_det = True

                # if both X and Z errors are detected
                else:
                    if Xerror_det:
                        # correct final state with corrX
                        corrX = ['I' for i in range(9)]
                        for i in surf17.Code.lookuptable['Zstabs'][synd_Zstabs]:
                            corrX[i] = 'X' 
                        corr_state = update_stabs(surf_stabs, 
                                                  surf_destabs,
                                                  ''.join(corrX))
                        surf_stabs = corr_state[0][:]
                        surf_destabs = corr_state[1][:] 
                        #print corrX
                        #print n_runs
                    else:
                        Xerror_det = True
                    
                    if Zerror_det:
                        # correct final state with corrZ
                        corrZ = ['I' for i in range(9)]
                        for i in surf17.Code.lookuptable['Xstabs'][synd_Xstabs]:
                            corrZ[i] = 'Z' 
                        corr_state = update_stabs(surf_stabs, 
                                                  surf_destabs,
                                                  ''.join(corrZ))
                        surf_stabs = corr_state[0][:]
                        surf_destabs = corr_state[1][:] 
                        #print corrZ
                        #print n_runs
                        #sys.exit(0)
                    else:
                        Zerror_det = True
                    

                # add something to the output_dictionary
           
            n_runs += 1
            #if n_runs == 10:  break   
 
    # write output_dict to hard drive
    return n_runs, total_error_dict



def run_several_sim_until_failure(chp_location, decoder, Is_after2q, Is_after1q,
                                  MS_heating, Stark, qubit_assignment, time_MS,
                                  time_1q, output_folder, stream_index, error_info, 
                                  n_sim, init_state='+Z', initial_I=False, 
                                  initial_trans=False, ion_trap=False,
                                  corrtable=None, init_vec=None, gamma=None, delta=None,
                                  save_error_dict=False):
    '''
    '''
    output_dict = {}

    if save_error_dict:
        err_folder = output_folder + 'error_dicts/'
        if not os.path.exists(err_folder):
            os.makedirs(err_folder)
    else:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)


    for sim_i in range(n_sim):
        n_runs, error_dict = run_surface17_until_failure(chp_location, decoder, Is_after2q,
                                                         Is_after1q, MS_heating, Stark,
                                                         qubit_assignment, time_MS, time_1q,
                                                         output_folder, stream_index,
                                                         error_info, init_state, initial_I,
                                                         initial_trans, ion_trap, corrtable,
                                                         init_vec, gamma, delta) 
        if type(stream_index) is int:
            output_dict[stream_index*n_sim + sim_i] = n_runs
            if save_error_dict:
                err_filename = err_folder + str(stream_index*n_sim + sim_i) + '.json'
        else:
            output_dict[sim_i] = n_runs       
            if save_error_dict:
                err_filename = err_folder + stream_index + '_' + str(sim_i) + '.json'
        
        if save_error_dict:
            err_file = open(err_filename, 'w')
            json.dump(error_dict, err_file, indent=4, separators=(',', ':'),
                      sort_keys=True)
            err_file.close()

    output_filename = output_folder + str(stream_index) + '.json'
    output_file = open(output_filename, 'w')
    json.dump(output_dict, output_file, indent=4, separators=(',', ':'),
              sort_keys=True)
    output_file.close()

    return



# def run_simulation_EC(n_runs, initi_stabs, initi_destabs, chp_location, add_error, error_info, n_errors, gate_indices,
#                       CHP_IO_files=False, first_run_index=0,
#                       kind, initial_trans=False, initial_state = '+Z', sampling = 'old'):
def run_simulation_EC(n_runs, initi_stabs, initi_destabs, chp_location, add_error, 
                      error_info, n_errors, gate_indices, CHP_IO_files, first_run_index,
                      kind, initial_I, initial_trans, initial_state, sampling, Is_after2q):
    '''
    '''
    output_dict = {}

    #print 'add error =', add_error
    #print 'error info =', error_info.dic
    #print 'n errors =', n_errors
    #print 'gate indices =', gate_indices[1]

    if kind[:7]=='Prepare':
       
        # first run the circuit without the errors to get the
        # perfect stabilizers. 
        perfect_stabs = FT_state_Steane(kind[7], chp_location)[1]

        for run in xrange(n_runs):
            results = FT_state_Steane(kind[7], chp_location, add_error,
                                      error_info, CHP_IO_files)
            num_tries, stabs, destabs, run_circ = results       
            if stabs != perfect_stabs:
                # do perfect EC to get rid of correctable errors.
                corr_stabs, corr_destabs = run_ShorEC(chp_location=chp_location,
                                                      initi_stabs=stabs,
                                                      initi_destabs=destabs,
                                                      Is_after2q=Is_after2q,
                                                      carry_run=True)
                res_dict = {'num_tries': num_tries,
                           'stabs':     stabs,
                           'destabs':   destabs,
                           'run_circ':  run_circ,
                           'corr_stabs': corr_stabs,
                           'corr_destabs': corr_destabs
                           }
                output_dict[run + first_run_index] = res_dict
            
    
    elif kind=='Shor' or kind=='Cross' or kind=='ShorEC_fivequbit':
        
        if kind=='Shor':  run_func = run_ShorEC
        elif kind=='Cross':  run_func = run_CrossEC
        elif kind=='ShorEC_fivequbit':  run_func = run_5qubit

        for run in xrange(n_runs):
            results = run_func(chp_location, initi_stabs, initi_destabs,
                               Is_after2q, add_error, error_info, 
                               CHP_IO_files, initial_state, False,
                               sampling, n_errors, gate_indices,
                               initial_I, initial_trans) 
            if results != None:
                stabs = results[1]['stabs']
                destabs = results[1]['destabs']
                # do perfect EC to get rid of correctable errors.
                corr_stabs, corr_destabs = run_func(chp_location=chp_location,
                                                    initi_stabs=stabs,
                                                    initi_destabs=destabs,
                                                    Is_after2q=Is_after2q,
                                                    carry_run=True) 

                results[1]['corr_stabs'] = corr_stabs
                results[1]['corr_destabs'] = corr_destabs
                output_dict[run + first_run_index] = results
    
    else:

        if kind=='Steane' or kind=='SteaneEC':   simulation_func = run_SteaneEC
        elif kind=='Knill' or kind=='KnillEC':  simulation_func = run_KnillEC

        perfect_stabs = simulation_func(chp_location, initi_stabs, initi_destabs)[1]
        
        for run in xrange(n_runs):
            results = simulation_func(chp_location, initi_stabs, initi_destabs,
                           add_error, error_info, CHP_IO_files, initial_trans)
            num_tries, stabs, destabs, run_circ = results
            if stabs != perfect_stabs:
                # do perfect EC to get rid of correctable errors.
                corr_stabs, corr_destabs = run_ShorEC(chp_location=chp_location,
                                                      initi_stabs=stabs,
                                                      initi_destabs=destabs,
                                                      carry_run=True)
                res_dict = {'num_tries': num_tries,
                           'stabs':     stabs,
                           'destabs':   destabs,
                           'run_circ':  run_circ,
                           'corr_stabs': corr_stabs,
                           'corr_destabs': corr_destabs
                           }
                output_dict[run + first_run_index] = res_dict

    
    return output_dict

# Alonzo's new functions #

def locate_error_gates(sub_circs, error_gate_names=['I']):
    error_gates = []
    for i in range(len(sub_circs)):
        #print 'running circuit %i' %i
        #subcirc_gates = sub_circs[i].gates[0].circuit_list[0].gates
        subcirc_gates = sub_circs[i].gates
        for j in range(len(subcirc_gates)):
            #print 'running gate %i' %j
            if subcirc_gates[j].gate_name in error_gate_names:
                error_gates += [(i,j)]

    return error_gates	

# calculate the total probability of all possible error subsets
def totalprob(ps,pt, ns, nt):
    result = 0
    for s in range(ns+1):
        for t in range(nt+1):
            result = result + scipy.misc.comb(ns, s, exact=False)*(ps**s)*((1-ps)**(ns-s))*scipy.misc.comb(nt,t,exact=False)*(pt**t)*((1-pt)**(nt-t))
    return result

# calculate the probability of a particular error subset (s,t)
def prob_for_subset(ps,pt,ns,nt,s,t):
    a = scipy.misc.comb(ns, s, exact=False)*(ps**s)*((1-ps)**(ns-s))*scipy.misc.comb(nt,t,exact=False)*(pt**t)*((1-pt)**(nt-t))
    return a


def prob_for_subset_general(n_gates_list, n_errors_list, n_ps):
    '''
    Generalized function for prob_for_subset
    calculates the probability of occurrence for a given general subset.
    n_gates_list:  a list with the number of gates for each kind.
    n_errors_list:  number of errors after each particular kind of gate
    '''
    a = 1.
    for i in range(len(n_gates_list)):
        n_g = n_gates_list[i]   # number of gates
        n_e = n_errors_list[i]  # number of errors
        p = n_ps[i]   # error p
        a *= scipy.misc.comb(n_g,n_e,exact=False)*(p**n_e)*((1-p)**(n_g-n_e))

    return a



def find_necessary_subsets(n_gates_list, n_ps, tolerance, all_subsets):
    '''
    For a given value of tolerance, it returns the necessary subsets.
    n_gates_list:  a list with the number of gates for each kind.
    n_ps:  a list with the error rate for each kind of gate.
    tolerance:  the probability of occurrence of the subsets not included.
    all_subsets:  all the subsets to include
    '''

    prob_list = []   # list of probabilities of occurrence of every subset
    for subset in all_subsets:
        prob_list += [prob_for_subset_general(n_gates_list, subset, n_ps)]    

    n_subsets = len(all_subsets)
    total_prob = sum(prob_list)
    #if total_prob < 1. - tolerance:
    #    print 'Not enough subsets.'
    #    print 'Total prob. of occurrence =', total_prob
    #    return [], [], n_subsets

    sorted_prob_list, sorted_indices = [], []
    stepwise_total_prob = 0.
    found_last_subset = False
    for i in range(n_subsets):
        max_prob = max(prob_list)
        sorted_prob_list += [max_prob]
        max_prob_index = prob_list.index(max_prob)
        sorted_indices += [max_prob_index]
        prob_list[max_prob_index] = -1.
   
        stepwise_total_prob += max_prob
        if (not found_last_subset) and (stepwise_total_prob >= 1.-tolerance):
            last_subset_i = i    
            found_last_subset = True

    if not found_last_subset:  last_subset_i = n_subsets

    return sorted_indices, sorted_prob_list, last_subset_i



def cardinality_subset(n_gates_list, n_errors_list):
    '''
    Calculates the cardinality (number of elements) of a given subset,
    specified by the list n_errors_list.
    n_gates_list:  a list with the number of gates for each kind.
    n_errors_list:  number of errors after each particular kind of gate
    '''
    a = 1.
    for i in range(len(n_gates_list)):
        n_g = n_gates_list[i]   # number of gates
        n_e = n_errors_list[i]  # number of errors
        a *= scipy.misc.comb(n_g,n_e,exact=False)

    return a


# given a tolerance value, find all error subsets that have probability larger than tolerance
def find_subsets(ps,pt,ns,nt,tol):
    tot = totalprob(ps,pt,ns,nt)
    stop_1 = 0
    stop_2 = 0
    stops_12 = []
    subsets = []
    
    for i in range(ns):
        if prob_for_subset(ps,pt,ns,nt,i,0)/tot > tol and prob_for_subset(ps,pt,ns,nt,i+1,0)/tot < tol:
            stop_1 = i
            break
    for i in range(stop_1 + 1):
        subsets.append((i,0))
    
    for j in range(nt):
        if prob_for_subset(ps,pt,ns,nt,0,j)/tot > tol and prob_for_subset(ps,pt,ns,nt,0,j+1)/tot < tol:
            stop_2 = j
            break
    for j in range(1,stop_2 + 1):
        subsets.append((0,j))
    
    for i in range(1,ns+1):
        for j in range(1, nt):
            if prob_for_subset(ps,pt,ns,nt,i,j)/tot > tol and prob_for_subset(ps,pt,ns,nt,i,j+1)/tot < tol:
                stops_12.append((i,j))
                break
            elif prob_for_subset(ps,pt,ns,nt,i,j)/tot < tol:
                break
    for pair in stops_12:
        for j in range(1,pair[1]+1):
            subsets.append((pair[0], j))

    subsets.remove((0,0))

    return subsets


# return a list of single qubit gates and a list of two qubit gates
def gates_list(circs, faulty_gates_names):
    single_qubit_gates = []
    two_qubit_gates = []
    high_w_gates = []
    gl = faulty_gates_names[:]

    #for i in range(len(circs)):
        #gate = circs.gates[i]
        #for j in range(len(gate.circuit_list[0].gates):
    #return

    for i in range(len(circs)):
        gates = circs[i].gates
        for j in range(len(gates)):
            if gates[j].gate_name in gl:
                if len(gates[j].qubits) == 1:
                    single_qubit_gates.append((i,j))
                elif len(gates[j].qubits) == 2:
                    #print gates[j].gate_name
                    two_qubit_gates.append((i,j))
                elif len(gates[j].qubits) > 2:
                    #print gates[j].gate_name
                    high_w_gates.append((i,j))

    return single_qubit_gates, two_qubit_gates, high_w_gates



def gates_list_general(circs, faulty_gates_names_grouped):
    '''
    More general than gates_list.
    faulty_gates_names_grouped:  a list of lists specifying how to group
    the faulty gates.  
    For example, [['IMS5', 'MS'], ['ImZ', 'ImX'], ['Ism'], ['Icool']]
    means that there are 4 groups.
    '''

    out_gates_lists = [[] for group in faulty_gates_names_grouped]
    
    for i in range(len(circs)):
        gates = circs[i].gates
        for j in range(len(gates)):
            for k in range(len(out_gates_lists)):
                if gates[j].gate_name in faulty_gates_names_grouped[k]:
                    out_gates_lists[k].append((i,j))

    return out_gates_lists



def gates_list_surface(flat_circ, faulty_gates_names):
    single_qubit_gates = []
    two_qubit_gates = []
    gates = flat_circ.gates
    gl = faulty_gates_names[:]

    for i in range(len(gates)):
        if gates[i].gate_name in gl:
            if len(gates[i].qubits) == 1:
                single_qubit_gates.append(i)
            elif len(gates[i].qubits) == 2:
                two_qubit_gates.append(i)

    return single_qubit_gates, two_qubit_gates



def create_EC_subcircs(operation, Is_after2q,
                       initial_I=True, initial_trans=False):
    '''
    Cross, Shor, 5-qubit for now
    '''
    
    redun = 3
    verify = False
    ancilla_parallel = True
    diVincenzo = True
    meas_errors = True

    # Cross code uses bare ancillae
    if operation == 'Cross':
        Cross_stabs = cross.Code.stabilizer_Colin[:]
        n_subcircs = 3
        EC_circ = cor.Bare_Correct.generate_rep_bare_meas(7,Cross_stabs, 
                                                          redun, initial_I,
                                                          meas_errors, 
                                                          Is_after2q,
                                                          initial_trans,
                                                          ancilla_parallel)
        #EC_circ = EC_circ.gates[0].circuit_list[0]
    
    # Steane and 5-qubit codes use 4-qubit cat states as ancillae
    else:
        if operation == 'Shor':
            initial_stabs = steane.Code.stabilizer[:]
            code = 'Steane' 
            n_subcircs = 6
        elif operation == 'ShorEC_fivequbit':
            initial_stabs = fivequbit.Code.stabilizer[:]
            code = '5qubit'
            n_subcircs = 3

        EC_circ = cor.Cat_Correct.cat_syndrome_4(initial_stabs,
                                                 redun,
                                                 verify,
                                                 ancilla_parallel,
                                                 diVincenzo,
                                                 initial_I,
                                                 initial_trans,
                                                 code,
                                                 meas_errors,
                                                 Is_after2q)    
        #EC_circ = EC_circ.gates[0].circuit_list[0]


    brow.from_circuit(EC_circ, True)
    sys.exit(0)

    subcircs = []
    for i in range(n_subcircs):
        subcirc = Circuit(gates=[EC_circ.gates[i]])
        subcirc.update_map()
        subcirc = subcirc.unpack()
        subcircs += [subcirc]
   
    return subcircs
 


def gates_list_for_operation(operation, faulty_gates_names,
                             Is_after2q, initial_I=True,
                             initial_trans=False):
    '''
    Builds the circuit from operation and returns 
    a list of 1-qubit and 2-qubit gates 
    '''
    subcircs = create_EC_subcircs(operation, Is_after2q,
                                  initial_I, initial_trans)
    gate_indices = gates_list(subcircs, faulty_gates_names)

    return gate_indices



def get_errors_dict(subcirc):
    '''
    '''
    carry_run = False
    subcirc_errors_dict = {}
    data_errors, anc_errors = [], []
    for gate in subcirc.gates:
        if gate.is_error:
            carry_run = True
            qubit = gate.qubits[0]
            gate_i = subcirc.gates.index(gate)
            after_gate = subcirc.gates[gate_i-1]
            if qubit.qubit_type == 'data':
                data_errors += [{'error': (gate.gate_name, qubit.qubit_id),
                                 'location': (after_gate.gate_name, gate_i)}]
            elif qubit.qubit_type == 'ancilla':
                anc_errors += [{'error': (gate.gate_name, qubit.qubit_id),
                                'location': (after_gate.gate_name, gate_i)}]

    if len(data_errors) > 0 or len(anc_errors) > 0:
        subcirc_errors_dict = {'d': data_errors, 'a': anc_errors}

    return subcirc_errors_dict, carry_run




def dict_for_error_model(error_model, p_1q, p_2q, p_meas, 
                         p_bath=0., heating_rate=0., 
                         Stark_rate=0., p_prep=0., 
                         p_sm=0., p_cool=0.,
                         p_cross=0., p_5q=0.):
    '''
    p_prep:  the error rate associated with a qubit state prep.
             We assume it's a bit-flip channel for Z prep and
             a phase-flip channel for X prep.
    p_sm:    the error rate associated with the idle time during 
             shuttling (s) and merging (m).  We assume it's only
             phase damping.
    p_cool:  the error rate associated with the idle time during
             the cooling.  Likewise, we assume it's only phase
             damping.
    p_cross: the error rate associated with the idle time during
             the crossing through the trap junction.  Likewise, we
             assume it's only phase damping.
    p_5q:  the error rate of the 5-qubit MS gate
    '''

    error_dict = {}
    error_dict['error_kind'] = 1
    error_dict['ImX'] = {'error_rate': p_meas, 'error_ratio': {'Z': 1}}
    error_dict['ImZ'] = {'error_rate': p_meas, 'error_ratio': {'X': 1}}

    if error_model == 'standard':
        
        Is_after_two_qubit = False
        Is_after_one_qubit = False
        # Depolarizing Pauli noise after 1-qubit gates
        one_qubit_gates = ['PrepareZ', 'PrepareZPlus', 'PrepareZMinus',
                           'PrepareX', 'PrepareXPlus', 'PrepareXMinus',  
                           'I', 'X', 'Y', 'Z', 'H', 'S']

        one_qubit_ratio = {'X': 1, 'Y': 1, 'Z': 1}
        one_qubit_dict = {
                    'error_rate': p_1q,
                    'error_ratio': one_qubit_ratio
                    }
    
        for g in one_qubit_gates:
            error_dict[g] = one_qubit_dict
        
        twoq_ratio = {}
        for prod in it.product('IXYZ', repeat=2):
            twoq_error = ''.join(prod)
            if twoq_error == 'II':
                continue
            twoq_ratio[twoq_error] = 1
        
        for g in ['CX','CY','CZ']:
            error_dict[g] = {'error_rate': p_2q, 'error_ratio': twoq_ratio} 


    elif error_model == 'ion_trap0_NN':
        # simple error model added as a quick approximation to realistic
        # ion-trap noise for the neural-network decoder
        # MGA: 2/25/2020        

        Is_after_two_qubit = False
        Is_after_one_qubit = False
        error_dict['PrepareXPlus'] = {'error_rate': p_meas, 'error_ratio': {'Z': 1}}
        error_dict['CX'] = {'error_rate': p_meas, 'error_ratio': {'ZX': 1}}
        error_dict['CY'] = {'error_rate': p_meas, 'error_ratio': {'ZY': 1}}
        error_dict['CZ'] = {'error_rate': p_meas, 'error_ratio': {'ZZ': 1}}
        
        error_dict['I_gate'] = {'error_rate': p_1q, 'error_ratio': {'Z': 1}}
        error_dict['I_idle'] = {'error_rate': p_2q, 'error_ratio': {'Z': 1}}



    elif error_model == 'ion_trap_simple':
        
        Is_after_two_qubit = True
        Is_after_one_qubit = False
        # Depolarizing Pauli noise after 1-qubit gates
        one_qubit_gates = ['PrepareZ', 'PrepareZPlus', 'PrepareZMinus',
                           'PrepareX', 'PrepareXPlus', 'PrepareXMinus',  
                           'I', 'X', 'Y', 'Z', 'H', 'S']

        one_qubit_ratio = {'X': 1, 'Y': 1, 'Z': 1}
        one_qubit_dict = {
                    'error_rate': p_1q,
                    'error_ratio': one_qubit_ratio
                    }
    
        for g in one_qubit_gates:
            error_dict[g] = one_qubit_dict
        
        # Over-rotation after CX, CY, or CZ.  Assume these are primitives.
        error_dict['CX'] = {'error_rate': p_2q, 'error_ratio': {'ZX': 1}}
        error_dict['CY'] = {'error_rate': p_2q, 'error_ratio': {'ZY': 1}}
        error_dict['CZ'] = {'error_rate': p_2q, 'error_ratio': {'ZZ': 1}}


    elif error_model == 'ion_trap0':
        # Assumes only gates are the primitive ones for ion traps:
        # preparations, measurements, 1-qubit rotations, MS gates.       
 
        Is_after_two_qubit = True
        Is_after_one_qubit = True

        # Over-rotations after 1-qubit gates
        error_dict['RX +'] = {'error_rate': p_1q, 'error_ratio': {'X': 1}}
        error_dict['RY +'] = {'error_rate': p_1q, 'error_ratio': {'Y': 1}}
        error_dict['RZ +'] = {'error_rate': p_1q, 'error_ratio': {'Z': 1}}
        error_dict['RX -'] = {'error_rate': p_1q, 'error_ratio': {'X': 1}}
        error_dict['RY -'] = {'error_rate': p_1q, 'error_ratio': {'Y': 1}}
        error_dict['RZ -'] = {'error_rate': p_1q, 'error_ratio': {'Z': 1}}

        # Depolarizing Pauli bath after preparations and I.
        one_qubit_gates = ['PrepareZ', 'PrepareZPlus', 'PrepareZMinus',
                           'PrepareX', 'PrepareXPlus', 'PrepareXMinus',  
                           'I']

        one_qubit_ratio = {'X': 1, 'Y': 1, 'Z': 1}
        one_qubit_dict = {
                    'error_rate': p_bath,
                    'error_ratio': one_qubit_ratio
                    }
    
        for g in one_qubit_gates:
            error_dict[g] = one_qubit_dict


        # Over-rotation after MS gates
        error_dict['MS'] = {'error_rate': p_2q, 'error_ratio': {'XX': 1}}
        
        # XX after II_heat.
        # II_heat is a 2-qubit I gate inserted after the MS.
        # The heating is not the error rate.
        # The error rate is obtained by multiplying the heating rate by the
        # duration of the MS gate, which in turns depends on the ions' separation.
        error_dict['II_heat'] = {'error_rate': heating_rate, 'error_ratio': {'XX': 1}}
            
        # Z after I_Stark
        # I_Stark is a 1-qubit I gate inserted after any 1-qubit rotation of MS gate
        # The error rate is obtained by multiplying the Stark_rate by the duration of
        # the gate.
        error_dict['I_stark'] = {'error_rate': Stark_rate, 'error_ratio': {'Z': 1}}

    
    elif error_model[:14] == 'ion_trap_eQual':
        
        Is_after_two_qubit = False
        Is_after_one_qubit = False

        # For now, no errors after 1-qubit gates to speed up the fast sampler
        # After preparations we add a flip (bit or phase)
        prep_gatesZ = ['PrepareZ', 'PrepareZPlus', 'PrepareZMinus']
        prep_dictZ = {'error_rate': p_prep, 'error_ratio': {'X': 1}}
        for g in prep_gatesZ:
            error_dict[g] = dict(prep_dictZ)
        
        prep_gatesX = ['PrepareX', 'PrepareXPlus', 'PrepareXMinus']
        prep_dictX = {'error_rate': p_prep, 'error_ratio': {'Z': 1}}
        for g in prep_gatesX:
            error_dict[g] = dict(prep_dictX)

        # All the prep gates and the measurements have the same error rate for
        # now, so we group them together.
        faulty_group1 = prep_gatesZ + prep_gatesX + ['ImZ', 'ImX']


        # Over-rotation after MS gates
        # For now the error after an MS gate will be the symmetric
        # depolarizing channel, not XX, in order to be on the same page
        # with arXiv:1705.02771.  Might change later on.
        
        # 2-qubit MS gate
        #error_dict['MS'] = {'error_rate': p_2q, 'error_ratio': {'XX': 1}}
        twoq_ratio = {}
        for prod in it.product('IXYZ', repeat=2):
            twoq_error = ''.join(prod)
            if twoq_error == 'II':
                continue
            twoq_ratio[twoq_error] = 1
        error_dict['MS'] = {'error_rate': p_2q, 'error_ratio': twoq_ratio} 

        # 5-qubit MS gate
        # Since we don't know yet how to compile a 5-qubit MS gate into
        # CHP, we are using this trick where we just use CXs and CZs and
        # afterwards add a gate name IMS5.
        # The error rate has to be twice p_2q, because when we measure a
        # stabilizer using 5-qubit MS gates, we really apply 2 of these.
        # Later on, we might have a different error rate for 2-qubit MS gates
        # and 5-qubit MS gates.
        #error_dict['IMS5'] = {'error_rate': 2*p_2q, 'error_ratio': 'XXXXX': 1}}
        fiveq_ratio = {}
        for prod in it.product('IXYZ', repeat=5):
            fiveq_error = ''.join(prod)
            if fiveq_error == 'IIIII':
                continue
            fiveq_ratio[fiveq_error] = 1

        
        if error_model[-1] == '2' or error_model[-1] == '3':
            # Second and thirdeQual error model
            
            # Errors after idle qubit and junction crossing.
            # We assume p_sm is the error rate associated with the idle qubits.
            error_dict['I_idle'] = {'error_rate': p_sm, 'error_ratio': {'Z': 1}}
            faulty_group3 = ['I_idle']
            error_dict['I_cross'] = {'error_rate': p_cross, 'error_ratio': {'Z': 1}}
            faulty_group4 = ['I_cross']
        
            # Depolarizing Pauli bath after 1-q rotations.
            one_qubit_gates = ['RX +', 'RX -', 'RY +', 'RY -', 'RZ +', 'RZ -']
            one_qubit_ratio = {'X': 1, 'Y': 1, 'Z': 1}
            one_qubit_dict = {
                    'error_rate': p_1q,
                    'error_ratio': one_qubit_ratio
                    }
    
            for g in one_qubit_gates:
                error_dict[g] = one_qubit_dict
            faulty_group5 = one_qubit_gates[:]
            
            if error_model[-1] == '2':
                error_dict['IMS5'] = {'error_rate': 2*p_2q, 'error_ratio': fiveq_ratio} 
                faulty_group2 = ['IMS5', 'MS']
                faulty_groups = [faulty_group1, faulty_group2, faulty_group3,
                                 faulty_group4, faulty_group5]            
            else:
                error_dict['IMS5'] = {'error_rate': 2*p_5q, 'error_ratio': fiveq_ratio} 
                faulty_group2 = ['MS']
                faulty_group6 = ['IMS5']
                faulty_groups = [faulty_group1, faulty_group2, faulty_group3,
                                 faulty_group4, faulty_group5, faulty_group6]            


        else:
            # First eQual error model

            # Errors after shuttling, merging and cooling
            error_dict['Ism'] = {'error_rate': p_sm, 'error_ratio': {'Z': 1}}
            faulty_group3 = ['Ism']
            error_dict['Icool'] = {'error_rate': p_cool, 'error_ratio': {'Z': 1}}
            faulty_group4 = ['Icool']
            faulty_groups = [faulty_group1, faulty_group2, faulty_group3, faulty_group4]

        return error_dict, Is_after_two_qubit, Is_after_one_qubit, faulty_groups


    else:
        raise NameError('This error model is not implemented.')

    return error_dict, Is_after_two_qubit, Is_after_one_qubit



def add_errors_fast_sampler(gate_indices, n_errors, subcircs, error_info):
    '''
    '''
    sampling = 'Muyalon'
    n_subcircs = len(subcircs)   # For now: 3 (Cross and 5-qubit)
                                 #          6 (Steane)

    # get list of indices for one qubit gates and two qubit gates
    one_q_gates = gate_indices[0]
    two_q_gates = gate_indices[1]
    
    # shuffle the gate indices
    rd.shuffle(one_q_gates)
    rd.shuffle(two_q_gates)

    selected_one_q_gates = one_q_gates[ : n_errors[0]]
    # print "selected one qubit gates are", selected_one_q_gates
    selected_two_q_gates = two_q_gates[ : n_errors[1]]
    # print "selected two qubit gates are", selected_two_q_gates
 
    carry_run = False   
    faulty_subcircs = []
    errors_dict = {}
    for i in range(n_subcircs):
        subcirc = subcircs[i]
        # local gates are the gates that we are adding errors after
        local_gates = []

        for pair in selected_one_q_gates:
            if i == pair[0]: 
                # the selected gate must be in the current redundency subcircuit
                local_gates += [pair[1]]

        for pair in selected_two_q_gates:
            if i == pair[0]:
                local_gates += [pair[1]]

        if len(local_gates) != 0:           
            # if no gates selected at all, then no need to perform the add error step at all
	        error.add_error_alternative(subcirc, error_info, sampling, local_gates)

        faulty_subcircs += [subcirc]
        
        subcirc_dict, local_carry_run = get_errors_dict(subcirc)
        if i >= int(2*n_subcircs/3): 
            # 2 for Cross and 5-qubit; 4 for Steane
            local_carry_run = False
        carry_run = carry_run or local_carry_run
        errors_dict[i] = subcirc_dict

    return errors_dict, carry_run



def add_errors_fast_sampler_new(gate_indices, n_errors, subcircs, error_info,
                                sub_circ_break=1):
    '''
    '''
    sampling = 'Muyalon'
    n_subcircs = len(subcircs)   # For now: 3 (Cross and 5-qubit)
                                 #          6 (Steane)

    # get list of indices for one qubit gates and two qubit gates
    one_q_gates = gate_indices[0]
    two_q_gates = gate_indices[1]

    # shuffle the gate indices
    rd.shuffle(one_q_gates)
    rd.shuffle(two_q_gates)

    selected_one_q_gates = one_q_gates[ : n_errors[0]]
    #print "selected one qubit gates are", selected_one_q_gates
    selected_two_q_gates = two_q_gates[ : n_errors[1]]
    #print "selected two qubit gates are", selected_two_q_gates
    
    sorted_one_q_gates = sorted(selected_one_q_gates, key=lambda gate: gate[0])
    sorted_two_q_gates = sorted(selected_two_q_gates, key=lambda gate: gate[0])

    if n_errors[0] > 0:
        subcirc_one_q_gates_indices = [gate[0] for gate in sorted_one_q_gates]      
    else:
        subcirc_one_q_gates_indices = [2]

    if n_errors[1] > 0:
        subcirc_two_q_gates_indices = [gate[0] for gate in sorted_two_q_gates]
    else:
        subcirc_two_q_gates_indices = [2]
    
    total_indices = subcirc_one_q_gates_indices + subcirc_two_q_gates_indices
    sorted_total_indices = sorted(total_indices)

    # if all the errors are added after the first round of QEC, then we don't
    # even run the circuit.
    # MGA 12/23/19.  Commented out this for the NN decoder.
    #if sorted_total_indices[0] > sub_circ_break:
    #    return [selected_one_q_gates, selected_two_q_gates], False, None


    carry_run = True
    faulty_subcircs = []
    errors_dict = {}
    errors_added_total = []
    for i in range(n_subcircs):
        subcirc = copy.deepcopy(subcircs[i])
        #subcirc = subcircs[i]
        # local gates are the gates that we are adding errors after
        local_gates = []

        for pair in selected_one_q_gates:
            if i == pair[0]: 
                # the selected gate must be in the current redundency subcircuit
                local_gates += [pair[1]]

        for pair in selected_two_q_gates:
            if i == pair[0]:
                local_gates += [pair[1]]

        if len(local_gates) != 0:           
            # if no gates selected at all, then no need to perform the add error step at all
            errors_added_local = error.add_error_alternative(subcirc, error_info, 
                                                             sampling, local_gates)

            #print i
            #print 'added errors locally =', errors_added_local
            errors_added_total += [[i,errors_added_local]]


        faulty_subcircs += [subcirc]
        
        #subcirc_dict, local_carry_run = get_errors_dict(subcirc)
        #if i >= int(2*n_subcircs/fraction_of_circ): 
        #    local_carry_run = False
        #carry_run = carry_run or local_carry_run
        #errors_dict[i] = subcirc_dict
        
    #print errors_added_total

    #return errors_dict, carry_run, faulty_subcircs
    return [selected_one_q_gates, selected_two_q_gates], carry_run, faulty_subcircs, errors_added_total



def add_errors_general_sampler_new(subcircs, error_info, sub_circ_break=1):
    '''
    '''
    sampling = 'old'
    n_subcircs = len(subcircs)   # For now: 3 (Cross and 5-qubit)
                                 #          6 (Steane)


    # if all the errors are added after the first round of QEC, then we don't
    # even run the circuit.
    # MGA 12/23/19.  Commented out this for the NN decoder.
    #if sorted_total_indices[0] > sub_circ_break:
    #    return [selected_one_q_gates, selected_two_q_gates], False, None

    carry_run = True
    faulty_subcircs = []
    errors_dict = {}
    errors_added_total = []
    for i in range(n_subcircs):
        subcirc = copy.deepcopy(subcircs[i])
        error.add_error_alternative(subcirc, error_info, sampling)
        faulty_subcircs += [subcirc]
        
    
    return carry_run, faulty_subcircs



def add_errors_fast_sampler_ion(gate_indices, n_errors, subcircs, error_info,
                                n_subcircs_first_round=12):
    '''
    gate_indices:  a list of lists, each one corresponding to the indices for 
                   each particular kind of gate.

    Right now the output variable carry_run is always True.

    '''
    sampling = 'Muyalon'
    n_subcircs = len(subcircs)

    # shuffle each list of gate indices and select the faulty gates
    total_selected_gates = []
    subcirc_gates_indices = []
    total_indices = []
    for i in range(len(n_errors)):
        rd.shuffle(gate_indices[i])  # shuffle one of the lists
        selected_gates = gate_indices[i][:n_errors[i]]  # select faulty gates
        sorted_selected_gates = sorted(selected_gates, key=lambda gate: gate[0])
        total_selected_gates += [sorted_selected_gates]

        # define a list with the subcirc indices for each faulty gate
        # if this particular faulty gate has no errors, then we add
        # 1 to the number of subcircs on the first round of QEC.
        # This is done just to make it easier to determine if the circuit
        # will be run at all.
        if n_errors[i] > 0:
            subcirc_gates_indices += [[gate[0] for gate in sorted_selected_gates]]
        else:
            subcirc_gates_indices += [[n_subcircs_first_round + 1]]
        
        # total_indices is just the flat-list version of subcirc_gates_indices
        total_indices += subcirc_gates_indices[-1]

    sorted_total_indices = sorted(total_indices)


    # Decide whether or not to run the circuit
    # First step:  determine which one is the first even index.
    # Even indices correspond to FT circuits, whereas odd indices
    # correspond to non-FT circuits.  If the first even index is
    # after the first round of QEC, then we don't run the circuit
    # because it means that no error would occur on the first round.
    

    # Commenting out these lines to run the new decoder
    #first_even_index = len(subcircs)
    #for index in sorted_total_indices:
    #    if index%2 == 0:
    #        first_even_index = index
    #        break

    #print first_even_index

    #if first_even_index >= n_subcircs_first_round:
    #    return total_selected_gates, False, None

    
    carry_run = True
    faulty_subcircs = []
    for i in range(n_subcircs):
        subcirc = copy.deepcopy(subcircs[i])
        
        # local gates are the gates that we are adding errors after
        local_gates = []

        # figure out if this particular subcircuit has errors
        for i_gate in range(len(n_errors)):
            for pair in total_selected_gates[i_gate]:
                if i == pair[0]:
                    # the selected gate must be in the current subcircuit
                    local_gates += [pair[1]]

        if len(local_gates) != 0:           
            # if no gates selected at all, then no need to perform the add error step at all
	        error.add_error_alternative(subcirc, error_info, sampling, local_gates)

        faulty_subcircs += [subcirc]
        
    return total_selected_gates, carry_run, faulty_subcircs



def add_errors_fast_sampler_ion_latt_surg(gate_indices, n_errors, latt_circs, error_info,
                                          n_subcircs_first_round=12):
    '''
    gate_indices:  a list of lists, each one corresponding to the indices for 
                   each particular kind of gate.

    Right now the output variable carry_run is always True.

    '''
    sampling = 'Muyalon'
    n_subcircs = len(latt_circs)

    # shuffle each list of gate indices and select the faulty gates
    total_selected_gates = []
    subcirc_gates_indices = []
    total_indices = []
    for i in range(len(n_errors)):
        rd.shuffle(gate_indices[i])  # shuffle one of the lists
        selected_gates = gate_indices[i][:n_errors[i]]  # select faulty gates
        sorted_selected_gates = sorted(selected_gates, key=lambda gate: gate[0])
        #total_selected_gates += [sorted_selected_gates]
        total_selected_gates += sorted_selected_gates
        
        # define a list with the subcirc indices for each faulty gate
        # if this particular faulty gate has no errors, then we add
        # 1 to the number of subcircs on the first round of QEC.
        # This is done just to make it easier to determine if the circuit
        # will be run at all.
        #if n_errors[i] > 0:
        #    subcirc_gates_indices += [[gate[0] for gate in sorted_selected_gates]]
        #else:
        #    subcirc_gates_indices += [[n_subcircs_first_round + 1]]
        
        # total_indices is just the flat-list version of subcirc_gates_indices
        #total_indices += subcirc_gates_indices[-1]

    #sorted_total_indices = sorted(total_indices)

    # group the selected gates
    gate_groups = []
    for gate in total_selected_gates:
        in_group = False
        for group in gate_groups:
            for g in group:
                if g[:-1] == gate[:-1]:
                    group.insert(0, gate)
                    in_group = True
                    break
        if not in_group:
            gate_groups += [[gate]]

    #print 'total selected gates =', total_selected_gates
    #print 'gate groups =', gate_groups

    # insert errors
    for group in gate_groups:
        local_gates = [g[-1] for g in group]
        if len(group[0]) >= 2:
            faulty_circ = latt_circs.gates[group[0][0]].circuit_list[0]
        if len(group[0]) >= 3:
            faulty_circ = faulty_circ.gates[group[0][1]].circuit_list[0]
        if len(group[0]) == 4:
            faulty_circ = faulty_circ.gates[group[0][2]].circuit_list[0]

        error.add_error_alternative(faulty_circ, error_info, 'Muyalon', local_gates)

    return



def add_errors_fast_sampler_color(gate_indices, n_errors, subcircs, error_info):
    '''
    '''
    sampling = 'Muyalon'
    n_subcircs = len(subcircs) 

    # get list of indices for one qubit gates and two qubit gates
    one_q_gates = gate_indices[0]
    two_q_gates = gate_indices[1]

    # shuffle the gate indices
    rd.shuffle(one_q_gates)
    rd.shuffle(two_q_gates)

    selected_one_q_gates = one_q_gates[ : n_errors[0]]
    sorted_one_q_gates = sorted(selected_one_q_gates, key=lambda gate: gate[0])
    selected_two_q_gates = two_q_gates[ : n_errors[1]]
    sorted_two_q_gates = sorted(selected_two_q_gates, key=lambda gate: gate[0])


    if n_errors[0] > 0:
        subcirc_one_q_gates_indices = [gate[0] for gate in sorted_one_q_gates]      
    else:
        subcirc_one_q_gates_indices = [33]

    if n_errors[1] > 0:
        subcirc_two_q_gates_indices = [gate[0] for gate in sorted_two_q_gates]
    else:
        subcirc_two_q_gates_indices = [33]

    #print sorted_one_q_gates, sorted_two_q_gates
    #print subcirc_one_q_gates_indices, subcirc_two_q_gates_indices
    
    total_indices = subcirc_one_q_gates_indices + subcirc_two_q_gates_indices
    sorted_total_indices = sorted(total_indices)
    
    #print sorted_total_indices

    first_even_index = 128
    for index in sorted_total_indices:
        if index%2 == 0:
            first_even_index = index
            break

    #print first_even_index

    if first_even_index > 31:
        return [selected_one_q_gates, selected_two_q_gates], False, None


    #if subcirc_one_q_gates_indices[0] > 31 and subcirc_two_q_gates_indices[0] > 31:
    #    return {}, False, None

    #all_subcircs_odd = True
    #for subcirc in subcirc_one_q_gates_indices + subcirc_two_q_gates_indices:
    #    if subcirc%2 == 0:
    #        all_subcircs_odd = False
    #        break

    #if all_subcircs_odd:
    #    return {}, False, None
    

    carry_run = True
    faulty_subcircs = []
    errors_dict = {}
    for i in range(n_subcircs):
        subcirc = copy.deepcopy(subcircs[i])
        #subcirc = subcircs[i]
        # local gates are the gates that we are adding errors after
        local_gates = []

        for pair in selected_one_q_gates:
            if i == pair[0]: 
                # the selected gate must be in the current redundency subcircuit
                local_gates += [pair[1]]

        for pair in selected_two_q_gates:
            if i == pair[0]:
                local_gates += [pair[1]]

        if len(local_gates) != 0:           
            # if no gates selected at all, then no need to perform the add error step at all
	        error.add_error_alternative(subcirc, error_info, sampling, local_gates)

        faulty_subcircs += [subcirc]
        
        #subcirc_dict, local_carry_run = get_errors_dict(subcirc)
        #if i >= int(2*n_subcircs/fraction_of_circ): 
        #    local_carry_run = False
        #carry_run = carry_run or local_carry_run
        #errors_dict[i] = subcirc_dict

    #return errors_dict, carry_run, faulty_subcircs
    return [selected_one_q_gates, selected_two_q_gates], carry_run, faulty_subcircs



def add_specific_error_configuration(circ_list, errors_to_add, gate_indexes):
    '''
    '''
    if len(errors_to_add) != len(gate_indexes):
        raise ValueError('The errors and the gates cannot have different length')

    #print 'Indexes =', gate_indexes
    #print 'Errors =', errors_to_add
    #print gate_indexes[0][0]

    faulty_circ_list = []
    for subcirc_i in range(len(circ_list)):
        faulty_circ = circ_list[subcirc_i]
        local_gates, local_errors = [], []
        for gate_i in range(len(gate_indexes)):
            if subcirc_i == gate_indexes[gate_i][0]:
                local_gates += [gate_indexes[gate_i][1]]
                local_errors += [errors_to_add[gate_i]]
        
        if len(local_gates) != 0:
            qfun.add_errors_after_gates(faulty_circ, local_gates, local_errors)

        faulty_circ_list += [faulty_circ]

    return faulty_circ_list



def add_errors_fast_sampler_surface(gate_indices, n_errors, circ, error_info):
    '''
    '''
    sampling = 'Muyalon'

    # get list of indices for one qubit gates and two qubit gates
    one_q_gates = gate_indices[0]
    two_q_gates = gate_indices[1]
    
    # shuffle the gate indices
    rd.shuffle(one_q_gates)
    rd.shuffle(two_q_gates)

    selected_one_q_gates = one_q_gates[ : n_errors[0]]
    print selected_one_q_gates
    # print "selected one qubit gates are", selected_one_q_gates
    selected_two_q_gates = two_q_gates[ : n_errors[1]]
    print selected_two_q_gates
    # print "selected two qubit gates are", selected_two_q_gates
 
    carry_run = False   
    faulty_subcircs = []
    errors_dict = {}
    # local gates are the gates that we are adding errors after
    local_gates = []

    for gate in selected_one_q_gates:
        # the selected gate must be in the current redundency subcircuit
        local_gates.append(gate)

    for gate in selected_two_q_gates:
        local_gates.append(gate)

    if len(local_gates) != 0:           
        # if no gates selected at all, then no need to perform the add error step at all
        error.add_error_alternative(circ, error_info, sampling, local_gates)

    # faulty_subcircs += [subcirc]
    
    circ_dict, local_carry_run = get_errors_dict(circ)

    carry_run = carry_run or local_carry_run
    errors_dict = circ_dict

    return errors_dict, carry_run



def gates_list_CNOT(CNOT_circuits, faulty_gates_names):
    '''
    improvised function to calculate the indices for 1-qubit and 2-qubit gates
    '''

    single_qubit_gates, two_qubit_gates = [], []

    for i in range(len(CNOT_circuits.gates)):
        supra_gate = CNOT_circuits.gates[i] 
        if supra_gate.gate_name == 'Logical_I' or supra_gate.gate_name == 'MeasureX':
            for j in range(len(supra_gate.circuit_list[0].gates)):
                in_gate1 = supra_gate.circuit_list[0].gates[j]
                if in_gate1.gate_name in faulty_gates_names:
                    if len(in_gate1.qubits) == 1:
                        single_qubit_gates.append((i,j))
                    elif len(in_gate1.qubits) == 2:
                        two_qubit_gates.append((i,j))
        
        elif supra_gate.gate_name[:8] == 'Measure2' or supra_gate.gate_name[:5] == 'Joint':
            for j in range(len(supra_gate.circuit_list[0].gates)):
                in_gate1 = supra_gate.circuit_list[0].gates[j]
                if in_gate1.gate_name[:7] == 'Partial':
                    for k in range(len(in_gate1.circuit_list[0].gates)):
                        in_gate2 = in_gate1.circuit_list[0].gates[k]
                        if in_gate2.gate_name in faulty_gates_names:
                            if len(in_gate2.qubits) == 1:
                                single_qubit_gates.append((i,j,k))
                            elif len(in_gate2.qubits) == 2:
                                two_qubit_gates.append((i,j,k))
                            
                elif in_gate1.gate_name[:2] == 'EC':
                    for k in range(len(in_gate1.circuit_list[0].gates)):
                        in_gate2 = in_gate1.circuit_list[0].gates[k]
                        for l in range(len(in_gate2.circuit_list[0].gates)):
                            in_gate3 = in_gate2.circuit_list[0].gates[l]
                            if in_gate3.gate_name in faulty_gates_names:
                                if len(in_gate3.qubits) == 1:
                                    single_qubit_gates.append((i,j,k,l))
                                elif len(in_gate3.qubits) == 2:
                                    two_qubit_gates.append((i,j,k,l))

        elif supra_gate.gate_name[:9]=='MeasureXX' or supra_gate.gate_name[:9]=='MeasureZZ':
            for j in range(len(supra_gate.circuit_list[0].gates)):
                in_gate1 = supra_gate.circuit_list[0].gates[j]
                clauseX = in_gate1.gate_name[0] == 'X'
                clauseZ = in_gate1.gate_name[0] == 'Z'
                clauseC = in_gate1.gate_name[0] == 'C'
                if (clauseX or clauseZ) or clauseC:
                    for k in range(len(in_gate1.circuit_list[0].gates)):
                        in_gate2 = in_gate1.circuit_list[0].gates[k]
                        if in_gate2.gate_name in faulty_gates_names:
                            if len(in_gate2.qubits) == 1:
                                single_qubit_gates.append((i,j,k))
                            elif len(in_gate2.qubits) == 2:
                                two_qubit_gates.append((i,j,k))
                
                elif in_gate1.gate_name[:3] == 'QEC':
                    for k in range(len(in_gate1.circuit_list[0].gates)):
                        in_gate2 = in_gate1.circuit_list[0].gates[k]
                        for l in range(len(in_gate2.circuit_list[0].gates)):
                            in_gate3 = in_gate2.circuit_list[0].gates[l]
                            if in_gate3.gate_name in faulty_gates_names:
                                if len(in_gate3.qubits) == 1:
                                    single_qubit_gates.append((i,j,k,l))
                                elif len(in_gate3.qubits) == 2:
                                    two_qubit_gates.append((i,j,k,l))
                                
    return single_qubit_gates, two_qubit_gates
    


def gates_list_CNOT_general(CNOT_circuits, faulty_gates_names_grouped):
    '''
    improvised function to calculate the indices for 1-qubit and 2-qubit gates
    '''
    
    out_gates_lists = [[] for group in faulty_gates_names_grouped] 

    for i in range(len(CNOT_circuits.gates)):
        supra_gate = CNOT_circuits.gates[i] 
        if supra_gate.gate_name == 'Logical_I' or supra_gate.gate_name == 'MeasureX':
            for j in range(len(supra_gate.circuit_list[0].gates)):
                in_gate1 = supra_gate.circuit_list[0].gates[j]
                for fgg in range(len(out_gates_lists)):
                    if in_gate1.gate_name in faulty_gates_names_grouped[fgg]:
                        out_gates_lists[fgg].append((i,j))
                 
        elif supra_gate.gate_name[:8] == 'Measure2' or supra_gate.gate_name[:5] == 'Joint':
            for j in range(len(supra_gate.circuit_list[0].gates)):
                in_gate1 = supra_gate.circuit_list[0].gates[j]
                if in_gate1.gate_name[:7] == 'Partial':
                    for k in range(len(in_gate1.circuit_list[0].gates)):
                        in_gate2 = in_gate1.circuit_list[0].gates[k]
                        for fgg in range(len(out_gates_lists)):
                            if in_gate2.gate_name in faulty_gates_names_grouped[fgg]:
                                out_gates_lists[fgg].append((i,j,k))
                      
                elif in_gate1.gate_name[:2] == 'EC':
                    for k in range(len(in_gate1.circuit_list[0].gates)):
                        in_gate2 = in_gate1.circuit_list[0].gates[k]
                        for l in range(len(in_gate2.circuit_list[0].gates)):
                            in_gate3 = in_gate2.circuit_list[0].gates[l]
                            for fgg in range(len(out_gates_lists)):
                                if in_gate3.gate_name in faulty_gates_names_grouped[fgg]:
                                    out_gates_lists[fgg].append((i,j,k,l))

        elif supra_gate.gate_name[:9]=='MeasureXX' or supra_gate.gate_name[:9]=='MeasureZZ':
            for j in range(len(supra_gate.circuit_list[0].gates)):
                in_gate1 = supra_gate.circuit_list[0].gates[j]
                clauseX = in_gate1.gate_name[0] == 'X'
                clauseZ = in_gate1.gate_name[0] == 'Z'
                clauseC = in_gate1.gate_name[0] == 'C'
                if (clauseX or clauseZ) or clauseC:
                    for k in range(len(in_gate1.circuit_list[0].gates)):
                        in_gate2 = in_gate1.circuit_list[0].gates[k]
                        for fgg in range(len(out_gates_lists)):
                            if in_gate2.gate_name in faulty_gates_names_grouped[fgg]:
                                out_gates_lists[fgg].append((i,j,k))
                        
                elif in_gate1.gate_name[:3] == 'QEC':
                    for k in range(len(in_gate1.circuit_list[0].gates)):
                        in_gate2 = in_gate1.circuit_list[0].gates[k]
                        for l in range(len(in_gate2.circuit_list[0].gates)):
                            in_gate3 = in_gate2.circuit_list[0].gates[l]
                            for fgg in range(len(out_gates_lists)):
                                if in_gate3.gate_name in faulty_gates_names_grouped[fgg]:
                                    out_gates_lists[fgg].append((i,j,k,l))
                                
    return out_gates_lists
    
