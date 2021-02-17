import sys
import os
import time
import json
import copy
import random as rd
import multiprocessing as mp
import circuit as cir
import steane
import surface17 as surf17
import correction as cor
import chper_wrapper as wrapper
import MC_functions as mc
import qcircuit_functions as qfun
import qcircuit_wrapper as qwrap
from visualizer import browser_vis as brow
import itertools


chp_loc = './chp_extended'
error_model = 'ion_trap0_NN'
p1, p2, p_meas = 0.001, 0.001, 0.001  # these don't matter for the fast sampler
error_dict, Is_after2q, Is_after_1q = wrapper.dict_for_error_model(error_model, p1, p2, p_meas)
error_info = mc.read_error_info(error_dict)




n_per_proc, n_proc = int(sys.argv[1]), int(sys.argv[2])
n_errors = [int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])]
QEC_kind = sys.argv[6]   # surface17
state = sys.argv[7]  # either 'Z' or 'X'
n_rounds = int(sys.argv[8])
even_spaced = True
if n_rounds%3 != 0:  only_NN = True
else:  only_NN = False


# Define initial state
if QEC_kind == 'surface17':
    init_stabs = surf17.Code.stabilizers_CHP[state][:]
    init_destabs = surf17.Code.destabilizers_CHP[state][:]
    init_state = [init_stabs, init_destabs]


if QEC_kind == 'surface17':
    output_folder = './MC_results/QECd3_surface17/'+error_model+'/'+str(n_rounds)+'/'+state+'/'




def generate_wait_circ(n_data):
    '''
    '''

    circ_I = cir.Circuit()
    for i in range(n_data):
        circ_I.add_gate_at([i], 'I_idle')

    return circ_I



def generate_bare_meas_ion(n_data, stabilizer_list):
    '''
    Generates stabilizer measurements with bare ancilla, no cat states.
    parallel = True means that we don't recicle the ancilla

    stabilizer_list is assumed to have the following format:
        [
          [('X',0), ('X',4)], ... 

        ]
    '''

    n_ancilla = len(stabilizer_list)
    
    bare_meas_circ = cir.Circuit()
    for i in range(n_data):
        bare_meas_circ.add_gate_at([i], 'I_gate')

    for i in range(n_data, n_data+n_ancilla):
        bare_meas_circ.add_gate_at([i], 'PrepareXPlus')
        
    for i in range(len(stabilizer_list)):
        gateprefix = 'C'
        for gate in stabilizer_list[i]:
            bare_meas_circ.add_gate_at([n_data+i,gate[1]], gateprefix+gate[0])
            for j in range(n_data+n_ancilla):
                bare_meas_circ.add_gate_at([j], 'I_gate')

    for i in range(len(stabilizer_list)):
        bare_meas_circ.add_gate_at([n_data+i], 'ImX')
        bare_meas_circ.add_gate_at([n_data+i], 'MeasureX')
       
    bare_meas_circ.to_ancilla(range(n_data, n_data+n_ancilla))

    return bare_meas_circ



def generate_rep_bare_meas_ion(n_data, stabilizer_list, n_rounds, even_spaced=True):
    '''
    '''

    n_stab = len(stabilizer_list)/2
    
    rep_meas_circ = cir.Circuit()
    circ_I = generate_wait_circ(n_data)
    circ_I = cir.Encoded_Gate('WAIT', [circ_I]).circuit_wrap()
    rep_meas_circ.join_circuit(circ_I, False)


    for i in range(n_rounds):
        gate_nameX = 'SX_%i'%i
        stab_circX = generate_bare_meas_ion(n_data, stabilizer_list[:n_stab])
        stab_circX = cir.Encoded_Gate(gate_nameX, [stab_circX]).circuit_wrap()
        gate_nameZ = 'SZ_%i'%i
        stab_circZ = generate_bare_meas_ion(n_data, stabilizer_list[n_stab:])
        stab_circZ = cir.Encoded_Gate(gate_nameZ, [stab_circZ]).circuit_wrap()
        rep_meas_circ.join_circuit(stab_circX, False)
        rep_meas_circ.join_circuit(stab_circZ, False)

        if even_spaced:
            circ_I = generate_wait_circ(n_data)
            circ_I = cir.Encoded_Gate('WAIT', [circ_I]).circuit_wrap()
            rep_meas_circ.join_circuit(circ_I, False)

    if not even_spaced:
        circ_I = generate_wait_circ(n_data)
        circ_I = cir.Encoded_Gate('WAIT', [circ_I]).circuit_wrap()
        rep_meas_circ.join_circuit(circ_I, False)

    return rep_meas_circ
    



# Define the circuit and the circuit_list
if QEC_kind == 'surface17':
    surface17_stabs = surf17.Code.stabilizers[:]
    QEC_circ = generate_rep_bare_meas_ion(9, surface17_stabs, n_rounds, even_spaced)
    QEC_circ_list = []
    for supra_gate in QEC_circ.gates:
        QEC_circ_list += [supra_gate.circuit_list[0]]


#brow.from_circuit(QEC_circ, True)

# Define the list of error-prone 1-q and 2-q gates
faulty_gates_grouped = [['PrepareXPlus','ImX','CX','CZ'], ['I_gate'], ['I_idle']]
regular_gates, I_gates, idle_gates = wrapper.gates_list_general(QEC_circ_list, faulty_gates_grouped)

#print len(regular_gates)
#print len(I_gates)
#print len(idle_gates)

#sys.exit(0)




# Define the number of all the 2-q gates, 1-q gates, and measurements
# for resource counting purposes
n_2q_gates, n_1q_gates, n_meas = [], [], []
for subcirc in QEC_circ_list:
    two_q, one_q, meas = 0, 0, 0
    for phys_gate in subcirc.gates:
        if len(phys_gate.qubits) > 1:
            two_q += 1
        else:
            if phys_gate.gate_name[0] != 'I':
                one_q += 1
            if phys_gate.gate_name[:4] == 'Meas':
                meas += 1
    n_2q_gates += [two_q]
    n_1q_gates += [one_q]
    n_meas += [meas]






def run_QEC_d3(init_state, QEC_circ_list, kind='flag'):
    '''
    kind: 'flag' or 'diVin'
    '''
    if kind=='surface17':
        QEC_object = qwrap.QEC_d3(init_state, QEC_circ_list[:], chp_loc)
        # MGA 12/23/19: we will use the old metadecoder.
        sim_output = QEC_object.run_fullQEC_CSS_d3_rounds_even('surface17', 3)
        n_X, n_Z, X_stab_outcomes, Z_stab_outcomes, Z_corr, X_corr, unc_stabs, unc_destabs = sim_output
        n_supra_gates = n_X + n_Z

    
    final_stabs, final_destabs = QEC_object.stabs[:], QEC_object.destabs[:]
    
    # Determine if there is an error (both failures and correctable errors)
    final_error = False
    for stab in final_stabs:
        if stab[0] != '+':
            final_error = True
            break

    # do perfect EC
    # for both the corrected state (original) and the uncorrected states (for the NN decoder)
    # MGA 1/29/2020
    if kind=='surface17':
        surface17_stabs = surf17.Code.stabilizers[:]
        corr_circ = cor.Bare_Correct.generate_rep_bare_meas(9, surface17_stabs, 2, False, True,
                                                            False, False, False, True)
        corr_circ_list = []
        for supra_gate in corr_circ.gates:
            corr_circ_list += [supra_gate.circuit_list[0]]

        corr_object = qwrap.QEC_d3([final_stabs[:],final_destabs[:]], corr_circ_list[:], chp_loc) 
        # don't matter 1 and don't matter 2. The others also don't matter.
        dm1, dm2, X_stab, Z_stab, blahZ, blahX, blah_stab, blah_destab = corr_object.run_fullQEC_CSS_d3('surface17', True, False)
        corr_stabs = corr_object.stabs[:]

        # New lines added for the NN decoder.  MGA 1/29/2020:
        corr_object_NN = qwrap.QEC_d3([unc_stabs[:],unc_destabs[:]], corr_circ_list[:], chp_loc) 
        # don't matter 1 and don't matter 2. The others also don't matter.
        dm1, dm2, X_stab, Z_stab, blahZ, blahX, blah_stab, blah_destab = corr_object_NN.run_fullQEC_CSS_d3('surface17', True, False)
        corr_stabs_NN = corr_object_NN.stabs[:]
        


    # Determine if a failure has occured (for the lookup table decoder).
    fail = False
    for stab in corr_stabs:
        if stab[0] != '+':
            fail = True
            break

    # Determine if a failure has occured (for the NN decoder).
    fail_NN = False
    for stab in corr_stabs_NN:
        if stab[0] != '+':
            fail_NN = True
            break
    
    
    return final_error, fail, fail_NN, n_supra_gates, X_stab_outcomes, Z_stab_outcomes, Z_corr, X_corr



def run_QEC_d3_onlyNN(init_state, QEC_circ_list, kind='flag'):
    '''
    kind: 'flag' or 'diVin'
    '''
    if kind=='surface17':
        QEC_object = qwrap.QEC_d3(init_state, QEC_circ_list[:], chp_loc)
        # MGA 12/23/19: we will use the old metadecoder.
        sim_output = QEC_object.run_fullQEC_CSS_d3_rounds_even_onlyNN('surface17', n_rounds)
        n_X, n_Z, X_stab_outcomes, Z_stab_outcomes = sim_output
        n_supra_gates = n_X + n_Z

    
    final_stabs, final_destabs = QEC_object.stabs[:], QEC_object.destabs[:]
    
    # Determine if there is an error (both failures and correctable errors)
    final_error = False
    for stab in final_stabs:
        if stab[0] != '+':
            final_error = True
            break

    # do perfect EC
    # Only for the uncorrected states (for the NN decoder)
    # MGA 6/9/2020
    if kind=='surface17':
        surface17_stabs = surf17.Code.stabilizers[:]
        corr_circ = cor.Bare_Correct.generate_rep_bare_meas(9, surface17_stabs, 2, False, True,
                                                            False, False, False, True)
        corr_circ_list = []
        for supra_gate in corr_circ.gates:
            corr_circ_list += [supra_gate.circuit_list[0]]

        corr_object_NN = qwrap.QEC_d3([final_stabs[:],final_destabs[:]], corr_circ_list[:], chp_loc) 
        # don't matter 1 and don't matter 2. The others also don't matter.
        dm1, dm2, X_stab, Z_stab, blahZ, blahX, blah_stab, blah_destab = corr_object_NN.run_fullQEC_CSS_d3('surface17', True, False)
        corr_stabs_NN = corr_object_NN.stabs[:]
        

    # Determine if a failure has occured (for the NN decoder).
    fail_NN = False
    for stab in corr_stabs_NN:
        if stab[0] != '+':
            fail_NN = True
            break
    
    
    return final_error, fail_NN, n_supra_gates, X_stab_outcomes, Z_stab_outcomes



def run_several_QEC_fast(error_info, n_runs_total, init_state, QEC_kind, QEC_circ_list):
    '''
    '''

    did_run = 0
    n_final_errors = 0
    n_fails = 0
    n_fails_NN = 0
    n_supra_gates = 0
    errors_added = []
    X_stab_outcomes, Z_stab_outcomes = [], []
    Z_corrections, X_corrections = [], []
    final_errors, failings, failings_NN = [], [], []

    for n_run in xrange(n_runs_total):

        if QEC_kind == 'surface17':
            
            fraction_of_circ = 4
            
            QEC_circ_list_copy = []
            for subcirc in QEC_circ_list:
                QEC_circ_list_copy += [copy.deepcopy(subcirc)]

            # Add the errors and decide to run
            errors_dict, carry_run, faulty_circs = wrapper.add_errors_fast_sampler_ion(
                                                        [regular_gates, I_gates, idle_gates],
                                                        n_errors,
                                                        QEC_circ_list_copy,
                                                        error_info)

            #print errors_dict
            #new_errors_dict = [[], []] 
            #for oneq_error in errors_dict[0]:
            #    print oneq_error
            #    subcirc = faulty_circs[oneq_error[0]]
            #    gate_index = oneq_error[1]
            #    while subcirc.gates[gate_index].is_error == False:
            #        gate_index += 1
  
            #    gate = faulty_circs[oneq_error[0]].gates[gate_index]
            #    print gate.is_error
            #    print gate.gate_name 
            #    new_errors_dict[0] += [gate.gate_name, oneq_error]

            #for twoq_error in errors_dict[1]:
            #    print twoq_error
            #    subcirc = faulty_circs[twoq_error[0]]
            #    gate_index = twoq_error[1]
                 
            #    gate1 = faulty_circs[twoq_error[0]].gates[gate_index+1]
            #    gate2 = faulty_circs[twoq_error[0]].gates[gate_index+2]
 
            #    print gate1.is_error
            #    print gate1.gate_name 
            #    print gate2.is_error
            #    print gate2.gate_name 
                
            #    new_errors_dict[1] += [(gate1.gate_name, gate2.gate_name), twoq_error]
            


            #print new_errors_dict
            #brow.from_circuit(faulty_circs[0], True)
            #time.sleep(5)
            #print errors_dict
            
            # MGA 12/23/19: just for now in order to test NN decoder.
            carry_run = True
            #errors_added += [err_add_circ]
        



        # Run
        if carry_run:
            did_run += 1
            final_error, fail, fail_NN, n_supra_local, X_out, Z_out, Z_correct, X_correct = run_QEC_d3(
                                                                                     init_state, 
                                                                                     faulty_circs, 
                                                                                     QEC_kind)
            X_stab_outcomes += [X_out]
            Z_stab_outcomes += [Z_out]
            Z_corrections += [Z_correct]
            X_corrections += [X_correct]
            n_supra_gates += n_supra_local
            final_errors += [final_error]
            failings += [fail]
            failings_NN += [fail_NN]
            if final_error:
                n_final_errors += 1
            if fail:  
                n_fails += 1
                #print errors_dict
                #for key in errors_dict:
                #    if errors_dict[key] != {}:
                #        brow.from_circuit(faulty_circs[key], True)
                #break
            if fail_NN:
                n_fails_NN += 1

        else:
            n_supra_gates += 2*len(QEC_circ_list)/fraction_of_circ


    return n_final_errors, n_fails, n_fails_NN, n_supra_gates, errors_added, X_stab_outcomes, Z_stab_outcomes, Z_corrections, X_corrections, final_errors, failings, failings_NN



def run_several_QEC_fast_onlyNN(error_info, n_runs_total, init_state, QEC_kind, QEC_circ_list):
    '''
    '''

    did_run = 0
    n_final_errors = 0
    n_fails = 0
    n_fails_NN = 0
    n_supra_gates = 0
    errors_added = []
    X_stab_outcomes, Z_stab_outcomes = [], []
    Z_corrections, X_corrections = [], []
    final_errors, failings, failings_NN = [], [], []

    for n_run in xrange(n_runs_total):

        if QEC_kind == 'surface17':
            
            fraction_of_circ = 4
            
            QEC_circ_list_copy = []
            for subcirc in QEC_circ_list:
                QEC_circ_list_copy += [copy.deepcopy(subcirc)]

            # Add the errors and decide to run
            errors_dict, carry_run, faulty_circs = wrapper.add_errors_fast_sampler_ion(
                                                        [regular_gates, I_gates, idle_gates],
                                                        n_errors,
                                                        QEC_circ_list_copy,
                                                        error_info)

            carry_run = True

        # Run
        if carry_run:
            did_run += 1
            final_error, fail_NN, n_supra_local, X_out, Z_out = run_QEC_d3_onlyNN(init_state, 
                                                                                  faulty_circs, 
                                                                                  QEC_kind)
            X_stab_outcomes += [X_out]
            Z_stab_outcomes += [Z_out]
            #Z_corrections += [Z_correct]
            #X_corrections += [X_correct]
            n_supra_gates += n_supra_local
            final_errors += [final_error]
            #failings += [fail]
            failings_NN += [fail_NN]
            if final_error:
                n_final_errors += 1
            #if fail:  
                #n_fails += 1
                #print errors_dict
                #for key in errors_dict:
                #    if errors_dict[key] != {}:
                #        brow.from_circuit(faulty_circs[key], True)
                #break
            if fail_NN:
                n_fails_NN += 1

        else:
            n_supra_gates += 2*len(QEC_circ_list)/fraction_of_circ


    return n_final_errors, n_fails_NN, n_supra_gates, errors_added, X_stab_outcomes, Z_stab_outcomes, final_errors, failings_NN



def run_parallel_QEC(error_info, n_runs_per_proc, n_proc, init_state, QEC_kind, QEC_circ_list, onlyNN):
    '''
    '''
    if QEC_kind == 'all_flags':
        sim_func = run_several_QEC_fast_all_flags
    elif not only_NN:
        sim_func = run_several_QEC_fast
    elif only_NN:
        sim_func = run_several_QEC_fast_onlyNN

    pool = mp.Pool()
    results = [pool.apply_async(sim_func, (error_info, n_runs_per_proc, 
                                           init_state, QEC_kind, QEC_circ_list[:]))
                    for proc in range(n_proc)]
    pool.close()
    pool.join()
    dicts = [r.get() for r in results]

    return dicts


#print run_several_QEC_fast(error_info, 10, init_state, QEC_kind, QEC_circ_list[:])
#print run_several_QEC_fast_all_flags(error_info, 10, init_state, QEC_kind)
#print n_per_proc, n_proc, QEC_kind
#sys.exit(0)

#print run_several_QEC_fast(error_info, 100, init_state, QEC_kind)
#out_list = run_parallel_QEC(error_info, n_per_proc, n_proc, init_state, 
#                            QEC_kind, QEC_circ_list)
n_total = n_per_proc*n_proc
#print out_list
#print n_total
#sys.exit(0)
#n_final_errors = sum([event[0] for event in out_list])
#n_fail = sum([event[1] for event in out_list])
#n_supra_gates = sum([event[2] for event in out_list])

##################### MGA 12/23/19:  for the NN decoder #########################
#output_list = run_several_QEC_fast(error_info, n_per_proc*n_proc, init_state, 
#                                   QEC_kind, QEC_circ_list)
#n_final_errors = output_list[0]
#n_fail = output_list[1]
#n_fail_NN = output_list[2]
#n_supra_gates = output_list[3]
#errors_added_total = output_list[4]
#X_stab_outcomes_total = output_list[5]
#Z_stab_outcomes_total = output_list[6]
#Z_corrections_total = output_list[7]
#X_corrections_total = output_list[8]
#final_errors_total = output_list[9]
#failings_total = output_list[10]
#failings_total_NN = output_list[11]
#################################################################################


if not only_NN:
    output_list = run_parallel_QEC(error_info, n_per_proc, n_proc, init_state,
       			       QEC_kind, QEC_circ_list, False)
    n_total = n_per_proc*n_proc
    n_final_errors = sum([event[0] for event in output_list])
    n_fail = sum([event[1] for event in output_list])
    n_fail_NN = sum([event[2] for event in output_list])
    n_supra_gates = sum([event[3] for event in output_list])
    errors_added_total = list(itertools.chain.from_iterable([event[4] for event in output_list]))
    X_stab_outcomes_total = list(itertools.chain.from_iterable([event[5] for event in output_list]))
    Z_stab_outcomes_total = list(itertools.chain.from_iterable([event[6] for event in output_list]))
    Z_corrections_total = list(itertools.chain.from_iterable([event[7] for event in output_list]))
    X_corrections_total = list(itertools.chain.from_iterable([event[8] for event in output_list]))
    final_errors_total = list(itertools.chain.from_iterable([event[9] for event in output_list]))
    failings_total = list(itertools.chain.from_iterable([event[10] for event in output_list]))
    failings_total_NN = list(itertools.chain.from_iterable([event[11] for event in output_list]))

    #n_even_gates = sum([event[3] for event in out_list])
    #n_odd_gates = sum([event[4] for event in out_list])

    #n_twoq_gates = n_even_gates*n_2q_gates[0] + n_odd_gates*n_2q_gates[1]
    #n_oneq_gates = n_even_gates*n_1q_gates[0] + n_odd_gates*n_1q_gates[1]
    #n_meas_gates = n_even_gates*n_meas[0] + n_odd_gates*n_meas[1]

    n_twoq_gates = n_supra_gates*n_2q_gates[0]
    n_oneq_gates = n_supra_gates*n_1q_gates[0]
    n_meas_gates = n_supra_gates*n_meas[0]

    n_correctable = n_final_errors - n_fail
    p_correctable = float(n_correctable)/float(n_total)
    p_fail = float(n_fail)/float(n_total)
    p_fail_NN = float(n_fail_NN)/float(n_total)
    p_supra_gates = float(n_supra_gates)/float(n_total)
    #p_even_supra = float(n_even_gates)/float(n_total)
    #p_odd_supra = float(n_odd_gates)/float(n_total)
    p_2q_gates = float(n_twoq_gates)/float(n_total)
    p_1q_gates = float(n_oneq_gates)/float(n_total)
    p_meas = float(n_meas_gates)/float(n_total)
    out_dict = {'n_total': n_total, 
                'n_correctable': n_correctable,
                'p_correctable': p_correctable,
                'n_fail': n_fail, 
                'p_fail': p_fail,
                'n_fail_NN': n_fail_NN,
                'p_fail_NN': p_fail_NN,
                'n_supra_gates': n_supra_gates,
                'p_supra_gates': p_supra_gates,
                #'p_even_supra': p_even_supra,
                #'p_odd_supra': p_odd_supra,
                'p_2q': p_2q_gates,
                'p_1q': p_1q_gates,
                'p_meas': p_meas,
                'error_added': errors_added_total,
                'X_stab_outcomes': X_stab_outcomes_total,
                'Z_stab_outcomes': Z_stab_outcomes_total,
                'Z_corrections': Z_corrections_total,
                'X_corrections': X_corrections_total,
                'final_errors': final_errors_total,
                'failings': failings_total,
                'failings_total_NN': failings_total_NN}
    #print out_dict
else:
    output_list = run_parallel_QEC(error_info, n_per_proc, n_proc, init_state,
       			       QEC_kind, QEC_circ_list, True)
    n_total = n_per_proc*n_proc
    n_final_errors = sum([event[0] for event in output_list])
    n_fail_NN = sum([event[1] for event in output_list])
    n_supra_gates = sum([event[2] for event in output_list])
    errors_added_total = list(itertools.chain.from_iterable([event[3] for event in output_list]))
    X_stab_outcomes_total = list(itertools.chain.from_iterable([event[4] for event in output_list]))
    Z_stab_outcomes_total = list(itertools.chain.from_iterable([event[5] for event in output_list]))
    final_errors_total = list(itertools.chain.from_iterable([event[6] for event in output_list]))
    failings_total_NN = list(itertools.chain.from_iterable([event[7] for event in output_list]))
    
    n_twoq_gates = n_supra_gates*n_2q_gates[0]
    n_oneq_gates = n_supra_gates*n_1q_gates[0]
    n_meas_gates = n_supra_gates*n_meas[0]

    p_fail_NN = float(n_fail_NN)/float(n_total)
    p_supra_gates = float(n_supra_gates)/float(n_total)
    #p_even_supra = float(n_even_gates)/float(n_total)
    #p_odd_supra = float(n_odd_gates)/float(n_total)
    p_2q_gates = float(n_twoq_gates)/float(n_total)
    p_1q_gates = float(n_oneq_gates)/float(n_total)
    p_meas = float(n_meas_gates)/float(n_total)
    out_dict = {'n_total': n_total, 
                'n_fail_NN': n_fail_NN,
                'p_fail_NN': p_fail_NN,
                'n_supra_gates': n_supra_gates,
                'p_supra_gates': p_supra_gates,
                'p_2q': p_2q_gates,
                'p_1q': p_1q_gates,
                'p_meas': p_meas,
                'error_added': errors_added_total,
                'X_stab_outcomes': X_stab_outcomes_total,
                'Z_stab_outcomes': Z_stab_outcomes_total,
                'final_errors': final_errors_total,
                'failings_total_NN': failings_total_NN}


if not os.path.exists(output_folder):
    os.makedirs(output_folder)

json_filename = str(n_errors[0]) + '_' + str(n_errors[1]) + '_' + str(n_errors[2]) + '.json'
abs_filename = output_folder + json_filename
json_file = open(abs_filename, 'w')
json.dump(out_dict, json_file, indent=4, separators=(',', ':'), sort_keys=True)
json_file.close()




