import sys
import copy
import circuit as cir
import steane as st
#import d5color
#import surface49 as surf49
import chper_extended as chper
import qcircuit_functions as qfun
from visualizer import browser_vis as brow



class Quantum_Operation(object):
    '''
    A quantum operation is a portion of a quantum circuit that
    consists of unitary operations and possibly measurements.
    Examples: - transversal logical gate
              - redundancy-3 QEC step for a distance-3 code
              - measurement of two logical operators for
                lattice surgery in the color code.
    '''

    def __init__(self, initial_state, circuits, chp_location, CHP_IO_files=False):
        '''
        Inputs:  - initial_state: generally without the ancilla;
                                  should be a list of lists:
                                  (1) stabs and (2) destabs.
                                  can also be an empty list in case the circuit
                                  already has prep gates.  This is the case for 
                                  the verification of the 8-qubit cat state.
                 - circuits: list of circuits that will be run
                             serially on chp.
        '''
            
        self.stabs = initial_state[0][:]
        self.destabs = initial_state[1][:]
        if len(self.stabs) == 0:
            self.n_d_q = len(circuits[0].data_qubits())
        else:
            self.n_d_q = len(self.stabs)  
            # number of total data qubits in the whole "supra circuit", not
            # necessarily in this particular circuit.

        self.circuits = circuits[:]
        self.chp_loc = chp_location
        self.CHP_IO_files = CHP_IO_files
        #print 'n_data =', self.n_d_q


    def run_one_circ(self, circuit): 
        '''
        - runs circuit on chp
        - updates the stabs and destabs
        - returns the dictionary of measurement outcomes
        Input: - circuit: either a circuit object or an index
        '''

        if type(circuit) == type(0):
            circuit = self.circuits[circuit]

        n_a_q = len(circuit.ancilla_qubits())
        #print 'n_data =', self.n_d_q
        #print 'n_anc =', n_a_q
        #print self.stabs

        circ_chp = chper.Chper(circ=circuit,
                               num_d_qub=self.n_d_q,
                               num_a_qub=n_a_q,
                               stabs=self.stabs[:],
                               destabs=self.destabs[:],
                               anc_stabs=[],
                               anc_destabs=[],
                               input_output_files=self.CHP_IO_files)

        circ_chp_output = circ_chp.run(self.chp_loc)
        dic = circ_chp_output[0]
        self.stabs = circ_chp_output[1][:]
        self.destabs = circ_chp_output[2][:]
        #print 'dict =', dic

        return dic
       
    
    
    def run_one_diVincenzo(self, circuit, code, stab_kind=None,
                           parity_oct=0):
        '''
        runs one round of X or Z stabilizers (for a CSS code)
        or the whole set of stabilizers (for a non-CSS code)
        assuming the ancillae are in an unverified cat state.
        circuit refers to the index of the circuit to be run.
        Based on "code", we then correct the hook errors and
        return the data errors.
        
        parity_oct refers to the parity of the weight-8 
        stabilizer, which is measured at a previous step.
        '''

        data_qs = self.circuits[circuit].data_qubits()
        data_q_ids = [q.qubit_id for q in data_qs]
        pre_ns = min(data_q_ids)
        pre_Is = ['I' for i in range(pre_ns)]
        post_ns = len(self.stabs) - max(data_q_ids) - 1
        post_Is = ['I' for i in range(post_ns)]

        output_dict = self.run_one_circ(circuit)
        n_first_anc = min(output_dict.keys())
        data_errors, hook_errors = qfun.stabs_QEC_diVin(output_dict,
                                                        n_first_anc,
                                                        code,
                                                        stab_kind,
                                                        False,
                                                        parity_oct)

        if hook_errors.count('I') != len(hook_errors):
            hook_errors = pre_Is + hook_errors + post_Is
            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           hook_errors)
       
            self.stabs, self.destabs = corr_state[0][:], corr_state[1][:]

        return data_errors





class Measure_2_logicals(Quantum_Operation):
    '''
    Measurement of two logical operators for lattice surgery
    in the color code.
    Inherits from class Quantum_Operation.
    '''
    
    def run_all(self, stab_kind):
        '''
        Run the whole circuit.
        Currently just for a distance-3 code.
        
        Basic idea: we need to measure two operators in a FT way.
        We measured them twice.  If the outcomes coincide, we stop.
        If they're different we measured that operator again.
        The assumption is that there are 6 circuits.
        '''
       
        QEC_func = self.run_one_diVincenzo
        code = 'Steane'

        #brow.from_circuit(self.circuits[1])
        n_subcircs = len(self.circuits)
        #print n_subcircs
        #sys.exit(0)


        
        # (1) Measure one time stabilizers to define them
        # (only used for Mxx if we start with anc in product state of |0>)
        if len(self.circuits) == 2:
            data_errors0 = QEC_func(0, code, stab_kind)
        
            data_qs = self.circuits[0].data_qubits()
            data_q_ids = [q.qubit_id for q in data_qs]
            pre_ns = min(data_q_ids)
            pre_Is = ['I' for i in range(pre_ns)]
            post_ns = len(self.stabs) - max(data_q_ids) - 1
            post_Is = ['I' for i in range(post_ns)]
            
            if data_errors0.count('I') != len(data_errors0):
                data_errors0 = pre_Is + data_errors0 + post_Is
                corr_state = qfun.update_stabs(self.stabs,
                                               self.destabs,
                                               data_errors0)
       
                self.stabs, self.destabs = corr_state[0][:], corr_state[1][:]
        
            
            # (2) Measure M one time
            only_M = self.run_one_circ(1).values()[0][0]
    
            return only_M, 'normal'

        elif len(self.circuits) == 1:
            
            return self.run_one_circ(0).values()[0][0], 'normal'



        # (1) Measure M first time
        first_M = self.run_one_circ(0).values()[0][0]
        #print 'First M =', first_M
        
        # (2) Measure stabilizers first time
        data_errors1 = QEC_func(1, code, stab_kind)
        #print 'data errors 1 =', data_errors1

        # (3) Measure M second time
        second_M = self.run_one_circ(2).values()[0][0]
        #print 'Second M =', second_M

        if data_errors1.count('I') == 7:
            if first_M == second_M:
                return first_M, 'normal'
            else:
                third_M = self.run_one_circ(4).values()[0][0]
                return third_M, 'normal'


        else:
            
            # (4) Measure stabilizers second time
            data_errors2 = QEC_func(3, code, stab_kind)
            #print 'data errors 2 =', data_errors2

            if data_errors1 == data_errors2:
                if first_M == second_M:
                    if stab_kind == 'X':
                        err_i = data_errors1.index('Z')
                        if err_i in [1,3,5]:  
                            corr_type = 'alternative'
                        else:  
                            corr_type = 'normal'

                    elif stab_kind == 'Z':
                        err_i = data_errors1.index('X')
                        if err_i in [0,3,4]:
                            corr_type = 'alternative'
                        else:
                            corr_type = 'normal'

                    return first_M, corr_type

                else:
                    third_M = self.run_one_circ(4).values()[0][0]
                    return third_M, 'normal'

            else:
                
                if first_M == second_M:
                    return first_M, 'normal'

                else:
                    third_M = self.run_one_circ(4).values()[0][0]
                    return third_M, 'normal'
                    
    
    
    def run_all_long(self, stab_kind):
        '''
        Run the whole circuit.
        Currently just for a distance-3 code.
        
        Basic idea: we need to measure two operators in a FT way.
        We measured them twice.  If the outcomes coincide, we stop.
        If they're different we measured that operator again.
        The assumption is that there are 6 circuits.
        '''
       
        code = 'Steane'
        QEC_func = self.run_one_diVincenzo
        if stab_kind == 'X':    
            tricky_indices = [1,3,5]
            Pauli_error = 'Z'
        elif stab_kind == 'Z':    
            tricky_indices = [0,3,4]
            Pauli_error = 'X'
        
        # default values
        corr_nonanc = 'normal'
        corr_anc = False

        #brow.from_circuit(self.circuits[0])
        #n_subcircs = len(self.circuits)
        #print n_subcircs
        #sys.exit(0)

        # list of outcomes of the low-weight and high-weight operators
        low_w, high_w = [], []

        first_subcirc_i = 0

        # (1) Measure one time stabilizers to define them
        # (only used for Mxx if we start with anc in product state of |0>)
        if len(self.circuits) > 12:
            data_errors0 = QEC_func(0, code, stab_kind)
        
            data_qs = self.circuits[0].data_qubits()
            data_q_ids = [q.qubit_id for q in data_qs]
            pre_ns = min(data_q_ids)
            pre_Is = ['I' for i in range(pre_ns)]
            post_ns = len(self.stabs) - max(data_q_ids) - 1
            post_Is = ['I' for i in range(post_ns)]
            
            if data_errors0.count('I') != len(data_errors0):
                data_errors0 = pre_Is + data_errors0 + post_Is
                corr_state = qfun.update_stabs(self.stabs,
                                               self.destabs,
                                               data_errors0)
       
                self.stabs, self.destabs = corr_state[0][:], corr_state[1][:]

            first_subcirc_i = 1


        #print 'State after first projection:'
        #print self.stabs

        # (2) Measure low-weight operator first time
        low_w += [self.run_one_circ(first_subcirc_i).values()[0][0]]
        #print 'low_w1 =', low_w[0]
        
        # (3) Measure high-weight operator first time
        high_w += [self.run_one_circ(first_subcirc_i+1).values()[0][0]]
        #print 'high_w1 =', high_w[0]

        #print 'State after high-w operator:'
        #print self.stabs

        # (4) Measure stabilizers for non-anc first time
        #print 'stab kind =', stab_kind
        data_errors_nonanc1 = QEC_func(first_subcirc_i+2, code, stab_kind)
        #print 'data errors non anc 1 =', data_errors_nonanc1

        # (5) Measure stabilizers for anc first time
        data_errors_anc1 = QEC_func(first_subcirc_i+3, code, stab_kind)
        #print 'data errors anc 1 =', data_errors_anc1

        #print 'State after QEC1:'
        #print self.stabs


        # (6) Measure low-weight operutor second time
        low_w += [self.run_one_circ(first_subcirc_i+4).values()[0][0]]
        #print 'low_w2 =', low_w[1]
        
        # (7) Measure high-weight operator second time
        high_w += [self.run_one_circ(first_subcirc_i+5).values()[0][0]]
        #print 'high_w2 =', high_w[1]

        #print 'State after second high-w operator:'
        #print self.stabs

        
        # if error on non-ancillary qubit
        if data_errors_nonanc1.count('I') != 7:
            
            # run QEC(non-anc) again
            data_errors_nonanc2 = QEC_func(first_subcirc_i+6, code, stab_kind)
            #print 'data errors nonanc 2 =', data_errors_nonanc2

            # if QEC(non-anc)(1) is equal to QEC(non-anc)(2)
            if data_errors_nonanc1 == data_errors_nonanc2:

                # if the operators are the same, then the error happened before
                if (low_w[0]==low_w[1] and high_w[0]==high_w[1]):
                    err_i = data_errors_nonanc1.index(Pauli_error)
                    # if the faulty qubit is on the boundary 
                    if err_i in tricky_indices:  corr_nonanc = 'alternative'


        # if error on ancillary qubit
        if data_errors_anc1.count('I') != 7:
                
            # run QEC(anc) again
            data_errors_anc2 = QEC_func(first_subcirc_i+7, code, stab_kind)
            #print 'data errors anc 2 =', data_errors_anc2

            # if QEC(anc)(1) is equal to QEC(anc)(2)
            if data_errors_anc1 == data_errors_anc2:
                    
                # if the operators are the same, then the error happened before
                if (low_w[0]==low_w[1] and high_w[0]==high_w[1]):
                    err_i = data_errors_anc1.index(Pauli_error)
                    # if the faulty qubit is on the boundary 
                    if err_i in tricky_indices:  corr_anc = True

        if low_w[0] != low_w[1]:
            low_w += [self.run_one_circ(first_subcirc_i+8).values()[0][0]]
            #print 'low_w3 =', low_w[2]
        if high_w[0] != high_w[1]:
            high_w += [self.run_one_circ(first_subcirc_i+9).values()[0][0]]
            #print 'high_w3 =', high_w[2]

        # for the parity, select the last outcome of each operator
        parity = (low_w[-1] + high_w[-1])%2
        return parity, (len(low_w), len(high_w)), corr_nonanc, corr_anc
          

    
    def run_boundary_oper_flags_short_ion(self, error_det=False):
        '''
        Short version:  we run X2 (Z2) first two or three times and then
        X4 (Z4) one (or two or three times).
        This is adapted for the scheduled circuit for ion traps
        '''
    
        first_subcirc_i = 0
        subcircs_run = []
        subcircs_run_QECnonanc = {}
        subcircs_run_QECanc = {}
        for i in range(10):
            subcircs_run_QECnonanc[i] = 0
            subcircs_run_QECanc[i] = 0
        

        # list of outcomes of the low-weight and high-weight operators
        low_w, high_w = [], []
            
        # (1) Measure the low-weight operator
        # We need to flip the outcome because there are 2 MS gates
        low_w += [(self.run_one_circ(first_subcirc_i).values()[0][0]+1)%2]
        subcircs_run += [first_subcirc_i] 
        
        if error_det:
            # (2) Transport ion back to ancilla arm
            self.run_one_circ(first_subcirc_i+2)
            subcircs_run += [first_subcirc_i+2] 

            # (3) Measure the high-weight operator
            # In this case, we don't flip the outcome because there are 4 MS gates
            high_w += [self.run_one_circ(first_subcirc_i+4).values()[0][0]]
            subcircs_run += [first_subcirc_i+4] 

            # Define the dictionary with the subcircs that were run
            subcircs_run_dict = {}
            for i in range(13):
                if i==9:
                    subcircs_run_dict[i] = subcircs_run_QECnonanc
                elif i==10:
                    subcircs_run_dict[i] = subcircs_run_QECanc
                else:
                    if i in subcircs_run:
                        subcircs_run_dict[i] = 1
                    else:
                        subcircs_run_dict[i] = 0

            # If an error has already been detected in a previous step of
            # the supra-circuit, we only need to measure the two operators
            # once.
            output = [low_w, high_w, [0,0,0], [0,0,0]] 
            output += [[0,0,0], [0,0,0], 'normal', 'normal', True, subcircs_run_dict]
            return tuple(output)

        
        # (2) Measure low-weight operator second time
        # We need to flip the outcome because there are 2 MS gates
        low_w += [(self.run_one_circ(first_subcirc_i+1).values()[0][0]+1)%2]
        subcircs_run += [first_subcirc_i+1] 

        # If the two outcomes are equal, we now measure the high-weight
        if low_w[0]==low_w[1]:
           
            # (3) Transport ion back to ancilla arm
            self.run_one_circ(first_subcirc_i+2)
            subcircs_run += [first_subcirc_i+2] 

            # (4) Measure the high-weight operator twice
            high_w += [self.run_one_circ(first_subcirc_i+5).values()[0][0]]
            subcircs_run += [first_subcirc_i+5] 
            high_w += [self.run_one_circ(first_subcirc_i+6).values()[0][0]]
            subcircs_run += [first_subcirc_i+6] 
            
            # If the two outcomes are equal, we now do QEC on both logical qubits
            if high_w[0]==high_w[1]:

                #print 'Im here!'

                # (5) Transport two ions back to ancilla arm
                self.run_one_circ(first_subcirc_i+7)
                subcircs_run += [first_subcirc_i+7] 

                # (6) First QEC (FT) on target (if X) or control (if Z)
                QEC_supracirc = self.circuits[first_subcirc_i+9]
                subcircs_run += [first_subcirc_i+9] 
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_nonanc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_nonanc = QEC_nonanc.run_QEC_FT_lattsurg_ion()
                s_targ1, f_targ1, corr_targ, error_det_targ, subc_nonanc = outputQEC_nonanc
                for subcirc_i in subc_nonanc:
                    subcircs_run_QECnonanc[subcirc_i] = 1
                self.stabs, self.destabs = QEC_nonanc.stabs[:], QEC_nonanc.destabs[:]
       
                # (6) First QEC (FT) on ancilla 
                QEC_supracirc = self.circuits[first_subcirc_i+10]
                subcircs_run += [first_subcirc_i+10] 
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg_ion()
                s_anc1, f_anc1, corr_anc, error_det_anc, subc_anc = outputQEC_anc
                for subcirc_i in subc_anc:
                    subcircs_run_QECanc[subcirc_i] = 1
                self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]
    
                #print 's_targ1 =', s_targ1
                #print 'corr targ =', corr_targ
                #print 's_anc1 =', s_anc1
                #print 'corr anc =', corr_anc 

                #print 'f_targ1 =', f_targ1
                #print 'f_anc1 =', f_anc1

                if corr_targ=='normal' and corr_anc=='normal':
                    error_det_total = error_det_targ or error_det_anc

                else:
                    error_det_total = True
                    
                    clause_targ = (s_targ1[1]==0 and s_targ1[2]==1)
                    clause_anc = (s_anc1[1]==0 and s_anc1[2]==1)
                    clause_flag = (sum(f_targ1)==0 and sum(f_anc1)==0)

                    #print clause_targ
                    #print clause_anc
                    #print clause_flag

                    # if the error happened on one of the qubits involved in X2, then
                    # only re-measure X2.
                    if (clause_targ or clause_anc) and clause_flag:
                        # (7) Measure low-weight operator third time
                        low_w_3rd = self.run_one_circ(first_subcirc_i+11).values()[0][0]
                        subcircs_run += [first_subcirc_i+11] 
                        low_w += [(low_w_3rd+1)%2]
                        if low_w[2]!=low_w[1]:
                            corr_targ, corr_anc = 'normal', 'normal'
                            # the correct eigenvalue of low_w was the first one.
                            low_w[2] = low_w[0]

                        else:
                            if corr_targ=='unknown':  corr_targ = 'alternative'
                            if corr_anc=='unknown':  corr_anc = 'alternative'
                    
                    # else re-measure only X4
                    else:
                        # (8) Measure high-weight operator third time
                        high_w += [self.run_one_circ(first_subcirc_i+12).values()[0][0]]
                        subcircs_run += [first_subcirc_i+12] 
                    
                        #print high_w
                        if high_w[2]!=high_w[1]:
                            corr_targ, corr_anc = 'normal', 'normal'
                            # the correct eigenvalue of low_w was the first one.
                            high_w[2] = high_w[0]
                        
                        else:
                            if corr_targ=='unknown':  corr_targ = 'alternative'
                            if corr_anc=='unknown':  corr_anc = 'alternative'
                

            # If X4 (Z4) changed, then re-measure X4 (Z4).
            elif high_w[0]!=high_w[1]:
                error_det_total = True
                corr_targ, corr_anc = 'normal', 'normal'
                high_w += [self.run_one_circ(first_subcirc_i+8).values()[0][0]]
                subcircs_run += [first_subcirc_i+8] 
            
                # if the second and the third are the same, we still cannot distinguish between 
                # a measurement error and a data error.  We measure the last stabilizer 
                # of both logical qubits nonFT
                if high_w[1]==high_w[2]:
                    # QEC on target or control
                    QEC_supracirc = self.circuits[first_subcirc_i+9]
                    subcircs_run += [first_subcirc_i+9] 
                    QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                    QEC_targ = QEC_with_flags([self.stabs, self.destabs], 
                                               QEC_circs, self.chp_loc)
                    outputQEC_targ = QEC_targ.run_QEC_FT_lattsurg_ion(stab=1)
                    s_targ1, f_targ1, corr_targ, error_det_targ, subc_targ = outputQEC_targ
                    for subcirc_i in subc_targ:
                        subcircs_run_QECnonanc[subcirc_i] = 1
                    self.stabs, self.destabs = QEC_targ.stabs[:], QEC_targ.destabs[:]
       
                    # QEC on ancilla 
                    QEC_supracirc = self.circuits[first_subcirc_i+10]
                    subcircs_run += [first_subcirc_i+10] 
                    QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                    QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                    outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg_ion(stab=1)
                    s_anc1, f_anc1, corr_anc, error_det_anc, subc_anc = outputQEC_anc
                    for subcirc_i in subc_anc:
                        subcircs_run_QECanc[subcirc_i] = 1
                    self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]

                    # if an error was detected, then trust the first high-weight value
                    if s_targ1[1]==1 or s_anc1[1]==1:
                        high_w[2] = high_w[0]
            
                    #print 'f_targ1 =', f_targ1
                    #print 'f_anc1 =', f_anc1
            
                else:
                    s_targ1, f_targ1 = [0,0,0], [0,0,0]
                    s_anc1, f_anc1 = [0,0,0], [0,0,0]

        
        # If X2 (Z2) changed, measure X2 (Z2) a third time and measure X4 (Z4) one time only 
        elif low_w[0]!=low_w[1]:
            error_det_total = True
            corr_targ, corr_anc = 'normal', 'normal'
            low_w += [(self.run_one_circ(first_subcirc_i+3).values()[0][0]+1)%2]
            subcircs_run += [first_subcirc_i+3] 
            
            # Measure the high-weight operator only once
            high_w += [self.run_one_circ(first_subcirc_i+4).values()[0][0]]
            subcircs_run += [first_subcirc_i+4] 

            # if the second and the third for the low-weight are the same, we still 
            # cannot distinguish between a measurement error and a data error.  
            # We measure the last stabilizer of both logical qubits nonFT
            if low_w[1]==low_w[2]:
                # QEC on target or control
                QEC_supracirc = self.circuits[first_subcirc_i+9]
                subcircs_run += [first_subcirc_i+9] 
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_targ = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_targ = QEC_targ.run_QEC_FT_lattsurg_ion(stab=2)
                s_targ1, f_targ1, corr_targ, error_det_targ, subc_targ = outputQEC_targ
                for subcirc_i in subc_targ:
                    subcircs_run_QECnonanc[subcirc_i] = 1
                self.stabs, self.destabs = QEC_targ.stabs[:], QEC_targ.destabs[:]
       
                # QEC on ancilla 
                QEC_supracirc = self.circuits[first_subcirc_i+10]
                subcircs_run += [first_subcirc_i+10] 
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg_ion(stab=2)
                s_anc1, f_anc1, corr_anc, error_det_anc, subc_anc = outputQEC_anc
                for subcirc_i in subc_anc:
                    subcircs_run_QECanc[subcirc_i] = 1
                self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]

                # if an error was detected, then trust the first X2 value
                if s_targ1[2]==1 or s_anc1[2]==1:
                    low_w[2] = low_w[0]
            
                #print 'f_targ1 =', f_targ1
                #print 'f_anc1 =', f_anc1

            else:
                s_targ1, f_targ1 = [0,0,0], [0,0,0]
                s_anc1, f_anc1 = [0,0,0], [0,0,0]
            
        
        # Define the dictionary with the subcircs that were run
        subcircs_run_dict = {}
        for i in range(13):
            if i==9:
                subcircs_run_dict[i] = subcircs_run_QECnonanc
            elif i==10:
                subcircs_run_dict[i] = subcircs_run_QECanc
            else:
                if i in subcircs_run:
                    subcircs_run_dict[i] = 1
                else:
                    subcircs_run_dict[i] = 0

        output = [low_w, high_w, [s_targ1], [s_anc1]] 
        output += [f_targ1, f_anc1, corr_targ, corr_anc, error_det_total, subcircs_run_dict]
        return tuple(output)



    def run_boundary_oper_flags_short(self, error_det=False):
        '''
        Short version:  we run X2 (Z2) first two or three times and then
        X4 (Z4) one (or two or three times)
        '''
    
        first_subcirc_i = 0

        # list of outcomes of the low-weight and high-weight operators
        low_w, high_w = [], []

        # (2) Measure low-weight operator first time
        low_w += [self.run_one_circ(first_subcirc_i).values()[0][0]]
        #print 'low_w1 =', low_w[0]
       
        if error_det:
            # (3) Measure high-weight operator first time
            high_w += [self.run_one_circ(first_subcirc_i+1).values()[0][0]]
            #print 'high_w1 =', high_w[0]
       
            # If an error has already been detected in a previous step of
            # the supra-circuit, we only need to measure the two operators
            # once.
            output = [low_w, high_w, [0,0,0], [0,0,0]] 
            output += [[0,0,0], [0,0,0], 'normal', 'normal', True]
            return tuple(output)

        # (3) Measure low-weight operator second time
        low_w += [self.run_one_circ(first_subcirc_i+2).values()[0][0]]

        # If the two outcomes are equal, we now measure the high-weight
        if low_w[0]==low_w[1]:
            
            # (4) Measure the high-weight operator twice
            high_w += [self.run_one_circ(first_subcirc_i+1).values()[0][0]]
            high_w += [self.run_one_circ(first_subcirc_i+3).values()[0][0]]
            
            # If the two outcomes are equal, we now do QEC on both logical qubits
            if high_w[0]==high_w[1]:

                #print 'Im here!'

                # (5) First QEC (FT) on target (if X) or control (if Z)
                QEC_supracirc = self.circuits[first_subcirc_i+4]
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_nonanc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_nonanc = QEC_nonanc.run_QEC_FT_lattsurg()
                s_targ1, f_targ1, corr_targ, error_det_targ = outputQEC_nonanc
                self.stabs, self.destabs = QEC_nonanc.stabs[:], QEC_nonanc.destabs[:]
       
                # (6) First QEC (FT) on ancilla 
                QEC_supracirc = self.circuits[first_subcirc_i+5]
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg()
                s_anc1, f_anc1, corr_anc, error_det_anc = outputQEC_anc
                self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]
    
                #print 's_targ1 =', s_targ1
                #print 'corr targ =', corr_targ
                #print 's_anc1 =', s_anc1
                #print 'corr anc =', corr_anc 

                #print 'f_targ1 =', f_targ1
                #print 'f_anc1 =', f_anc1

                if corr_targ=='normal' and corr_anc=='normal':
                    error_det_total = error_det_targ or error_det_anc

                else:
                    error_det_total = True
                    
                    clause_targ = (s_targ1[1]==0 and s_targ1[2]==1)
                    clause_anc = (s_anc1[1]==0 and s_anc1[2]==1)
                    clause_flag = (sum(f_targ1)==0 and sum(f_anc1)==0)

                    #print clause_targ
                    #print clause_anc
                    #print clause_flag

                    # if the error happened on one of the qubits involved in X2, then
                    # only re-measure X2.
                    if (clause_targ or clause_anc) and clause_flag:
                        # (7) Measure low-weight operator third time
                        low_w += [self.run_one_circ(first_subcirc_i+6).values()[0][0]]
                        if low_w[2]!=low_w[1]:
                            corr_targ, corr_anc = 'normal', 'normal'
                            # the correct eigenvalue of low_w was the first one.
                            low_w[2] = low_w[0]

                        else:
                            if corr_targ=='unknown':  corr_targ = 'alternative'
                            if corr_anc=='unknown':  corr_anc = 'alternative'
                    
                    # else re-measure only X4
                    else:
                        # (8) Measure high-weight operator third time
                        high_w += [self.run_one_circ(first_subcirc_i+7).values()[0][0]]
                    
                        #print high_w
                        if high_w[2]!=high_w[1]:
                            corr_targ, corr_anc = 'normal', 'normal'
                            # the correct eigenvalue of low_w was the first one.
                            high_w[2] = high_w[0]
                        
                        else:
                            if corr_targ=='unknown':  corr_targ = 'alternative'
                            if corr_anc=='unknown':  corr_anc = 'alternative'
                

            # If X4 (Z4) changed, then re-measure X4 (Z4).
            elif high_w[0]!=high_w[1]:
                error_det_total = True
                corr_targ, corr_anc = 'normal', 'normal'
                high_w += [self.run_one_circ(first_subcirc_i+7).values()[0][0]]
            
                # if the second and the third are the same, we still cannot distinguish between 
                # a measurement error and a data error.  We measure the last stabilizer 
                # of both logical qubits nonFT
                if high_w[1]==high_w[2]:
                    # QEC on target or control
                    QEC_supracirc = self.circuits[first_subcirc_i+4]
                    QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                    QEC_targ = QEC_with_flags([self.stabs, self.destabs], 
                                               QEC_circs, self.chp_loc)
                    outputQEC_targ = QEC_targ.run_QEC_FT_lattsurg(stab=1)
                    s_targ1, f_targ1, corr_targ, error_det_targ = outputQEC_targ
                    self.stabs, self.destabs = QEC_targ.stabs[:], QEC_targ.destabs[:]
       
                    # QEC on ancilla 
                    QEC_supracirc = self.circuits[first_subcirc_i+5]
                    QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                    QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                    outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg(stab=1)
                    s_anc1, f_anc1, corr_anc, error_det_anc = outputQEC_anc
                    self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]

                    # if an error was detected, then trust the first high-weight value
                    if s_targ1[1]==1 or s_anc1[1]==1:
                        high_w[2] = high_w[0]
            
                    #print 'f_targ1 =', f_targ1
                    #print 'f_anc1 =', f_anc1
            
                else:
                    s_targ1, f_targ1 = [0,0,0], [0,0,0]
                    s_anc1, f_anc1 = [0,0,0], [0,0,0]

        
        # If X2 (Z2) changed, measure X2 (Z2) a third time and measure X4 (Z4) one time only 
        elif low_w[0]!=low_w[1]:
            error_det_total = True
            corr_targ, corr_anc = 'normal', 'normal'
            low_w += [self.run_one_circ(first_subcirc_i+6).values()[0][0]]
            
            # Measure the high-weight operator only once
            high_w += [self.run_one_circ(first_subcirc_i+1).values()[0][0]]

            # if the second and the third for the low-weight are the same, we still 
            # cannot distinguish between a measurement error and a data error.  
            # We measure the last stabilizer of both logical qubits nonFT
            if low_w[1]==low_w[2]:
                # QEC on target or control
                QEC_supracirc = self.circuits[first_subcirc_i+4]
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_targ = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_targ = QEC_targ.run_QEC_FT_lattsurg(stab=2)
                s_targ1, f_targ1, corr_targ, error_det_targ = outputQEC_targ
                self.stabs, self.destabs = QEC_targ.stabs[:], QEC_targ.destabs[:]
       
                # QEC on ancilla 
                QEC_supracirc = self.circuits[first_subcirc_i+5]
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg(stab=2)
                s_anc1, f_anc1, corr_anc, error_det_anc = outputQEC_anc
                self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]

                # if an error was detected, then trust the first X2 value
                if s_targ1[2]==1 or s_anc1[2]==1:
                    low_w[2] = low_w[0]
            
                #print 'f_targ1 =', f_targ1
                #print 'f_anc1 =', f_anc1

            else:
                s_targ1, f_targ1 = [0,0,0], [0,0,0]
                s_anc1, f_anc1 = [0,0,0], [0,0,0]

        output = [low_w, high_w, [s_targ1], [s_anc1]] 
        output += [f_targ1, f_anc1, corr_targ, corr_anc, error_det_total]
        return tuple(output)



    def run_boundary_oper_flags(self, error_det=False):
        '''
        runs the circuit to measure the XX or ZZ between two logical qubits in
        a FT fashion.
        '''    
        
        #tricky_indices = [4,5,6]  # using Alejandro's indexing
        
        # default values
        #corr_nonanc = 'normal'
        #corr_anc = False

        print 'error det =', error_det 

        # For now we will assume that there is no QEC step on the ancilla before
        first_subcirc_i = 0

        # list of outcomes of the low-weight and high-weight operators
        low_w, high_w = [], []

        # (2) Measure low-weight operator first time
        low_w += [self.run_one_circ(first_subcirc_i).values()[0][0]]
        print 'low_w1 =', low_w[0]
        
        # (3) Measure high-weight operator first time
        high_w += [self.run_one_circ(first_subcirc_i+1).values()[0][0]]
        print 'high_w1 =', high_w[0]
       
        # If an error has already been detected in a previous step of
        # the supra-circuit, we only need to measure the two operators
        # once.
        if error_det:
            output = [low_w, high_w, [0,0,0], [0,0,0]] 
            output += [[0,0,0], [0,0,0], 'normal', 'normal', True]
            return tuple(output)

        # (4) Measure low-weight operator second time
        low_w += [self.run_one_circ(first_subcirc_i+2).values()[0][0]]
        print 'low_w2 =', low_w[1]

        # (5) Measure high-weight operator second time
        high_w += [self.run_one_circ(first_subcirc_i+3).values()[0][0]]
        print 'high_w2 =', high_w[1]

        #print 'State after high-w operator:'
        #print self.stabs

        # We only measure the X stabilizers if the operators' eigenvalues did 
        # not change.
        if (low_w[0]==low_w[1]) and (high_w[0]==high_w[1]):
    
            #print 'Im here!'

            # (6) First QEC (FT) on target (if X) or control (if Z)
            QEC_supracirc = self.circuits[first_subcirc_i+4]
            QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
            QEC_nonanc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
            outputQEC_nonanc = QEC_nonanc.run_QEC_FT_lattsurg()
            s_targ1, f_targ1, corr_targ, error_det_targ = outputQEC_nonanc
            self.stabs, self.destabs = QEC_nonanc.stabs[:], QEC_nonanc.destabs[:]
       
            # (7) First QEC (FT) on ancilla 
            QEC_supracirc = self.circuits[first_subcirc_i+5]
            QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
            QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
            outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg()
            s_anc1, f_anc1, corr_anc, error_det_anc = outputQEC_anc
            self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]
    
            print 's_targ1 =', s_targ1
            #print 'corr targ =', corr_targ
            print 's_anc1 =', s_anc1
            #print 'corr anc =', corr_anc 

            print 'f_targ1 =', f_targ1
            print 'f_anc1 =', f_anc1

            if corr_targ=='normal' and corr_anc=='normal':
                error_det_total = error_det_targ or error_det_anc

            else:
                error_det_total = True
                    
                clause_targ = (s_targ1[1]==0 and s_targ1[2]==1)
                clause_anc = (s_anc1[1]==0 and s_anc1[2]==1)
                clause_flag = (sum(f_targ1)==0 and sum(f_anc1)==0)

                #print clause_targ
                #print clause_anc
                #print clause_flag

                # if the error happened on one of the qubits involved in X2, then
                # only re-measure X2.
                if (clause_targ or clause_anc) and clause_flag:
                    # (8) Measure low-weight operator third time
                    low_w += [self.run_one_circ(first_subcirc_i+6).values()[0][0]]
                    if low_w[2]!=low_w[1]:
                        corr_targ, corr_anc = 'normal', 'normal'
                        # the correct eigenvalue of low_w was the first one.
                        low_w[2] = low_w[0]

                    else:
                        if corr_targ=='unknown':  corr_targ = 'alternative'
                        if corr_anc=='unknown':  corr_anc = 'alternative'
                    
                # else re-measure only X4
                else:
                    # (9) Measure high-weight operator third time
                    high_w += [self.run_one_circ(first_subcirc_i+7).values()[0][0]]
                    
                    print high_w
                    if high_w[2]!=high_w[1]:
                        corr_targ, corr_anc = 'normal', 'normal'
                        # the correct eigenvalue of low_w was the first one.
                        high_w[2] = high_w[0]
                        
                    else:
                        if corr_targ=='unknown':  corr_targ = 'alternative'
                        if corr_anc=='unknown':  corr_anc = 'alternative'
            
            #print corr_targ
            #print corr_anc
                
        # If the parities of the boundary operators changed, we do not have to do
        # QEC because we can safely assume that the error occurred after the 
        # first round of boundary measurements.
       

        # If only X2 changed, then only re-measure X2.
        elif (low_w[0]!=low_w[1]) and (high_w[0]==high_w[1]):
            error_det_total = True
            corr_targ, corr_anc = 'normal', 'normal'
            low_w += [self.run_one_circ(first_subcirc_i+6).values()[0][0]]
            
            # if the second and the third are the same, we still cannot distinguish between 
            # a measurement error and a data error.  We measure the last stabilizer 
            # of both logical qubits nonFT
            if low_w[1]==low_w[2]:
                # QEC on target or control
                QEC_supracirc = self.circuits[first_subcirc_i+4]
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_targ = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_targ = QEC_targ.run_QEC_FT_lattsurg(stab=2)
                s_targ1, f_targ1, corr_targ, error_det_targ = outputQEC_targ
                self.stabs, self.destabs = QEC_targ.stabs[:], QEC_targ.destabs[:]
       
                # QEC on ancilla 
                QEC_supracirc = self.circuits[first_subcirc_i+5]
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg(stab=2)
                s_anc1, f_anc1, corr_anc, error_det_anc = outputQEC_anc
                self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]

                # if an error was detected, then trust the first X2 value
                if s_targ1[2]==1 or s_anc1[2]==1:
                    low_w[2] = low_w[0]
            
                #print 'f_targ1 =', f_targ1
                #print 'f_anc1 =', f_anc1

            else:
                s_targ1, f_targ1 = [0,0,0], [0,0,0]
                s_anc1, f_anc1 = [0,0,0], [0,0,0]


        # If only X4 changed, then only re-measure X4.
        elif (low_w[0]==low_w[1]) and (high_w[0]!=high_w[1]):
            error_det_total = True
            corr_targ, corr_anc = 'normal', 'normal'
            high_w += [self.run_one_circ(first_subcirc_i+7).values()[0][0]]
            
            # if the second and the third are the same, we still cannot distinguish between 
            # a measurement error and a data error.  We measure the last stabilizer 
            # of both logical qubits nonFT
            if high_w[1]==high_w[2]:
                # QEC on target or control
                QEC_supracirc = self.circuits[first_subcirc_i+4]
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_targ = QEC_with_flags([self.stabs, self.destabs], 
                                           QEC_circs, self.chp_loc)
                outputQEC_targ = QEC_targ.run_QEC_FT_lattsurg(stab=1)
                s_targ1, f_targ1, corr_targ, error_det_targ = outputQEC_targ
                self.stabs, self.destabs = QEC_targ.stabs[:], QEC_targ.destabs[:]
       
                # QEC on ancilla 
                QEC_supracirc = self.circuits[first_subcirc_i+5]
                QEC_circs = [gate.circuit_list[0] for gate in QEC_supracirc.gates]
                QEC_anc = QEC_with_flags([self.stabs, self.destabs], QEC_circs, self.chp_loc)
                outputQEC_anc = QEC_anc.run_QEC_FT_lattsurg(stab=1)
                s_anc1, f_anc1, corr_anc, error_det_anc = outputQEC_anc
                self.stabs, self.destabs = QEC_anc.stabs[:], QEC_anc.destabs[:]

                # if an error was detected, then trust the first X2 value
                if s_targ1[1]==1 or s_anc1[1]==1:
                    high_w[2] = high_w[0]
            
                #print 'f_targ1 =', f_targ1
                #print 'f_anc1 =', f_anc1
            
            else:
                s_targ1, f_targ1 = [0,0,0], [0,0,0]
                s_anc1, f_anc1 = [0,0,0], [0,0,0]

        # if both X2 and X4 changed, then this was caused by a w-2 error.
        # We don't do anything else.
        else:
            error_det_total = True
            corr_targ, corr_anc = 'normal', 'normal'
            s_targ1, s_anc1 = [0,0,0], [0,0,0]
            f_targ1, f_anc1 = [0,0,0], [0,0,0]

        
        #print 'f_targ1 =', f_targ1
        #print 'f_anc1 =', f_anc1
        output = [low_w, high_w, [s_targ1], [s_anc1]] 
        output += [f_targ1, f_anc1, corr_targ, corr_anc, error_det_total]
        return tuple(output)
            


class QEC_d3(Quantum_Operation):
    '''
    Quantum Error Correction for a distance-3 code.
    Inherits from class Quantum_Operation
    '''

    def run_one_bare_anc(self, circuit, code, stab_kind):
        '''
        runs one round of QEC for the whole set of stabilizers
        assuming bare ancillae.
        circuit refers to the index of the circuit to be run.
        '''

        # need to add stab_kind to include surface-17 and
        # large-distance tological codes.

        output_dict = self.run_one_circ(circuit)
        anc_qubit_list = sorted(output_dict.keys())
        n_first_anc_X = min(anc_qubit_list)
        #print n_first_anc_X
        #print anc_qubit_list
        if stab_kind == 'both':
            
            # X stabilizers come first by convention
            X_dict = {key: output_dict[key] 
                           for key in anc_qubit_list[:3]}
            Z_errors = qfun.stabs_QEC_bare_anc(X_dict,
                                               n_first_anc_X,
                                               code)
            Z_errors = ['Z' if oper=='E' else oper for oper in Z_errors]
        
            # Z stabilizers come second
            n_first_anc_Z = min(anc_qubit_list[3:])
            Z_dict = {key: output_dict[key] 
                           for key in anc_qubit_list[3:]} 
            
            X_errors = qfun.stabs_QEC_bare_anc(Z_dict,
                                               n_first_anc_Z,
                                               code)
            X_errors = ['X' if oper=='E' else oper for oper in X_errors]
       
            data_errors = Z_errors, X_errors

        
        elif stab_kind == 'X':
            
            data_errors = qfun.stabs_QEC_bare_anc(output_dict,
                                               n_first_anc_X,
                                               code,
                                               stab_kind)
            data_errors = ['Z' if oper=='E' else oper for oper in data_errors]
        

        elif stab_kind == 'Z':
            #print output_dict
            data_errors = qfun.stabs_QEC_bare_anc(output_dict,
                                               n_first_anc_X,
                                               code,
                                               stab_kind)
            data_errors = ['X' if oper=='E' else oper for oper in data_errors]

       
        # Added by MGA 12/23/19 to get stabilizer outcome info for the NN decoder        
        stabilizer_outcomes = [output_dict[anc][0] for anc in anc_qubit_list]

        return data_errors, stabilizer_outcomes

    
    
    def run_one_bare_anc_new(self, circuit, code, stab_kind, errors_det=0,
                             list_syndromesX=[], list_syndromesZ=[],
                             added_data_err_previous=[False,False]):
        '''
        runs one round of QEC for the whole set of stabilizers
        assuming bare ancillae.
        circuit refers to the index of the circuit to be run.
        '''

        output_dict = self.run_one_circ(circuit)
        anc_qubit_list = sorted(output_dict.keys())
        n_first_anc = min(anc_qubit_list)
        
        # These next lines are taken from the high indeterminacy function from
        # QEC_with_flags:
        # If a new data errror was detected by the stabilizers, add 1 to errors_det
        # We have to be very careful (conservative).  We only add 1 if the syndrome
        # is different from the last 2 syndromes, if no flags were triggered in the 2
        # previous steps and if no data error was detected on the previous step.
        # I'm sure we can refine it even further.
        
        out_keys = output_dict.keys()[:]
        out_keys.sort()
        syndromes = [output_dict[key][0] for key in out_keys]

        return syndromes

        #if len(list_syndromesX) == 0:
        #    last_syndromeX = [0 for i in range(len(anc_qubit_list))]
        #else:
        #    last_syndromeX = list_syndromesX[-1][:]
        #if len(list_syndromesZ) == 0:
        #    last_syndromeZ = [0 for i in range(len(anc_qubit_list))]
        #else:
        #    last_syndromeZ = list_syndromesZ[-1][:]
        
        # Define logical clauses
        #data_error_det = False
        #syn_clause = (syndromes != last_syndromeX) and (syndromes != last_syndromeZ)
        #data_err_clause = (not added_data_err_previous[0]) and (not added_data_err_previous[1])
        #if (syn_clause) and (data_err_clause):
        #    errors_det += 1
        #    data_error_det = True

        #return syndromes, errors_det, data_error_det



    def run_fullQEC_nonCSS_d3(self, code, bare=True):
        '''
        runs 2 or 3 rounds of QEC for a distance-3 non-CSS code,
        like the 5-qubit code or the Cross code.
        At the end, it applies a correction.
        It returns the number of QEC rounds.
        '''
        
        if bare:  QEC_func = self.run_one_bare_anc
        else:     QEC_func = self.run_one_diVincenzo
        
        data_qs = self.circuits[0].data_qubits()
        data_q_ids = [q.qubit_id for q in data_qs]
        pre_ns = min(data_q_ids)
        pre_Is = ['I' for i in range(pre_ns)]
        post_ns = len(self.stabs) - max(data_q_ids) - 1
        post_Is = ['I' for i in range(post_ns)]

        data_errors = []
        for i in range(2):
            #print 'stabs =', self.stabs
            data_errors += [QEC_func(i, code)]
            #print 'errors =', data_errors

        #print 'stabs after 2 =', self.stabs
        if data_errors[0] != data_errors[1]:
            data_errors += [QEC_func(2, code)]
            #print 'stabs after 3 =', self.stabs

        correction = data_errors[-1]
        #print 'correction =', correction

        # update the final states only if a correction
        # is needed, to save some time
        if correction.count('I') != len(correction):
            correction = pre_Is + correction + post_Is
            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           correction)
            self.stabs = corr_state[0][:]
            self.destabs = corr_state[1][:]

        #print 'stabs after correction =', self.stabs
        
        return len(data_errors)



    def run_fullQEC_CSS_d3(self, code, bare=True, old_dec=True):
        '''
        runs 2 or 3 rounds of QEC for a distance-3 CSS code,
        like the Steane code or the surface-17.
        At the end, it applies a correction.
        It returns the number of QEC rounds.
        It assumes the X stabilizers come first.
        pre_n:  number of physical qubits before the physical
                qubits onto which this QEC is acting.
        post_n: number of physical qubits after the physical
                qubits onto which this QEC is acting.

        old_dec:  True if we are using the old meta-decoder,
                  where we always measure the stabilizers at
                  least twice.
                  False if we are using the new meta-decoder,
                  where we measure the stabilizers at most
                  twice.
        '''

        if bare:  QEC_func = self.run_one_bare_anc
        else:     QEC_func = self.run_one_diVincenzo

        data_qs = self.circuits[0].data_qubits()
        data_q_ids = [q.qubit_id for q in data_qs]
        pre_ns = min(data_q_ids)
        pre_Is = ['I' for i in range(pre_ns)]
        post_ns = len(self.stabs) - max(data_q_ids) - 1
        post_Is = ['I' for i in range(post_ns)]
        #print 'pre =', pre_ns
        #print 'post =', post_ns

        Z_data_errors, X_data_errors = [], []
        X_stab_outcomes, Z_stab_outcomes = [], []

        # this is the case if we are just trying to
        # perfect EC to distinguish between correctable
        # and uncorrectable errors.
        if len(self.circuits) == 1:
            Z_corr, X_corr = QEC_func(0, code, 'both')
            Z_data_errors = [Z_corr]
            X_data_errors = [X_corr]

        else:

            if old_dec:
                # run first 4 subcircuits (X stabs first)
                for i in range(2):
                    #print 'stabs =', self.stabs
                    subcirc_outcome = QEC_func(2*i, code, 'X')
                    Z_data_errors += [subcirc_outcome[0]]
                    X_stab_outcomes += [subcirc_outcome[1]]
                    #print 'Z errors =', Z_data_errors
                    
                    subcirc_outcome = QEC_func(2*i+1, code, 'Z')
                    X_data_errors += [subcirc_outcome[0]]
                    Z_stab_outcomes += [subcirc_outcome[1]]
                    #print 'X errors =', X_data_errors

                #print 'stabs after 2 =', self.stabs

                # if the outcomes of the 2 X stabs measurements
                # don't coincide, do it a third time
                if Z_data_errors[0] != Z_data_errors[1]:
                    subcirc_outcome = QEC_func(4, code, 'X')
                    Z_data_errors += [subcirc_outcome[0]]
                    X_stab_outcomes += [subcirc_outcome[1]]
                  
                # same for the Z stabs
                if X_data_errors[0] != X_data_errors[1]:
                    subcirc_outcome = QEC_func(5, code, 'Z')
                    X_data_errors += [subcirc_outcome[0]]
                    Z_stab_outcomes += [subcirc_outcome[1]]

                #print 'X_errors =', X_data_errors
                #print 'Z_errors =', Z_data_errors

        
            else:
                # run the first 2 subcircuits (X stabs first)
                Z_data_errors += [QEC_func(0, code, 'X')[0]]
                X_data_errors += [QEC_func(1, code, 'Z')[0]]
           
                #print Z_data_errors
                #print X_data_errors

                # only run last 2 if there was an error.
                if 'Z' in Z_data_errors[0] or 'X' in X_data_errors[0]:
                    Z_data_errors += [QEC_func(2, code, 'X')[0]]
                    X_data_errors += [QEC_func(3, code, 'Z')[0]]
                
                    #print Z_data_errors
                    #print X_data_errors

            Z_corr, X_corr = Z_data_errors[-1], X_data_errors[-1]


        # Added to send the uncorrected stabilizers to the NN decoder
        # Since lists are mutable, we do a deep copy.
        # MGA 1/29/2020
        uncorrected_stabs = self.stabs[:]
        uncorrected_destabs = self.destabs[:]

        # update the final states only if a correction
        # is needed, to save some time
        if 'Z' in Z_corr:
            #print 'stabs =', self.stabs
            #print 'Z_corr =', Z_corr
            Z_corr = pre_Is[:] + Z_corr[:] + post_Is[:]
            #print 'Z_corr =', Z_corr
            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           Z_corr)
            self.stabs = corr_state[0][:]
            self.destabs = corr_state[1][:]
        
        if 'X' in X_corr:
            X_corr = pre_Is[:] + X_corr[:] + post_Is[:]
            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           X_corr)
            #print 'X_corr =', X_corr
            self.stabs = corr_state[0][:]
            self.destabs = corr_state[1][:]

        return len(Z_data_errors), len(X_data_errors), X_stab_outcomes, Z_stab_outcomes, Z_corr, X_corr, uncorrected_stabs, uncorrected_destabs

    
    
    def run_fullQEC_CSS_d3_rounds_even(self, code, n_rounds):
        '''
        runs n_rounds(multiple of 3) rounds of QEC for the surface-17.
        Applies a correction after 3 or 2 rounds
        It returns the number of QEC rounds.
        It assumes the X stabilizers come first.
        pre_n:  number of physical qubits before the physical
                qubits onto which this QEC is acting.
        post_n: number of physical qubits after the physical
                qubits onto which this QEC is acting.

        For the evenly spaced, there are Identity blocks in between
        each round of QEC.

        MGA: 2/25/2020.
        
        '''

        QEC_func = self.run_one_bare_anc

        data_qs = self.circuits[0].data_qubits()
        data_q_ids = [q.qubit_id for q in data_qs]
        pre_ns = min(data_q_ids)
        pre_Is = ['I' for i in range(pre_ns)]
        post_ns = len(self.stabs) - max(data_q_ids) - 1
        post_Is = ['I' for i in range(post_ns)]
        #print 'pre =', pre_ns
        #print 'post =', post_ns

        Z_data_errors, X_data_errors = [], []
        X_stab_outcomes, Z_stab_outcomes = [], []


        # Run first subcircuit (Identity)
        first_I_outcome = self.run_one_circ(0)

        # run first 4 subcircuits (X stabs first)
        for i in [1,4]:
            #print 'stabs =', self.stabs
            subcirc_outcome = QEC_func(i, code, 'X')
            Z_data_errors += [subcirc_outcome[0]]
            X_stab_outcomes += [subcirc_outcome[1]]
            #print 'Z errors =', Z_data_errors
                    
            subcirc_outcome = QEC_func(i+1, code, 'Z')
            X_data_errors += [subcirc_outcome[0]]
            Z_stab_outcomes += [subcirc_outcome[1]]
            #print 'X errors =', X_data_errors

            intermediate_I_outcome = self.run_one_circ(i+2)


        #print 'stabs after 2 =', self.stabs

        # if the outcomes of the 2 X stabs measurements
        # don't coincide, do it a third time
        if Z_data_errors[0] != Z_data_errors[1]:
            subcirc_outcome = QEC_func(7, code, 'X')
            Z_data_errors += [subcirc_outcome[0]]
            X_stab_outcomes += [subcirc_outcome[1]]
                  
        # same for the Z stabs
        if X_data_errors[0] != X_data_errors[1]:
            subcirc_outcome = QEC_func(8, code, 'Z')
            X_data_errors += [subcirc_outcome[0]]
            Z_stab_outcomes += [subcirc_outcome[1]]

        #print 'X_errors =', X_data_errors
        #print 'Z_errors =', Z_data_errors
        
        Z_corr, X_corr = Z_data_errors[-1], X_data_errors[-1]

        last_I_outcome = self.run_one_circ(9)

        # Added to send the uncorrected stabilizers to the NN decoder
        # Since lists are mutable, we do a deep copy.
        # MGA 1/29/2020
        uncorrected_stabs = self.stabs[:]
        uncorrected_destabs = self.destabs[:]

        # update the final states only if a correction
        # is needed, to save some time
        if 'Z' in Z_corr:
            #print 'stabs =', self.stabs
            #print 'Z_corr =', Z_corr
            Z_corr = pre_Is[:] + Z_corr[:] + post_Is[:]
            #print 'Z_corr =', Z_corr
            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           Z_corr)
            self.stabs = corr_state[0][:]
            self.destabs = corr_state[1][:]
        
        if 'X' in X_corr:
            X_corr = pre_Is[:] + X_corr[:] + post_Is[:]
            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           X_corr)
            #print 'X_corr =', X_corr
            self.stabs = corr_state[0][:]
            self.destabs = corr_state[1][:]

        return len(Z_data_errors), len(X_data_errors), X_stab_outcomes, Z_stab_outcomes, Z_corr, X_corr, uncorrected_stabs, uncorrected_destabs


    
    def run_fullQEC_CSS_d3_rounds_even_onlyNN(self, code, n_rounds):
        '''
        runs n_rounds(not necessariliy multiple of 3) rounds of QEC for the surface-17.
        Does not apply a correction.
                It assumes the X stabilizers come first.
        pre_n:  number of physical qubits before the physical
                qubits onto which this QEC is acting.
        post_n: number of physical qubits after the physical
                qubits onto which this QEC is acting.

        For the evenly spaced, there are Identity blocks in between
        each round of QEC.

        MGA: 6/9/2020.
        
        '''

        QEC_func = self.run_one_bare_anc

        data_qs = self.circuits[0].data_qubits()
        data_q_ids = [q.qubit_id for q in data_qs]
        pre_ns = min(data_q_ids)
        pre_Is = ['I' for i in range(pre_ns)]
        post_ns = len(self.stabs) - max(data_q_ids) - 1
        post_Is = ['I' for i in range(post_ns)]
        #print 'pre =', pre_ns
        #print 'post =', post_ns

        Z_data_errors, X_data_errors = [], []
        X_stab_outcomes, Z_stab_outcomes = [], []


        # Run first subcircuit (Identity)
        first_I_outcome = self.run_one_circ(0)

        # run all subcircuits (X stabs first)
        for i in range(n_rounds):
            #print 'stabs =', self.stabs
            subcirc_outcome = QEC_func(3*i+1, code, 'X')
            Z_data_errors += [subcirc_outcome[0]]
            X_stab_outcomes += [subcirc_outcome[1]]
            #print 'Z errors =', Z_data_errors
                    
            subcirc_outcome = QEC_func(3*i+2, code, 'Z')
            X_data_errors += [subcirc_outcome[0]]
            Z_stab_outcomes += [subcirc_outcome[1]]
            #print 'X errors =', X_data_errors

            intermediate_I_outcome = self.run_one_circ(3*i+3)

        return len(Z_data_errors), len(X_data_errors), X_stab_outcomes, Z_stab_outcomes    

    

    def run_jointQEC(self, stab_kind):
        '''
        '''

        code = 'Steane'
        QEC_func = self.run_one_diVincenzo
        if stab_kind == 'X':  
            Pauli_error = 'Z'
            error_index = 3
        elif stab_kind == 'Z':  
            Pauli_error = 'X'
            error_index = 1

        data_errors1, data_errors2 = [], []

        for i in range(2):
        
            data_errors1 += [QEC_func(i, code, stab_kind)]
            #print 'errors1 =', data_errors1
            data_errors2 += [QEC_func(i+3, code, stab_kind)]
            #print 'errors2 =', data_errors2

        #if stab_kind == 'Z':
        #    print 'errors1 =', data_errors1
        #    print 'errors2 =', data_errors2

        # if the outcomes of the 2 stabs measurements
        # don't coincide, do it a third time
        if data_errors1[0] != data_errors1[1]:
            data_errors1 += [QEC_func(2, code, stab_kind)]
        
        if data_errors2[0] != data_errors2[1]:
            data_errors2 += [QEC_func(5, code, stab_kind)]

        data_error1, data_error2 = data_errors1[-1], data_errors2[-1]
        w1 = 7 - data_error1.count('I')
        w2 = 7 - data_error2.count('I')

        if (w1+w2) == 2:
            #if stab_kind == 'Z':
            #    print 'Im here'
            data_error1 = ['I' for i in range(7)]
            data_error1[error_index] = Pauli_error
            data_error2 = ['I' for i in range(7)]
            data_error2[error_index] = Pauli_error
            
            #if data_error1[error_index] == 'I':
            #    bin_error1 = [0 if oper=='I' else 1 for oper in data_error1]
            #    bin_error1[error_index] = 1
            #    print bin_error1
            #    print st.decode_meas_Steane(bin_error1)

            if stab_kind == 'Z':
                log_corr = ['I' for i in range(7)] + data_error1 + data_error2
            elif stab_kind == 'X':
                log_corr = data_error1 + ['I' for i in range(7)] + data_error2

            #print 'State before JointQEC%s:' %stab_kind
            #print self.stabs

            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           log_corr)
            #print 'log_corr =', log_corr
            self.stabs = corr_state[0][:]
            self.destabs = corr_state[1][:]

        return len(data_error1), len(data_error2)





class Verify_CatState(Quantum_Operation):
    '''
    '''

    def create_and_verify_8cat_state(self):
        '''
        rep_cat_gate should be a logical gate composed with
        n preps and verifications of an 8-qubit cat state
        '''

        for ver_i in range(len(self.circuits)):
            self.stabs, self.destabs = [], []
            #local_circ = copy.deepcopy(ver_gates[ver_i].circuit_list[0])
            #local_circ = local_circ.replace_qubit_ids(range(8))
            #local_circ.update_map()
            #brow.from_circuit(local_circ, True)
            output_dict = self.run_one_circ(ver_i)
            outcomes = [outcome[0] for outcome in output_dict.values()]
            if outcomes[0]==0 and outcomes[1]==0:
                break

        # number of times we had to prep the cat state
        n_cat = ver_i + 1

        return n_cat



class QEC_d5(QEC_d3):
    '''
    QEC_d5 inherits from QEC_d3, which in turn inherits from Quantum_Operation
    This is not ideal, but I'm being careful about not over-writing methods.
    '''


    def measure_octagon(self, weight8_gate):
        '''
        measures the weight-8 stabilizer
        weight8_circ should be a gate with 2 logical gates:
        (1) the prep and verification of 8-qubit cat state
        (2) its coupling to the data qubits

        Outputs:  (1) the parity of the stabilizer measurement.
                  (2) the number of times the cat-state verf was done.
        '''
       
        prep_8cat = weight8_gate.circuit_list[0].gates[0]
        ver_circs = [g.circuit_list[0] for g in prep_8cat.circuit_list[0].gates]
        ver_oper = Verify_CatState([[],[]], ver_circs, self.chp_loc)
        n_cat = ver_oper.create_and_verify_8cat_state()
       
        # copy the circuit where we couple the cat state to the data
        couple_8cat = weight8_gate.circuit_list[0].gates[1]
        couple_circ = copy.deepcopy(couple_8cat.circuit_list[0])
        # convert the last 8 qubits to ancillary qubits (only on the copy)
        couple_circ.to_ancilla(range(17,17+8))

        self.stabs, self.destabs = qfun.combine_stabs(
                                            [self.stabs, ver_oper.stabs],
                                            [self.destabs, ver_oper.destabs])
        
        oct_dict = self.run_one_circ(couple_circ)
        oct_parity = sum([outcome[0] for outcome in oct_dict.values()])%2

        return oct_parity, n_cat
        


    def measure_stabilizers_one_kind(self, stab_gate, stab_kind):
        '''
        '''

        octagon_gate = stab_gate.circuit_list[0].gates[0]
        other_stabs_gate = stab_gate.circuit_list[0].gates[1]

        oct_par, n_cat = self.measure_octagon(octagon_gate)
        
        other_stabs_oper = Quantum_Operation([self.stabs[:], self.destabs[:]],
                                             other_stabs_gate.circuit_list[:],
                                             self.chp_loc)
        
        dat_err = other_stabs_oper.run_one_diVincenzo(0, 'd5color', 
                                                      stab_kind, oct_par)
        self.stabs = other_stabs_oper.stabs[:]
        self.destabs = other_stabs_oper.destabs[:]
            
        return dat_err


    
    def run_bare_anc_d5(self, redun=3):
        '''
        with bare anc
        '''
        
        Z_data_errors, X_data_errors = [], []
        do_X, do_Z = True, True

        if redun==3:
            for circ_ind in range(2):
                if do_X:
                    Xcirc = self.circuits[0].gates[2*circ_ind].circuit_list[0]
                    Xcirc = Xcirc.gates[0].circuit_list[0].gates[0].circuit_list[0]
                    Xcirc_oper = Quantum_Operation([self.stabs[:], self.destabs[:]],
                                                   [Xcirc], self.chp_loc)
                    out_dict = Xcirc_oper.run_one_circ(0)
                    self.stabs, self.destabs = Xcirc_oper.stabs[:], Xcirc_oper.destabs[:]
                    X_stabs = [out_dict[i][0] for i in range(17,17+8)]
                    X_stabs_dec = qfun.binary_to_decimal(X_stabs)
                    Z_err = d5color.Code.lookuptable_str[str(X_stabs_dec)]
                    Z_err = [i if i=='I' else 'Z' for i in Z_err]
                    Z_data_errors += [Z_err]
                    if Z_err.count('Z') == 0:  do_X = False

                if do_Z:
                    Zcirc = self.circuits[0].gates[2*circ_ind+1].circuit_list[0]
                    Zcirc = Zcirc.gates[0].circuit_list[0].gates[0].circuit_list[0]
                    Zcirc_oper = Quantum_Operation([self.stabs[:], self.destabs[:]],
                                                   [Zcirc], self.chp_loc)
                    out_dict = Zcirc_oper.run_one_circ(0)
                    self.stabs, self.destabs = Zcirc_oper.stabs[:], Zcirc_oper.destabs[:]
                    Z_stabs = [out_dict[i][0] for i in range(17,17+8)]
                    Z_stabs_dec = qfun.binary_to_decimal(Z_stabs)
                    X_err = d5color.Code.lookuptable_str[str(Z_stabs_dec)]
                    X_err = [i if i=='I' else 'X' for i in X_err]
                    X_data_errors += [X_err]
                    if X_err.count('X') == 0:  do_Z = False
                


            # if the outcomes of the 2 X stabs measurements
            # don't coincide, do it a third time
            if (do_X) and (Z_data_errors[0] != Z_data_errors[1]):
                Xcirc = self.circuits[0].gates[4].circuit_list[0]
                Xcirc = Xcirc.gates[0].circuit_list[0].gates[0].circuit_list[0]
                Xcirc_oper = Quantum_Operation([self.stabs[:], self.destabs[:]],
                                               [Xcirc], self.chp_loc)
                out_dict = Xcirc_oper.run_one_circ(0)
                self.stabs, self.destabs = Xcirc_oper.stabs[:], Xcirc_oper.destabs[:]
                X_stabs = [out_dict[i][0] for i in range(17,17+8)]
                X_stabs_dec = qfun.binary_to_decimal(X_stabs)
                Z_err = d5color.Code.lookuptable_str[str(X_stabs_dec)]
                Z_err = [i if i=='I' else 'Z' for i in Z_err]
                Z_data_errors += [Z_err]
                                


            # same for the Z stabs
            if (do_Z) and (X_data_errors[0] != X_data_errors[1]):
                Zcirc = self.circuits[0].gates[5].circuit_list[0]
                Zcirc = Zcirc.gates[0].circuit_list[0].gates[0].circuit_list[0]
                Zcirc_oper = Quantum_Operation([self.stabs[:], self.destabs[:]],
                                               [Zcirc], self.chp_loc)
                out_dict = Zcirc_oper.run_one_circ(0)
                self.stabs, self.destabs = Zcirc_oper.stabs[:], Zcirc_oper.destabs[:]
                Z_stabs = [out_dict[i][0] for i in range(17,17+8)]
                Z_stabs_dec = qfun.binary_to_decimal(Z_stabs)
                X_err = d5color.Code.lookuptable_str[str(Z_stabs_dec)]
                X_err = [i if i=='I' else 'X' for i in X_err]
                X_data_errors += [X_err]

            #print 'X_errors =', X_data_errors
            #print 'Z_errors =', Z_data_errors

            Z_corr, X_corr = Z_data_errors[-1], X_data_errors[-1]
        

            # update the final states only if a correction
            # is needed, to save some time
            if 'Z' in Z_corr:
                #print 'stabs =', self.stabs
                #print 'Z_corr =', Z_corr
                #Z_corr = pre_Is[:] + Z_corr[:] + post_Is[:]
                #print 'Z_corr =', Z_corr
                corr_state = qfun.update_stabs(self.stabs,
                                               self.destabs,
                                               Z_corr)
                self.stabs = corr_state[0][:]
                self.destabs = corr_state[1][:]
        
            if 'X' in X_corr:
                #X_corr = pre_Is[:] + X_corr[:] + post_Is[:]
                corr_state = qfun.update_stabs(self.stabs,
                                               self.destabs,
                                               X_corr)
                #print 'X_corr =', X_corr
                self.stabs = corr_state[0][:]
                self.destabs = corr_state[1][:]

            return len(Z_data_errors), len(X_data_errors)



    def run_fullQEC_CSS_color(self):
        '''
        runs QEC for the distance-5 4.8.8 color code.
        At the end, it applies a correction.
        It returns the number of QEC rounds.
        It assumes the X stabilizers come first.
        pre_n:  number of physical qubits before the physical
                qubits onto which this QEC is acting.
        post_n: number of physical qubits after the physical
                qubits onto which this QEC is acting.
        '''

        circ_to_run = self.circuits[0].gates[0].circuit_list[0]
        QEC_func = self.measure_stabilizers_one_kind

        #data_qs = self.circuits[0].data_qubits()
        #data_q_ids = [q.qubit_id for q in data_qs]
        #pre_ns = min(data_q_ids)
        #pre_Is = ['I' for i in range(pre_ns)]
        #post_ns = len(self.stabs) - max(data_q_ids) - 1
        #post_Is = ['I' for i in range(post_ns)]
        #print 'pre =', pre_ns
        #print 'post =', post_ns
        
        # to make things quicker; change this later
        pre_Is, post_Is = [], []

        
        Z_data_errors, X_data_errors = [], []

        list_subs = range(10)
        while len(list_subs) > 0:
            next_sub = list_subs.pop(0)
            if next_sub%2 == 0:
                residue = 0
                stab_kind = 'X'
                data_errors = Z_data_errors
            else:
                residue = 1
                stab_kind = 'Z'
                data_errors = X_data_errors
        
            data_errors += [self.measure_stabilizers_one_kind(
                                    circ_to_run.gates[next_sub],
                                    stab_kind)]

            if len(data_errors) < 3:  continue

            if (data_errors[-3]==data_errors[-2]) and (data_errors[-2]==data_errors[-1]):
                list_subs = qfun.remove_given_parity(list_subs, residue)

        Z_corr = Z_data_errors[-1]
        X_corr = X_data_errors[-1]

        # update the final states only if a correction
        # is needed, to save some time
        if 'Z' in Z_corr:
            #print 'stabs =', self.stabs
            #print 'Z_corr =', Z_corr
            Z_corr = pre_Is[:] + Z_corr[:] + post_Is[:]
            #print 'Z_corr =', Z_corr
            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           Z_corr)
            self.stabs = corr_state[0][:]
            self.destabs = corr_state[1][:]
        
        if 'X' in X_corr:
            X_corr = pre_Is[:] + X_corr[:] + post_Is[:]
            corr_state = qfun.update_stabs(self.stabs,
                                           self.destabs,
                                           X_corr)
            #print 'X_corr =', X_corr
            self.stabs = corr_state[0][:]
            self.destabs = corr_state[1][:]

        return len(Z_data_errors), len(X_data_errors)



    def run_QEC_d5(self):
        '''
        Run the whole thing
        Specifically thought for a CSS d-5 code with bare ancillae (surface-49 for now)
        The meta-decoder is essentially the same as for the d-5 color code, with the 
        exception that there are no flags.
        '''
        
        trivial_syn = [0 for i in range(12)]
        list_syndromesX = [trivial_syn[:], trivial_syn[:]]
        list_syndromesZ = [trivial_syn[:], trivial_syn[:]]
        list_sub_indices = []
        n_QEC = 0
        QEC_to_break = 8
        decided_to_break = False
        errors_det = 0
        added_error_previousX = False
        added_error_previousZ = False

        for rep in range(4):
        
            # First X stabilizers
            #print 'X stabs %i' %rep
            #print 'stabs before =', self.stabs
            index_X = 2*rep
            syndromeX = self.run_one_bare_anc_new(index_X, 'surface49', 'X') 
            #outputX = self.run_one_bare_anc_new(index_X, 'surface49', 'X',
            #                                    errors_det, list_syndromesX,
            #                                    list_syndromesZ, data_error_previous_list)
            #syndromeX, errors_det, data_error_previous = outputX
            #print 'syndrome X =', syndromeX
            #print 'errors determined =', errors_det
            
            list_syndromesX += [syndromeX]
            #list_flagsX += [flagsX]
            list_sub_indices += [index_X]
            n_QEC += 1
            #data_error_previous_list = [data_error_previous_list[1], data_error_previous]
            
            # New criteria to determine if an error has occured
            clause1 = list_syndromesZ[-2] == list_syndromesZ[-1]
            clause2 = list_syndromesX[-2] != list_syndromesX[-1]
            if clause1 and clause2 and (not added_error_previousX):
                errors_det += 1
                added_error_previousX = True
            else:
                added_error_previousX = False
                

            #print 'stabs after =', self.stabs

            #print 'Errors detected on this round =', data_error_previous
            #print 'Errors detected so far =', errors_det

            # Decide when to break
            if n_QEC == QEC_to_break:  break
            if not decided_to_break and errors_det >= 2:
                decided_to_break = True
                QEC_to_break = n_QEC + 2

            # Then Z stabilizers
            #print 'Z stabs %i' %rep
            #print 'stabs before =', self.stabs
            index_Z = 2*rep + 1
            syndromeZ = self.run_one_bare_anc_new(index_Z, 'surface49', 'Z')
            #outputZ = self.run_one_bare_anc_new(index_Z, 'surface49', 'Z',
            #                                    errors_det, list_syndromesX,
            #                                    list_syndromesZ, data_error_previous_list)
            #syndromeZ, errors_det, data_error_previous = outputZ
            #print 'syndrome Z =', syndromeZ
            #print 'errors determined =', errors_det
            list_syndromesZ += [syndromeZ]
            #list_flagsZ += [flagsZ]
            list_sub_indices += [index_Z]
            n_QEC += 1
            #data_error_previous_list = [data_error_previous_list[1], data_error_previous]
            
            # New criteria to determine if an error has occured
            clause1 = list_syndromesX[-2] == list_syndromesX[-1]
            clause2 = list_syndromesZ[-2] != list_syndromesZ[-1]
            if clause1 and clause2 and (not added_error_previousZ):
                errors_det += 1
                added_error_previousZ = True
            else:
                added_error_previousZ = False

            #print 'stabs after =', self.stabs
            
            #print 'Errors detected on this round =', data_error_previous
            #print 'Errors detected so far =', errors_det

            # Decide when to break
            if n_QEC == QEC_to_break:  break
            if not decided_to_break and errors_det >= 2:
                decided_to_break = True
                QEC_to_break = n_QEC + 2
      
            #print 'Decided to break?', decided_to_break
            #print 'QEC to break = ', QEC_to_break

            # Extra conditions to break
            if n_QEC == 2 and errors_det == 0:  break
            elif n_QEC == 6 and errors_det == 1:  break
            elif n_QEC == 8 and errors_det == 2:  break


        #print 'n_QEC = ', n_QEC
        #print 'X syndromes =', list_syndromesX
        #print 'Z syndromes =', list_syndromesZ
        #print 'X flags =', list_flagsX
        #print 'Z flags =', list_flagsZ
        #print list_sub_indices
        #sys.exit(0)

        # Combine or add the flags
        #combined_flagsX = qfun.combine_flags(list_flagsX)
        #combined_flagsZ = qfun.combine_flags(list_flagsZ)

        #print 'X flags combined =', combined_flagsX
        #print 'Z flags combined =', combined_flagsZ

        # Take last syndromes to be the "correct" ones
        last_syndromeX = ''.join(map(str,list_syndromesX[-1]))
        last_syndromeZ = ''.join(map(str,list_syndromesZ[-1]))

        # Read the correction from the lookup tables.
        surf49_lookup = surf49.Code.lookuptable
        # correction of X errors uses flagsX and syndromeZ
        corrX = surf49_lookup['Zstabs'][last_syndromeZ]
        # correction of Z errors viceversa
        corrZ = surf49_lookup['Xstabs'][last_syndromeX]

        #print 'corrX =', corrX
        #print 'corrZ =', corrZ

        # Perform the correction on the final state
        final_corrX = ['I' if oper==0 else 'X' for oper in corrX]
        final_corrZ = ['I' if oper==0 else 'Z' for oper in corrZ]
        
        #print final_corrX
        #print final_corrZ


        self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                     self.destabs[:],
                                                     final_corrX)
        self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                     self.destabs[:],
                                                     final_corrZ)

        return list_sub_indices




class QEC_with_flags(Quantum_Operation):
    '''
    '''

    def run_one_roundCSS(self, circuit, previous_flag_outcomes,
                         n_flags=[1,1,1]):
        '''
        previous_flags_outcomes:  the triggering pattern 
        n_flags:  the number of flags for each stabilizer
        QEC d3:   [1,1,1]
        QEC d5:   [2,1,1,1,1,1,1,1]
        '''

        out_dict = self.run_one_circ(circuit)
        corr, flag_outcomes = qfun.get_syn_with_flags(out_dict,
                                                      previous_flag_outcomes,
                                                      n_flags)
        
        return corr, flag_outcomes

    def run_first_round_d5(self, init_state, QEC_circs):
        '''
        Makeshift function only to test that we can correct every hook error
        '''
        corrX, flags_outcomesX = self.run_one_roundCSS(0, ((0,0),0,0,0,0,0,0,0), 
                                                       [2,1,1,1,1,1,1,1])
        print corrX, flags_outcomesX
        corrZ, flags_outcomesZ = self.run_one_roundCSS(1, flags_outcomesX,
                                                       [2,1,1,1,1,1,1,1])
        print corrZ, flags_outcomesZ                   
        
        return corrX, flags_outcomesX, corrZ, flags_outcomesZ



    def run_one_type_stabilizers_high_indet_d5(self, index_first_subcirc,
                                               list_syndromesX, list_syndromesZ,
                                               list_flagsX, list_flagsZ, 
                                               n_QEC=0, errors_det=0,
                                               added_data_err_previous=[False,False]):
        '''
        Runs one of round of 8 X or 8 Z stabilizers for the 4.8.8 d-5 color code
        with high indeterminancy.
        This means that as soon as we detect 2 errors, we stop using flags and
        start using bare ancillae.  This can happen after several events:
        (1) 2 flags get triggered (on separate stabilizers because the 2 flags
            of the octagon can get triggered by a single error)
        (2) we obtain a syndrome corresponding to 2 errors
        (3) 1 flag gets triggered and another stabilizer gets triggered as well.

        Inputs:
        (1) index_first_subcirc:  the index of the first sub-circuit to be run
                                  because the circuits are stored as a list.
        (2) list_syndromesX(Z):   the list of the previous X(Z) stab syndromes.
        (3) list_flagsX(Z):  the list of the previous X(Z) stab syndromes.
        (4) n_QEC:  the number of QEC steps (X stabs and Z stabs) run so far.
        (5) errors_det:     the number of errors detected so far. It it's larger than
                            or equal to 2 we don't use flags.
        (6) added_data_err_previous:  whether a data error was detected on the previous
                                      2 steps.
        '''

        syndromes, flags, subcircs_indices = [], [], []

        if errors_det >= 2:
            
            for i in range(8):
                i_subcirc = index_first_subcirc + 2*i + 1
                out_dict = self.run_one_circ(i_subcirc)
                syndromes += [out_dict.values()[0][0]]
                if i_subcirc%16 == 0 or i_subcirc%16 == 1:
                    flags += [(0,0)]
                else:  
                    flags += [0]
                subcircs_indices += [i_subcirc]

        else:

            #stab_light = False
            for i in range(8):
               
                #print errors_det

                if errors_det >= 2:
                    i_subcirc = index_first_subcirc + 2*i + 1
                    out_dict = self.run_one_circ(i_subcirc)
                    syn = out_dict.values()[0][0]
                    if i_subcirc%16 == 0 or i_subcirc%16 == 1:
                        flag = (0,0)
                    else:
                        flag = 0
                    flags += [flag]
    
                else:
                    i_subcirc = index_first_subcirc + 2*i
                    out_dict = self.run_one_circ(i_subcirc)
                    out_keys = out_dict.keys()[:]
                    out_keys.sort()
                    syn = out_dict[out_keys[0]][0]
                    flag = [out_dict[out_keys[i]][0] for i in range(1,len(out_keys))]
                    if len(out_keys) > 2:
                        flags += [tuple(flag)]
                    else:
                        flags += [flag[0]]
                    
                    #if syn == 1:
                    #    if not stab_light:
                    #        stab_light = True
                    #        errors_det += 1
                    #else:
                    #    if sum(flag) > 0:
                    #        errors_det += 1
                    if sum(flag) > 0:
                        errors_det += 1

                syndromes += [syn]
                subcircs_indices += [i_subcirc]
       
        #print 'Syndromes =', syndromes
        #print 'Flags =', flags

        # If a new data errror was detected by the stabilizers, add 1 to errors_det
        # We have to be very careful (conservative).  We only add 1 if the syndrome
        # is different from the last 2 syndromes, if no flags were triggered in the 2
        # previous steps and if no data error was detected on the previous step.
        # I'm sure we can refine it even further.
        if len(list_syndromesX) == 0:
            last_syndromeX = [0 for i in range(8)]
        else:
            last_syndromeX = list_syndromesX[-1][:]
        if len(list_syndromesZ) == 0:
            last_syndromeZ = [0 for i in range(8)]
        else:
            last_syndromeZ = list_syndromesZ[-1][:]
        
        if len(list_flagsX) == 0:
            last_flagX = [(0,0)] + [0 for i in range(7)]
        else:
            last_flagX = list_flagsX[-1][:]
        if len(list_flagsZ) == 0:
            last_flagZ = [(0,0)] + [0 for i in range(7)]
        else:
            last_flagZ = list_flagsZ[-1][:]

        sum_flagX = sum(last_flagX[0]) + sum(last_flagX[1:])
        sum_flagZ = sum(last_flagZ[0]) + sum(last_flagZ[1:])
        sum_current_flag = sum(flags[0]) + sum(flags[1:])
        
        # Define logical clauses
        data_error_det = False
        syn_clause = (syndromes != last_syndromeX) and (syndromes != last_syndromeZ)
        flag_clause = (sum_flagX == 0) and (sum_flagZ == 0) and (sum_current_flag == 0)
        data_err_clause = (not added_data_err_previous[0]) and (not added_data_err_previous[1])
        if (syn_clause) and (flag_clause) and (data_err_clause):
            errors_det += 1
            data_error_det = True


        return syndromes, flags, subcircs_indices, errors_det, data_error_det


    
    def run_QEC_d5(self):
        '''
        Run the whole thing
        '''
        
        list_syndromesX, list_syndromesZ = [], []
        list_flagsX, list_flagsZ = [], []
        list_sub_indices = []
        n_QEC = 0
        QEC_to_break = 8
        decided_to_break = False
        errors_det = 0
        data_error_previous_list = [False,False]

        for rep in range(4):
        
            # First X stabilizers
            #print 'X stabs %i' %rep
            #print 'stabs before =', self.stabs
            index_firstX = 32*rep
            outputX = self.run_one_type_stabilizers_high_indet_d5(index_firstX,
                                                                  list_syndromesX,
                                                                  list_syndromesZ,
                                                                  list_flagsX,
                                                                  list_flagsZ,
                                                                  n_QEC,
                                                                  errors_det,
                                                                  data_error_previous_list)
            syndromeX, flagsX, sub_indices, errors_det, data_error_previous = outputX
            list_syndromesX += [syndromeX]
            list_flagsX += [flagsX]
            list_sub_indices += sub_indices
            n_QEC += 1
            data_error_previous_list = [data_error_previous_list[1], data_error_previous]

            #print 'stabs after =', self.stabs

            #print 'Errors detected so far =', errors_det

            # Decide when to break
            if n_QEC == QEC_to_break:  break
            if not decided_to_break and errors_det >= 2:
                decided_to_break = True
                QEC_to_break = n_QEC + 2


            # Then Z stabilizers
            #print 'Z stabs %i' %rep
            #print 'stabs before =', self.stabs
            index_firstZ = 32*rep + 16
            outputZ = self.run_one_type_stabilizers_high_indet_d5(index_firstZ,
                                                                  list_syndromesX,
                                                                  list_syndromesZ,
                                                                  list_flagsX,
                                                                  list_flagsZ,
                                                                  n_QEC,
                                                                  errors_det,
                                                                  data_error_previous_list)
            syndromeZ, flagsZ, sub_indices, errors_det, data_error_previous = outputZ
            list_syndromesZ += [syndromeZ]
            list_flagsZ += [flagsZ]
            list_sub_indices += sub_indices
            n_QEC += 1
            data_error_previous_list = [data_error_previous_list[1], data_error_previous]
            
            #print 'stabs after =', self.stabs
            
            #print 'Errors detected so far =', errors_det

            # Decide when to break
            if n_QEC == QEC_to_break:  break
            if not decided_to_break and errors_det >= 2:
                decided_to_break = True
                QEC_to_break = n_QEC + 2
            
            
            # Extra conditions to break
            if n_QEC == 2 and errors_det == 0:  break
            elif n_QEC == 6 and errors_det == 1:  break
            elif n_QEC == 8 and errors_det == 2:  break


        #print 'n_QEC = ', n_QEC
        #print 'X syndromes =', list_syndromesX
        #print 'Z syndromes =', list_syndromesZ
        #print 'X flags =', list_flagsX
        #print 'Z flags =', list_flagsZ
        #print list_sub_indices
        #sys.exit(0)

        # Combine or add the flags
        combined_flagsX = qfun.combine_flags(list_flagsX)
        combined_flagsZ = qfun.combine_flags(list_flagsZ)

        #print 'X flags combined =', combined_flagsX
        #print 'Z flags combined =', combined_flagsZ

        # Take last syndromes to be the "correct" ones
        last_syndromeX = ''.join(map(str,list_syndromesX[-1]))
        last_syndromeZ = ''.join(map(str,list_syndromesZ[-1]))

        # Read the correction from the lookup tables.
        d5_lookups = d5color.Code.all_lookups
        # correction of X errors uses flagsX and syndromeZ
        corrX = d5_lookups[combined_flagsX][last_syndromeZ]
        # correction of Z errors viceversa
        corrZ = d5_lookups[combined_flagsZ][last_syndromeX]

        #print 'corrX =', corrX
        #print 'corrZ =', corrZ

        # Perform the correction on the final state
        final_corrX = ['I' if oper==0 else 'X' for oper in corrX]
        final_corrZ = ['I' if oper==0 else 'Z' for oper in corrZ]

        self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                     self.destabs[:],
                                                     final_corrX)
        self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                     self.destabs[:],
                                                     final_corrZ)

        return list_sub_indices



    def run_one_type_stabilizers_high_indet(self, index_first_subcirc,
                                            previous_flags=(0,0,0),
                                            error_det=False,
                                            change_error_det='flag'):
        '''
        Runs one round of 3 X or 3 Z stabilizers with high indeterminancy.
        This means that was soon as one flag is triggered or one error is
        detected, we stop using flags and start using bare ancillae.
        Inputs:
        (1) index_first_subcirc:  the index of the first sub-circuit to be run
                                  because the circuits are stored as a list.
        (2) previous_flag:  the outcome of the previous round.  If no flag was
                            previously triggered, then (0,0,0).  Else, the index
                            of the stabilizer whose flag was triggered: (0,1,0)
                            if the second stabilizer was triggered.
        (3) error_det:   whether an error has been detected previously.  If so,
                         we don't use any flags.
        '''
        
        syn, flags, subcircs_indices = [], [], []
        for i in range(3):
            #print 'Error det:', error_det
            if error_det:
                i_subcirc = index_first_subcirc + 2*i + 1
                out_dict = self.run_one_circ(i_subcirc)
                syn += [out_dict.values()[0][0]]
                flags += [0]
            else:
                i_subcirc = index_first_subcirc + 2*i 
                out_dict = self.run_one_circ(i_subcirc)
                out_keys = out_dict.keys()[:]
                out_keys.sort()
                syn += [out_dict[out_keys[0]][0]]
                flag = out_dict[out_keys[1]][0]
                if flag == 1:
                    error_det = True
                elif change_error_det == 'any':
                    if syn[-1] == 1:
                        error_det = True
                flags += [flag]
            
            subcircs_indices += [i_subcirc]

        corr = st.Code.total_lookup_table_one_flag[previous_flags][tuple(syn)]

        return corr, flags, error_det, subcircs_indices



    def run_stabilizers_high_indet_ion(self, index_first_subcirc, alternating=True,
                                       change_error_det='any', decoder='old'):
        '''
        Runs 1 or 2 rounds of the d-3 color code stabilizers.
        The first round has both the FT flagged stabilizer measurements and the
        non-FT unflagged measurement.  The second round has only the non-FT
        version because it is only run if an error occurred during the first round.
        
        Inputs:
        (1) index_first_subcirc:  always 0 for now.
        (2) alternaring: refers to whether the stabilizers alternate between
        X and Z: Sx1, Sz1, Sx2, Sz2, Sx3, Sz3; or not.
        
        The alternating version requires less shuttling.
        '''
        
        # to avoid adding an extra X on the syndrome and the flag qubits, 
        # we just re-interpret the outcome: in this case '1' means no error
        # and '0' means error.
        flag_no_trig = 1
        syn_no_error = 1
        syndromes, flags = [], []
        subcircs_indices = []
        
        error_det = False

        # New decoder
        if decoder == 'new':
            
            # First run the FT subcircuits
            for i in range(9):
                i_subcirc = index_first_subcirc + i
                out_dict = self.run_one_circ(i_subcirc)
                
                # The subcircs 2, 5, 8 are shuttling.
                if i_subcirc%3 != 2:
                    out_keys = out_dict.keys()[:]
                    out_keys.sort()
                    syndromes += [out_dict[out_keys[0]][0]]
                    flag = out_dict[out_keys[1]][0]
                    if flag != flag_no_trig:
                        error_det = True
                    elif change_error_det == 'any':
                        if syndromes[-1] != syn_no_error:
                            error_det = True
                    flags += [flag]
                
                subcircs_indices += [i_subcirc]
                
                # as soon as an error is detected, we break.
                if error_det:  break


            # If no error was detected during the first round
            # we stop and don't correct.
            if not error_det:  return subcircs_indices
            
            # Add trivial flags to the flags list
            if len(subcircs_indices) < 3:
                subcircs_run = len(subcircs_indices)
                case = 1
                next_subcircs = range(9,18)
            elif len(subcircs_indices) < 6:
                subcircs_run = len(subcircs_indices) - 1
                next_subcircs = range(12,18) + [9,10]
                case = 2
            elif len(subcircs_indices) < 9:
                subcircs_run = len(subcircs_indices) - 2
                next_subcircs = [15,16,14,12,13,11,9,10]
                case = 3
            n_extra_flags = 6 - subcircs_run
            extra_flags = [flag_no_trig for i in range(n_extra_flags)]
            flags += extra_flags

            # Now we run the second round (only non-FT stabs)
            syndromes = []
            for i_subcirc in next_subcircs:
                out_dict = self.run_one_circ(i_subcirc)
                
                # The subcircs 11, 14, 17 are shuttling.
                if i_subcirc%3 != 2:
                    syndromes += [out_dict.values()[0][0]]
                
                subcircs_indices += [i_subcirc]

            if case == 2:
                syndromes = [syndromes[4], syndromes[5], syndromes[0],
                             syndromes[1], syndromes[2], syndromes[3]]
            if case == 3:
                syndromes = [syndromes[4], syndromes[5], syndromes[2],
                             syndromes[3], syndromes[0], syndromes[1]]



        # Old decoder
        else:

            for i in range(6):

                # if an error has been determined, run the odd subcircuit
                if error_det:
                    i_subcirc = index_first_subcirc + 2*i + 1
                    out_dict = self.run_one_circ(i_subcirc)
                    syndromes += [out_dict.values()[0][0]]
                    flags += [flag_no_trig]
                else:
                    i_subcirc = index_first_subcirc + 2*i
                    out_dict = self.run_one_circ(i_subcirc)
                    out_keys = out_dict.keys()[:]
                    out_keys.sort()
                    syndromes += [out_dict[out_keys[0]][0]]
                    flag = out_dict[out_keys[1]][0] 
                    if flag != flag_no_trig:
                        error_det = True
                    elif change_error_det == 'any':
                        if syndromes[-1] != syn_no_error:
                            error_det = True
                    flags += [flag]

                subcircs_indices += [i_subcirc]

            
            # if an error happened on the first round, we need a second round
            if not error_det:  return subcircs_indices
            else:
                syndromes = []
                for i in range(6):
                    i_subcirc = 12 + i
                    out_dict = self.run_one_circ(i_subcirc)
                    syndromes += [out_dict.values()[0][0]]
                    subcircs_indices += [i_subcirc]

                #print 'syndromes extra =', syndromes

            #print 'syndromes =', syndromes
            #print 'flags =', flags
            #print 'subcircs indices =', subcircs_indices
        


        if alternating:
            X_is, Z_is = [0,2,4], [1,3,5]
        else:
            X_is, Z_is = [0,1,2], [3,4,5]
        
        #print syndromes
        #print flags
        syndromesX = tuple([(syndromes[i]+syn_no_error)%2 for i in X_is])
        flagsX = tuple([(flags[i]+flag_no_trig)%2 for i in X_is])
        syndromesZ = tuple([(syndromes[i]+syn_no_error)%2 for i in Z_is])
        flagsZ = tuple([(flags[i]+flag_no_trig)%2 for i in Z_is])

        #print 'synX =', syndromesX
        #print 'synZ =', syndromesZ

        # Perform the correction on the final state
        corrX = st.Code.total_lookup_table_one_flag[flagsX][syndromesZ]
        corrZ = st.Code.total_lookup_table_one_flag[flagsZ][syndromesX]
       
        corrX = [oper if oper=='I' else 'X' for oper in corrX]
        corrZ = [oper if oper=='I' else 'Z' for oper in corrZ]

        self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                     self.destabs[:],
                                                     corrX)
        self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                     self.destabs[:],
                                                     corrZ)

        return subcircs_indices



    def run_Reichardt_d3_one_flag_stab(self, metadecoder='cheap',
                                       change_error_det='flag'):
        '''
        metadecoder:  'cheap': 2.2 from Notes
                      'standard':  more standard form
        '''
        
        initial_flags = (0,0,0)
        error_det = False
        
        # We start with X stabilizers
        output_run = self.run_one_type_stabilizers_high_indet(0,
                                                         initial_flags,
                                                         error_det,
                                                         change_error_det)
        corrX, flagsX, error_det, subcircs_i = output_run       
        
        #print corrX
        #print flagsX
        
        # If no errors were detected and no flags were triggered 
        if corrX.count('E') == 0 and flagsX.count(1) == 0:
        
            # Now we measure the Z stabilizers
            output_run = self.run_one_type_stabilizers_high_indet(6,
                                                            tuple(flagsX),
                                                            error_det,
                                                            change_error_det)
            corrZ, flagsZ, error_det, subcircs_i1 = output_run
            subcircs_i += subcircs_i1[:]

            # If no errors were detected and no flags were triggered
            if corrZ.count('E') == 0 and flagsZ.count(1) == 0:
                final_corrX, final_corrZ = corrX[:], corrZ[:]

            # Error detected but no flag triggered:
            # re-measure only Z stabilizers
            elif corrZ.count('E') > 0 and flagsZ.count(1) == 0:
                output_run = self.run_one_type_stabilizers_high_indet(18,
                                                                tuple(flagsX),
                                                                error_det,
                                                                change_error_det)
                corrZ2, flagsZ2, error_det, subcircs_i3 = output_run
                subcircs_i += subcircs_i3[:]
                final_corrX, final_corrZ = corrX[:], corrZ2[:]

            # No error detected but flag triggered:
            # re-measure the X stabilizers to identify the Z hook error
            elif corrZ.count('E') == 0 and flagsZ.count(1) > 0:
                output_run = self.run_one_type_stabilizers_high_indet(12,
                                                                tuple(flagsZ),
                                                                error_det,
                                                                change_error_det)
                corrX2, flagsX2, error_det, subcircs_i2 = output_run
                subcircs_i += subcircs_i2[:]
                final_corrX, final_corrZ = corrX2[:], corrZ[:]

            # Error detected and flag triggered:
            # re-measure the X stabilizers to correctly identify Z error
            # re-measure the Z stabilizers to identify the X error.
            else:
                output_run = self.run_one_type_stabilizers_high_indet(12,
                                                                tuple(flagsZ),
                                                                error_det,
                                                                change_error_det)
                corrX2, flagsX2, error_det, subcircs_i2 = output_run
                subcircs_i += subcircs_i2[:]
                output_run = self.run_one_type_stabilizers_high_indet(18,
                                                                   tuple(flagsX),
                                                                   error_det,
                                                                   change_error_det)
                corrZ2, flagsZ2, error_det, subcircs_i3 = output_run
                subcircs_i += subcircs_i3[:]
                final_corrX, final_corrZ = corrX2[:], corrZ2[:]


        # Error detected but no flag triggered:
        elif corrX.count('E') > 0 and flagsX.count(1) == 0:
            
            # Measure the Z stabilizers:
            output_run = self.run_one_type_stabilizers_high_indet(6,
                                                               tuple(flagsX),
                                                               error_det,
                                                               change_error_det)
            corrZ, flagsZ, error_det, subcircs_i1 = output_run
            subcircs_i += subcircs_i1[:]
            # If the Z syndrome indicates an error, apply a Y correction based on
            # this syndrome (not the X syndrome)
            if corrZ.count(1) > 0:
                final_corrX, final_corrZ = corrZ[:], corrZ[:]
            # If not, then re-measure the X stabilizers to correctly identify the Z
            # error
            else:
                output_run = self.run_one_type_stabilizers_high_indet(12,
                                                                tuple(flagsZ),
                                                                error_det,
                                                                change_error_det)
                corrX2, flagsX2, error_det, subcircs_i2 = output_run
                subcircs_i += subcircs_i2[:] 
                final_corrX, final_corrZ = corrX2[:], corrZ[:]
                

        # No error detected but flag triggered:
        # measure the Z stabilizers to identify the X hook error 
        elif corrX.count('E') == 0 and flagsX.count(1) > 0:
            output_run = self.run_one_type_stabilizers_high_indet(6,
                                                               tuple(flagsX),
                                                               error_det,
                                                               change_error_det)
            corrZ, flagsZ, error_det, subcircs_i1 = output_run
            subcircs_i += subcircs_i1
            final_corrX, final_corrZ = corrX[:], corrZ[:] 


        # Error detected and flag triggered
        # measure the Z stabilizers to identify the X hook error
        # measure the X stabilizers to correctly identify the Z error
        else:
            output_run = self.run_one_type_stabilizers_high_indet(6,
                                                               tuple(flagsX),
                                                               error_det,
                                                               change_error_det)
            corrZ, flagsZ, error_det, subcircs_i1 = output_run
            subcircs_i += subcircs_i1
            output_run = self.run_one_type_stabilizers_high_indet(12,
                                                                tuple(flagsZ),
                                                                error_det,
                                                                change_error_det)
            corrX2, flagsX2, error_det, subcircs_i2 = output_run
            subcircs_i += subcircs_i2
            final_corrX, final_corrZ = corrX2[:], corrZ[:]
            
      
        # Perform the correction on the final state
        final_corrX = [oper if oper=='I' else 'Z' for oper in final_corrX]
        final_corrZ = [oper if oper=='I' else 'X' for oper in final_corrZ]

        self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                     self.destabs[:],
                                                     final_corrX)
        self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                     self.destabs[:],
                                                     final_corrZ)

        return subcircs_i



    def run_one_round_Reichardt_d3(self, circuit, previous_flag_outcome,
                                   with_flag=True):
        '''
        One round of the circuit from Reichardt paper (Figure 8b) to measure
        the 3 stabilizers of the d3 color code with only one flag.
        '''

        out_dict = self.run_one_circ(circuit)
        out_keys = out_dict.keys()[:]
        out_keys.sort()
        if with_flag:
            syn = [out_dict[i][0] for i in out_keys[:-1]]
            flag = out_dict[out_keys[-1]][0]
        else:
            syn = [out_dict[i][0] for i in out_keys[:]]
            flag = None
            
        #print 'syn =', syn
        #print 'flag =', flag
        Steane_lookup = st.Code.total_lookup_table
        corr = Steane_lookup[previous_flag_outcome][tuple(syn)]
        #print 'corr = ', corr

        return corr, flag

   

    def run_all_Reichardt_d3(self, init_state):
        '''
        Makeshift function to make sure we can correct all the hook errors
        '''
        
        corrX, flagX = self.run_one_round_Reichardt_d3(0, 0)
        corrZ, flagZ = self.run_one_round_Reichardt_d3(1, flagX)
        n_errors = corrX.count('E') + corrZ.count('E') + flagX + flagZ
        
        if n_errors > 0:
            # flagX2 will always be None
            corrX2, flagX2 = self.run_one_round_Reichardt_d3(2, flagZ, False)
            # flagZ2 will always be None as well.
            corrZ2, flagZ2 = self.run_one_round_Reichardt_d3(3, flagX, False)
            
            corrX2 = [oper if oper=='I' else 'Z' for oper in corrX2]
            corrZ2 = [oper if oper=='I' else 'X' for oper in corrZ2]

            self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                         self.destabs[:],
                                                         corrX2)
            self.stabs, self.destabs = qfun.update_stabs(self.stabs[:],
                                                         self.destabs[:],
                                                         corrZ2)
            #print corrX, flagX, corrZ, flagZ, corrX2, corrZ2

        return None
            
    
    
    def run_QEC_FT_lattsurg_ion(self, stab='all'):
        '''
        output:  syndrome, flags, type of correction, error detected
        if stab=='all', we measure all the stabilizers the usual way
        if stab==2, we only measure the last stabilizer nonFT
        if stab==1, we only measure the second stabilizer nonFT
        '''

        corr_type = 'normal'
        error_det = True
        subcircs_run = []

        # All the FT circuits have a number of MS, n, gates such that
        # n(mod4) = 2, that is 6 for the s qubit and 2 for the f qubit.
        # Therefore we have to flip the outcomes.
        # We don't have to flip the outcomes for the nonFT one, which
        # has 4 MS gates.
    
        if stab==2:
            syn3 = self.run_one_circ(7).values()[0][0]
            subcircs_run += [7] 
            # reorder from stab 2 to 0
            self.run_one_circ(9)
            subcircs_run += [9] 
            return [0,0,syn3], [0,0,0], corr_type, error_det, subcircs_run
            
        elif stab==1:
            syn2 = self.run_one_circ(6).values()[0][0]
            subcircs_run += [6] 
            # reorder from stab 1 to 0
            self.run_one_circ(8)
            subcircs_run += [8] 
            return [0,syn2,0], [0,0,0], corr_type, error_det, subcircs_run


        # first we measure the stabilizer that doesn't "touch" the boundary
        output_first = self.run_one_circ(0).values()
        subcircs_run += [0] 
        syn1, flag1 = (output_first[0][0]+1)%2, (output_first[1][0]+1)%2
   
        #print 'syn1 =', syn1
        #print 'flag1 =', flag1

        # if an error is detected, we're done
        if syn1==1 or flag1==1:
            total_syn, total_flag = [syn1,0,0], [flag1,0,0]
        
        else:
            # measure the second stabilizer FT
            output_second = self.run_one_circ(1).values()
            subcircs_run += [1] 
            syn2, flag2 = (output_second[0][0]+1)%2, (output_second[1][0]+1)%2
           
            #print 'syn2 =', syn2
            #print 'flag2 =', flag2
            
            # if the flag was triggered, we're done
            if flag2==1:
                # reorder from stab 1 to 0
                self.run_one_circ(3)
                subcircs_run += [3] 
                total_syn, total_flag = [syn1,syn2,0], [flag1,flag2,0]
            
            else:
                # if an error was detected, measure the first stab nonFT to
                # make distinguish between (1) an error that happened on one of
                # qubits of the first and second stabs after the first stab was
                # measured and (2) an error that happened on one of the qubits
                # of the second stab not shared with the first stab.
                if syn2==1:
                
                    # Reorder from 1 to 0
                    self.run_one_circ(3)
                    subcircs_run += [3] 
                    syn1 = self.run_one_circ(5).values()[0][0]
                    subcircs_run += [5] 
                    #print 'syn1 (second time) =', syn1
                    
                    # if syn1 is 1, then it was the first option and we're done
                    if syn1==1:
                        total_syn, total_flag = [syn1,syn2,0], [flag1,flag2,0]
                    
                    # else: measure the second stab again, this time nonFT
                    else:
                        syn2 = self.run_one_circ(6).values()[0][0]
                        subcircs_run += [6] 
                        #print 'syn2 (second time) =', syn2
                    
                        # if the syndrome is 0 this time, that means it was a
                        # measurement error and we're done.
                        if syn2==0:
                            # reorder from 1 to 0
                            self.run_one_circ(8)
                            subcircs_run += [8] 
                            total_syn, total_flag = [syn1,syn2,0], [flag1,flag2,0]

                        # if the syndrome was 1 the second time, it was a data error
                        # and we need to measure the third stabilizer nonFT.
                        else:
                            syn3 = self.run_one_circ(7).values()[0][0]
                            subcircs_run += [7] 
                            #print 'syn3 =', syn3
                            # reorder from 2 to 0
                            self.run_one_circ(9)
                            subcircs_run += [9] 
                            total_syn, total_flag = [syn1,syn2,syn3], [flag1,flag2,0]
                            corr_type = 'unknown'

                # if no error was detected, measure the third stabilizer FT
                else:
                    output_third = self.run_one_circ(2).values()
                    subcircs_run += [2] 
                    syn3, flag3 = (output_third[0][0]+1)%2, (output_third[1][0]+1)%2
                    total_syn, total_flag = [syn1,syn2,syn3], [flag1,flag2,flag3]
                    # reorder from 2 to 0
                    self.run_one_circ(4)
                    subcircs_run += [4] 

                    # if the third flag wasn't triggered
                    if flag3==0:
                        
                        # if the flag is 0, but the third syndrome is 1, we need
                        # to measure the stabilizers again
                        if syn3==1:
                            syn1 = self.run_one_circ(5).values()[0][0]
                            subcircs_run += [5] 
                            
                            # if syn1 is 1, the error was on a non-boundary data qubit
                            # after the first stabilizer
                            if syn1==1:
                                total_syn = [syn1,syn2,syn3]
                                total_flag = [flag1,flag2,flag3]
                            else:
                                syn2 = self.run_one_circ(6).values()[0][0]
                                subcircs_run += [6] 
                                # if syn2 is 1, the error was on a boundary qubit
                                if syn2==1:
                                    # reorder from 1 to 0
                                    self.run_one_circ(8)
                                    subcircs_run += [8] 
                                    total_syn = [syn1,syn2,syn3]
                                    total_flag = [flag1,flag2,flag3]
                                    corr_type = 'unknown'
                                else:
                                    syn3 = self.run_one_circ(7).values()[0][0]
                                    subcircs_run += [7] 
                                    # reorder from 2 to 0
                                    self.run_one_circ(9)
                                    subcircs_run += [9] 
                                    total_syn = [syn1,syn2,syn3]
                                    total_flag = [flag1,flag2,flag3]
                                    # if the syndrome is 1 again, it was a data error
                                    if syn3==1:
                                        corr_type = 'unknown'
                                        

                        else:
                            error_det = False

        return total_syn, total_flag, corr_type, error_det, subcircs_run
    
   

    def run_QEC_FT_lattsurg(self, stab='all'):
        '''
        output:  syndrome, flags, type of correction, error detected
        if stab=='all', we measure all the stabilizers the usual way
        if stab==2, we only measure the last stabilizer nonFT
        if stab==1, we only measure the second stabilizer nonFT
        '''

        corr_type = 'normal'
        error_det = True

        
        if stab==2:
            syn3 = self.run_one_circ(5).values()[0][0]
            return [0,0,syn3], [0,0,0], corr_type, error_det
            
        elif stab==1:
            syn2 = self.run_one_circ(4).values()[0][0]
            return [0,syn2,0], [0,0,0], corr_type, error_det


        # first we measure the stabilizer that doesn't "touch" the boundary
        output_first = self.run_one_circ(0).values()
        syn1, flag1 = output_first[0][0], output_first[1][0]
   
        print 'syn1 =', syn1
        print 'flag1 =', flag1

        # if an error is detected, we're done
        if syn1==1 or flag1==1:
            total_syn, total_flag = [syn1,0,0], [flag1,0,0]
        
        else:
            # measure the second stabilizer FT
            output_second = self.run_one_circ(1).values()
            syn2, flag2 = output_second[0][0], output_second[1][0]
           
            print 'syn2 =', syn2
            print 'flag2 =', flag2
            
            # if the flag was triggered, we're done
            if flag2==1:
                total_syn, total_flag = [syn1,syn2,0], [flag1,flag2,0]
            
            else:
                # if an error was detected, measure the first stab nonFT to
                # make distinguish between (1) an error that happened on one of
                # qubits of the first and second stabs after the first stab was
                # measured and (2) an error that happened on one of the qubits
                # of the second stab not shared with the first stab.
                if syn2==1:
                
                    syn1 = self.run_one_circ(3).values()[0][0]
                    #print 'syn1 (second time) =', syn1
                    
                    # if syn1 is 1, then it was the first option and we're done
                    if syn1==1:
                        total_syn, total_flag = [syn1,syn2,0], [flag1,flag2,0]
                    
                    # else: measure the second stab again, this time nonFT
                    else:
                        syn2 = self.run_one_circ(4).values()[0][0]
                        #print 'syn2 (second time) =', syn2
                    
                        # if the syndrome is 0 this time, that means it was a
                        # measurement error and we're done.
                        if syn2==0:
                            total_syn, total_flag = [syn1,syn2,0], [flag1,flag2,0]

                        # if the syndrome was 1 the second time, it was a data error
                        # and we need to measure the third stabilizer nonFT.
                        else:
                            syn3 = self.run_one_circ(5).values()[0][0]
                            #print 'syn3 =', syn3
                            total_syn, total_flag = [syn1,syn2,syn3], [flag1,flag2,0]
                            corr_type = 'unknown'

                # if no error was detected, measure the third stabilizer FT
                else:
                    output_third = self.run_one_circ(2).values()
                    syn3, flag3 = output_third[0][0], output_third[1][0]
                    total_syn, total_flag = [syn1,syn2,syn3], [flag1,flag2,flag3]

                    # if the third flag wasn't triggered
                    if flag3==0:
                        
                        # if the flag is 0, but the third syndrome is 1, we need
                        # to measure the stabilizers again
                        if syn3==1:
                            syn1 = self.run_one_circ(3).values()[0][0]
                            
                            # if syn1 is 1, the error was on a non-boundary data qubit
                            # after the first stabilizer
                            if syn1==1:
                                total_syn = [syn1,syn2,syn3]
                                total_flag = [flag1,flag2,flag3]
                            else:
                                syn2 = self.run_one_circ(4).values()[0][0]
                                # if syn2 is 1, the error was on a boundary qubit
                                if syn2==1:
                                    total_syn = [syn1,syn2,syn3]
                                    total_flag = [flag1,flag2,flag3]
                                    corr_type = 'unknown'
                                else:
                                    syn3 = self.run_one_circ(5).values()[0][0]
                                    total_syn = [syn1,syn2,syn3]
                                    total_flag = [flag1,flag2,flag3]
                                    # if the syndrome is 1 again, it was a data error
                                    if syn3==1:
                                        corr_type = 'unknown'
                                        

                        else:
                            error_det = False

        return total_syn, total_flag, corr_type, error_det



    def run_QECX_nonFT(self, previous_syn):
        '''
        We only measure the 2 stabilizers on the boundary.
        If the syndrome doesn't coincide with the first
        one, we know the error occurred after the XX measurement
        and the correction is normal.
        '''

        # measure the first stabilizer
        syn2 = self.run_one_circ(0).values()[0][0]

        # if this does not coincide with the previous one, we're done
        if syn2 != previous_syn[1]:
            return 'normal'

        else:
            syn3 = self.run_one_circ(1).values()[0][0]
            if syn3 != previous_syn[2]:
                return 'normal'
            else:
                return 'alternative'



    def run_jointQECZ(self, error_det, inflags1, inflags2):
        '''
        error_det:  True if an error has already been detected
                    False otherwise.
                    If error_det is True, then we just run the
                    stabilizers 1 time nonFT.
        flags1:  the previous flags for the first logical qubit
        flags2:  the previous flags for the second logical qubit
        '''

        error_index = 1
        Pauli_error = 'X'

        #print error_det, inflags1, inflags2
        
        # undefined stab is the index of the stabilizer that is undefined.
        # With the current numbering, it's the second stab (stab 1)
        undefined_stab = 1

        outflags1, outflags2 = [], []
        if error_det:
            # In this case, we just run the non-FT circuits
        
            # Because the order of the stabilizer measurements in the QECx during
            # the merging step is different from the standard order, we need to
            # change the flag orderings
            inflags1 = [inflags1[1], inflags1[2], inflags1[0]]
            inflags2 = [inflags2[1], inflags2[2], inflags2[0]]

            # First logical qubit
            syn1 = []
            for i in range(3,6):
                syn1 += [self.run_one_circ(i).values()[0][0]]
            
            # Second logical qubit
            syn2 = []
            for i in range(12,15):
                syn2 += [self.run_one_circ(i).values()[0][0]]
            
            outflags1, outflags2 = [0,0,0], [0,0,0]


        else:
            # In this case, both inflags are [0,0,0]

            ##########################################################################
            # Alternative:  only measure the undefined stabilizers
            error_det1 = False
            syn1 = []
            out_run1 = self.run_one_circ(1).values()
            syn1 += [out_run1[0][0]]
            outflags1 = out_run1[1][0]
            if outflags1==1:  error_det1 = True

            #brow.from_circuit(self.circuits[1], True)

            #print 'syn1 =', syn1
            #print 'outflags1 =', outflags1

            error_det2 = False
            syn2 = []
            out_run2 = self.run_one_circ(10).values()
            syn2 += [out_run2[0][0]]
            outflags2 = out_run2[1][0]
            if outflags2==1:  error_det2 = True
            
            #brow.from_circuit(self.circuits[7], True)
            #sys.exit(0)

            #print 'syn2 =', syn2
            #print 'outflags2 =', outflags2
            
            # If a flag was triggered in either one of the stabs, error_det = True
            error_det = error_det1 or error_det2
            
            # If both syndromes were 1, correct both stabilizers
            if syn1[0]==1 and syn2[0]==1:
                total_corr1 = ['X' if i==1 else 'I' for i in range(7)]
                total_corr2 = ['X' if i==1 else 'I' for i in range(7)]
    
                total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                self.stabs = corr_state[0][:]
                self.destabs = corr_state[1][:]

            # If the 2 syndromes are different, error_det = True and we measure
            # the other 2 stabilizer pairs with no flags.
            elif syn1[0] != syn2[0]:
                error_det = True
                syn1 += [self.run_one_circ(3).values()[0][0]]
                #syn1 += [self.run_one_circ(5).values()[0][0]]
                syn2 += [self.run_one_circ(12).values()[0][0]]
                #syn2 += [self.run_one_circ(11).values()[0][0]]

                #print 'syn1 =', syn1
                #print 'syn2 =', syn2

                #brow.from_circuit(self.circuits[3], True)
                #brow.from_circuit(self.circuits[9], True)

                # If an error has been detected on both logical qubits, we
                # apply a correction on both qubits 1.
                if (1 in syn1) and (1 in syn2):
                    total_corr1 = ['X' if i==1 else 'I' for i in range(7)]
                    total_corr2 = ['X' if i==1 else 'I' for i in range(7)]
    
                    total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                    self.stabs = corr_state[0][:]
                    self.destabs = corr_state[1][:]
                
                else:
                    syn1 += [self.run_one_circ(5).values()[0][0]]
                    syn2 += [self.run_one_circ(14).values()[0][0]]
                   
                    #print 'syn1 =', syn1
                    #print 'syn2 =', syn2

                    if (1 in syn1) and (1 in syn2):
                        total_corr1 = ['X' if i==1 else 'I' for i in range(7)]
                        total_corr2 = ['X' if i==1 else 'I' for i in range(7)]
    
                        total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                        self.stabs = corr_state[0][:]
                        self.destabs = corr_state[1][:]


                #if (syn1[0]==0 and syn1[1]==0) or (syn2[0]==0 and syn2[1] == 0):
                #    return error_det, [0,outflags1,0], [0,outflags2,0]

                #else:
                #    total_corr1 = ['X' if i==1 else 'I' for i in range(7)]
                #    total_corr2 = ['X' if i==1 else 'I' for i in range(7)]
    
                #    total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                #    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                #    self.stabs = corr_state[0][:]
                #    self.destabs = corr_state[1][:]
                    
            # If the measurement outcomes of both undefined stabilizers is 0, we're done.
           
            #print 'outflags =', outflags1, outflags2

            total_corr1 = ['I' for i in range(7)]
            if outflags1==1:
                syn_flags1 = []
                syn_flags1 += [self.run_one_circ(7).values()[0][0]]
                syn_flags1 += [self.run_one_circ(8).values()[0][0]]
                
                #print 'synflags1 =', syn_flags1

                if sum(syn_flags1) > 0:
                    if syn_flags1[0] == 0:
                        total_corr1[1] = 'Z'
                        total_corr1[2] = 'Z'
                    else:
                        if syn_flags1[1] == 0:
                            total_corr1[1] = 'Z'
                        else:
                            total_corr1[6] = 'Z'

            total_corr2 = ['I' for i in range(7)]
            if outflags2==1:
                syn_flags2 = []
                syn_flags2 += [self.run_one_circ(16).values()[0][0]]
                syn_flags2 += [self.run_one_circ(17).values()[0][0]]
                if sum(syn_flags2) > 0:
                    if syn_flags2[0] == 0:
                        total_corr2[1] = 'Z'
                        total_corr2[2] = 'Z'
                    else:
                        if syn_flags2[1] == 0:
                            total_corr2[1] = 'Z'
                        else:
                            total_corr2[6] = 'Z'
           
            if (outflags1==1) or (outflags2==1):
                if (total_corr1.count('Z') > 0) or (total_corr2.count('Z') > 0):
                    total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                    #print 'total corr =', total_corr
                    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                    self.stabs = corr_state[0][:]
                    self.destabs = corr_state[1][:]
                
            return error_det
                
            ##########################################################################


            # First logical qubit
            error_det1 = False
            syn1 = []
            first_elem = 0
            for i in range(first_elem,first_elem+3):
                out_run = self.run_one_circ(i).values()
                syn1 += [out_run[0][0]]
                outflags1 += [out_run[1][0]]
                if (syn1[-1]==1 and i!=(first_elem+undefined_stab)) or outflags1[-1]==1:
                    # Really, if a flag gets triggered, we could STOP right there.
                    # Change it later on
                    error_det1 = True
                    break
            if error_det1:
                syn1 = []
                for i in range(3,6):
                    syn1 += [self.run_one_circ(i).values()[0][0]]
                   
            # Second logical qubit
            error_det2 = False
            syn2, outflags2 = [], []
            first_elem = 6
            for i in range(first_elem,first_elem+3):
                out_run = self.run_one_circ(i).values()
                syn2 += [out_run[0][0]]
                outflags2 += [out_run[1][0]]
                if (syn2[-1]==1 and i!=(first_elem+undefined_stab)) or outflags2[-1]==1:
                    error_det2 = True
                    break
            if error_det2:
                syn2 = []
                for i in range(9,12):
                    syn2 += [self.run_one_circ(i).values()[0][0]]

            error_det = error_det1 or error_det2
        
        outflags1 += [0 for i in range(3-len(outflags1))]
        outflags2 += [0 for i in range(3-len(outflags2))]
            
        #print 'error det 1 =', error_det1
        #print 'error det 2 =', error_det2

        #print 'error_det =', error_det
        #print 'syns =', syn1, syn2
        #print 'outflags =', outflags1, outflags2
    
        base_corr1 = ['I' for i in range(7)]
        base_corr2 = ['I' for i in range(7)]
        
        if sum(syn1)>0 and sum(syn2)>0:
            # if both syndromes are non-trivial, we always apply a correction
            # on qubit 1 to account for the fact that the measurement outcomes
            # on stabilizers (1,2,5,6) are random.
            base_corr1[error_index] = Pauli_error
            base_corr2[error_index] = Pauli_error
    
            syn1[1] = (syn1[1]+1)%2
            syn2[1] = (syn2[1]+1)%2

        # if the target flag was triggered in the previous step
        if sum(inflags1)>0:
            # if the target syndrome doesn't show an error, but the ancilla syndrome does,
            # then flip the random stabilizers.
            if sum(syn1) < sum(syn2):
                syn1[error_index] = 1
                syn2[error_index] = 0
        
        elif sum(inflags2)>0:
            if sum(syn1) > sum(syn2):
                syn1[error_index] = 0
                syn2[error_index] = 1
            
        new_corr1 = st.Code.total_lookup_table_one_flag[tuple(inflags1)][tuple(syn1)]
        new_corr1 = [Pauli_error if oper=='E' else oper for oper in new_corr1]
        total_corr1 = ['I' if base_corr1[i]==new_corr1[i] else Pauli_error 
                       for i in range(7)]
        new_corr2 = st.Code.total_lookup_table_one_flag[tuple(inflags2)][tuple(syn2)]
        new_corr2 = [Pauli_error if oper=='E' else oper for oper in new_corr2]
        total_corr2 = ['I' if base_corr2[i]==new_corr2[i] else Pauli_error 
                       for i in range(7)]
    
        #print 'total corr 1 =', total_corr1
        #print 'total corr 2 =', total_corr2

        total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
        self.stabs = corr_state[0][:]
        self.destabs = corr_state[1][:]

        #return error_det, outflags1, outflags2
        return error_det
    
    
    
    def run_jointQECZ_ion(self, error_det, inflags1, inflags2):
        '''
        error_det:  True if an error has already been detected
                    False otherwise.
                    If error_det is True, then we just run the
                    stabilizers 1 time nonFT.
        flags1:  the previous flags for the first logical qubit
        flags2:  the previous flags for the second logical qubit
        '''

        error_index = 1
        Pauli_error = 'X'
        n_subcirc = 14
        
        subcircs_run_QECnonanc = {}
        subcircs_run_QECanc = {}
        for i in range(n_subcirc):
            subcircs_run_QECnonanc[i] = 0
            subcircs_run_QECanc[i] = 0

        #print error_det, inflags1, inflags2
        
        # undefined stab is the index of the stabilizer that is undefined.
        # With the current numbering, it's the second stab (stab 1)
        undefined_stab = 1

        outflags1, outflags2 = [], []
        if error_det:
            # In this case, we just run the non-FT circuits
        
            # Because the order of the stabilizer measurements in the QECx during
            # the merging step is different from the standard order, we need to
            # change the flag orderings
            inflags1 = [inflags1[1], inflags1[2], inflags1[0]]
            inflags2 = [inflags2[1], inflags2[2], inflags2[0]]

            # First logical qubit
            syn1 = []
            for i in range(5,8):
                syn1 += [self.run_one_circ(i).values()[0][0]]
                subcircs_run_QECnonanc[i] = 1
            # Reorder 2 to 0
            self.run_one_circ(13)
            subcircs_run_QECnonanc[13] = 1

            # Second logical qubit
            syn2 = []
            for i in range(n_subcirc+5,n_subcirc+8):
                syn2 += [self.run_one_circ(i).values()[0][0]]
                subcircs_run_QECanc[i] = 1
            # Reorder 2 to 0
            self.run_one_circ(n_subcirc+13)
            subcircs_run_QECanc[13] = 1
            
            outflags1, outflags2 = [0,0,0], [0,0,0]


        else:
            # In this case, both inflags are [0,0,0]

            ##########################################################################
            # Alternative:  only measure the undefined stabilizers
            error_det1 = False
            syn1 = []
            out_run1 = self.run_one_circ(1).values()
            subcircs_run_QECnonanc[1] = 1
            # Because the sub-circuit has 6 MS gates, we flip the outcomes
            syn1 += [(out_run1[0][0]+1)%2]
            outflags1 = (out_run1[1][0]+1)%2
            if outflags1==1:  error_det1 = True

            #brow.from_circuit(self.circuits[1], True)

            #print 'syn1 =', syn1
            #print 'outflags1 =', outflags1

            error_det2 = False
            syn2 = []
            out_run2 = self.run_one_circ(n_subcirc+1).values()
            subcircs_run_QECanc[1] = 1
            # Because the sub-circuit has 6 MS gates, we flip the outcomes
            syn2 += [(out_run2[0][0]+1)%2]
            outflags2 = (out_run2[1][0]+1)%2
            if outflags2==1:  error_det2 = True
            
            #brow.from_circuit(self.circuits[7], True)
            #sys.exit(0)

            #print 'syn2 =', syn2
            #print 'outflags2 =', outflags2
            
            # If a flag was triggered in either one of the stabs, error_det = True
            error_det = error_det1 or error_det2
            
            # If both syndromes were 1, correct both stabilizers
            if syn1[0]==1 and syn2[0]==1:
                total_corr1 = ['X' if i==1 else 'I' for i in range(7)]
                total_corr2 = ['X' if i==1 else 'I' for i in range(7)]
    
                total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                self.stabs = corr_state[0][:]
                self.destabs = corr_state[1][:]

                # if that flag was not triggered we reorder from stab 1 to stab 0
                # and we're done.
                if outflags1==0:
                    self.run_one_circ(12)
                    subcircs_run_QECnonanc[12] = 1
                if outflags2==0:
                    self.run_one_circ(n_subcirc+12)
                    subcircs_run_QECanc[12] = 1

            elif syn1[0]==0 and syn2[0]==0:
                # if that flag was not triggered we reorder from stab 1 to stab 0
                # and we're done.
                if outflags1==0:
                    self.run_one_circ(12)
                    subcircs_run_QECnonanc[12] = 1
                if outflags2==0:
                    self.run_one_circ(n_subcirc+12)
                    subcircs_run_QECanc[12] = 1

            # If the 2 syndromes are different, error_det = True and we measure
            # the other 2 stabilizer pairs with no flags.
            # First we measure the other boundary stabilizer
            elif syn1[0] != syn2[0]:
                error_det = True
                # reorder from 1 to 0
                self.run_one_circ(3)
                subcircs_run_QECnonanc[3] = 1
                syn1 += [self.run_one_circ(5).values()[0][0]]
                subcircs_run_QECnonanc[5] = 1
                # reorder from 1 to 0
                self.run_one_circ(n_subcirc+3)
                subcircs_run_QECanc[3] = 1
                syn2 += [self.run_one_circ(n_subcirc+5).values()[0][0]]
                subcircs_run_QECanc[5] = 1

                #print 'syn1 =', syn1
                #print 'syn2 =', syn2

                #brow.from_circuit(self.circuits[3], True)
                #brow.from_circuit(self.circuits[9], True)

                # If an error has been detected on both logical qubits, we
                # apply a correction on both qubits 1.
                if (1 in syn1) and (1 in syn2):
                    total_corr1 = ['X' if i==1 else 'I' for i in range(7)]
                    total_corr2 = ['X' if i==1 else 'I' for i in range(7)]
    
                    total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                    self.stabs = corr_state[0][:]
                    self.destabs = corr_state[1][:]
                
                else:
                    syn1 += [self.run_one_circ(7).values()[0][0]]
                    subcircs_run_QECnonanc[7] = 1
                    syn2 += [self.run_one_circ(n_subcirc+7).values()[0][0]]
                    subcircs_run_QECanc[7] = 1
                   
                    #print 'syn1 =', syn1
                    #print 'syn2 =', syn2

                    if (1 in syn1) and (1 in syn2):
                        total_corr1 = ['X' if i==1 else 'I' for i in range(7)]
                        total_corr2 = ['X' if i==1 else 'I' for i in range(7)]
    
                        total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                        self.stabs = corr_state[0][:]
                        self.destabs = corr_state[1][:]
                        
                        # if the flag was triggered, we reorder from 2 to 1
                        if outflags1==1:
                            self.run_one_circ(8)
                            subcircs_run_QECnonanc[8] = 1
                        if outflags2==1:
                            self.run_one_circ(n_subcirc+8)
                            subcircs_run_QECanc[8] = 1

                    
            # If the measurement outcomes of both undefined stabilizers is 0, we're done.
           
            #print 'outflags =', outflags1, outflags2

            total_corr1 = ['I' for i in range(7)]
            if outflags1==1:
                syn_flags1 = []
                syn_flags1 += [self.run_one_circ(10).values()[0][0]]
                subcircs_run_QECnonanc[10] = 1
                syn_flags1 += [self.run_one_circ(11).values()[0][0]]
                subcircs_run_QECnonanc[11] = 1
                # Reorder from stab 2 to 0
                self.run_one_circ(13)
                subcircs_run_QECnonanc[13] = 1

                #print 'synflags1 =', syn_flags1

                if sum(syn_flags1) > 0:
                    if syn_flags1[0] == 0:
                        total_corr1[1] = 'Z'
                        total_corr1[2] = 'Z'
                    else:
                        if syn_flags1[1] == 0:
                            total_corr1[1] = 'Z'
                        else:
                            total_corr1[6] = 'Z'

            total_corr2 = ['I' for i in range(7)]
            if outflags2==1:
                syn_flags2 = []
                syn_flags2 += [self.run_one_circ(n_subcirc+10).values()[0][0]]
                subcircs_run_QECanc[10] = 1
                syn_flags2 += [self.run_one_circ(n_subcirc+11).values()[0][0]]
                subcircs_run_QECanc[11] = 1
                # Reorder from stab 2 to 0
                self.run_one_circ(n_subcirc+13)
                subcircs_run_QECanc[13] = 1
                if sum(syn_flags2) > 0:
                    if syn_flags2[0] == 0:
                        total_corr2[1] = 'Z'
                        total_corr2[2] = 'Z'
                    else:
                        if syn_flags2[1] == 0:
                            total_corr2[1] = 'Z'
                        else:
                            total_corr2[6] = 'Z'
           
            if (outflags1==1) or (outflags2==1):
                if (total_corr1.count('Z') > 0) or (total_corr2.count('Z') > 0):
                    total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                    #print 'total corr =', total_corr
                    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                    self.stabs = corr_state[0][:]
                    self.destabs = corr_state[1][:]
                
            return error_det, subcircs_run_QECnonanc, subcircs_run_QECanc
                
            ##########################################################################


        
        outflags1 += [0 for i in range(3-len(outflags1))]
        outflags2 += [0 for i in range(3-len(outflags2))]
            
        #print 'error det 1 =', error_det1
        #print 'error det 2 =', error_det2

        #print 'error_det =', error_det
        #print 'syns =', syn1, syn2
        #print 'outflags =', outflags1, outflags2
    
        base_corr1 = ['I' for i in range(7)]
        base_corr2 = ['I' for i in range(7)]
        
        if sum(syn1)>0 and sum(syn2)>0:
            # if both syndromes are non-trivial, we always apply a correction
            # on qubit 1 to account for the fact that the measurement outcomes
            # on stabilizers (1,2,5,6) are random.
            base_corr1[error_index] = Pauli_error
            base_corr2[error_index] = Pauli_error
    
            syn1[1] = (syn1[1]+1)%2
            syn2[1] = (syn2[1]+1)%2

        # if the target flag was triggered in the previous step
        if sum(inflags1)>0:
            # if the target syndrome doesn't show an error, but the ancilla syndrome does,
            # then flip the random stabilizers.
            if sum(syn1) < sum(syn2):
                syn1[error_index] = 1
                syn2[error_index] = 0
        
        elif sum(inflags2)>0:
            if sum(syn1) > sum(syn2):
                syn1[error_index] = 0
                syn2[error_index] = 1
            
        new_corr1 = st.Code.total_lookup_table_one_flag[tuple(inflags1)][tuple(syn1)]
        new_corr1 = [Pauli_error if oper=='E' else oper for oper in new_corr1]
        total_corr1 = ['I' if base_corr1[i]==new_corr1[i] else Pauli_error 
                       for i in range(7)]
        new_corr2 = st.Code.total_lookup_table_one_flag[tuple(inflags2)][tuple(syn2)]
        new_corr2 = [Pauli_error if oper=='E' else oper for oper in new_corr2]
        total_corr2 = ['I' if base_corr2[i]==new_corr2[i] else Pauli_error 
                       for i in range(7)]
    
        #print 'total corr 1 =', total_corr1
        #print 'total corr 2 =', total_corr2

        total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
        self.stabs = corr_state[0][:]
        self.destabs = corr_state[1][:]

        #return error_det, outflags1, outflags2
        return error_det, subcircs_run_QECnonanc, subcircs_run_QECanc
   


    def run_jointQECX(self, error_det, inflags1, inflags2):
        '''
        error_det:  True if an error has already been detected
                    False otherwise.
                    If error_det is True, then we just run the
                    stabilizers 1 time nonFT.
        flags1:  the previous flags for the first logical qubit
        flags2:  the previous flags for the second logical qubit
        '''

        error_index = 3
        Pauli_error = 'Z'

        #print error_det, inflags1, inflags2
        
        # undefined stab is the index of the stabilizer that is undefined.
        # With the current numbering, it's the first stab (stab 0)
        undefined_stab = 0

        outflags1, outflags2 = [], []
        
        #print 'inflags1 =', inflags1
        #print 'inflags2 =', inflags2

        if error_det:
            # In this case, we just run the non-FT circuits
        
            # Because the order of the stabilizer measurements in the QECz during
            # the merging step is different from the standard order, we need to
            # change the flag orderings
            inflags1 = [inflags1[2], inflags1[0], inflags1[1]]
            inflags2 = [inflags2[2], inflags2[0], inflags2[1]]

            # First logical qubit
            syn1 = []
            for i in range(3,6):
                syn1 += [self.run_one_circ(i).values()[0][0]]
            
            # Second logical qubit
            syn2 = []
            for i in range(12,15):
                syn2 += [self.run_one_circ(i).values()[0][0]]
            
            outflags1, outflags2 = [0,0,0], [0,0,0]

            print 'syn1 =', syn1
            print 'syn2 =', syn2
            

        else:
            # In this case, both inflags are [0,0,0]

            ##########################################################################
            # Alternative:  only measure the undefined stabilizers
            error_det1 = False
            syn1 = []
            out_run1 = self.run_one_circ(0).values()
            syn1 += [out_run1[0][0]]
            outflags1 = out_run1[1][0]
            if outflags1==1:  error_det1 = True

            #brow.from_circuit(self.circuits[1], True)

            #print 'syn1 =', syn1
            #print 'outflags1 =', outflags1

            error_det2 = False
            syn2 = []
            out_run2 = self.run_one_circ(9).values()
            syn2 += [out_run2[0][0]]
            outflags2 = out_run2[1][0]
            if outflags2==1:  error_det2 = True
            
            #brow.from_circuit(self.circuits[7], True)
            #sys.exit(0)

            #print 'syn2 =', syn2
            #print 'outflags2 =', outflags2
            
            # If a flag was triggered in either one of the stabs, error_det = True
            error_det = error_det1 or error_det2
            
            # If both syndromes were 1, correct both stabilizers
            if syn1[0]==1 and syn2[0]==1:
                total_corr1 = ['Z' if i==3 else 'I' for i in range(7)]
                total_corr2 = ['Z' if i==3 else 'I' for i in range(7)]
    
                total_corr = total_corr1[:] + ['I' for i in range(7)] + total_corr2[:]
                corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                self.stabs = corr_state[0][:]
                self.destabs = corr_state[1][:]

            # If the 2 syndromes are different, error_det = True and we measure
            # the other 2 stabilizer pairs with no flags.
            elif syn1[0] != syn2[0]:
                error_det = True
                syn1 += [self.run_one_circ(5).values()[0][0]]
                #syn1 += [self.run_one_circ(5).values()[0][0]]
                syn2 += [self.run_one_circ(14).values()[0][0]]
                #syn2 += [self.run_one_circ(11).values()[0][0]]

                #print 'syn1 =', syn1
                #print 'syn2 =', syn2

                #brow.from_circuit(self.circuits[3], True)
                #brow.from_circuit(self.circuits[9], True)

                # If an error has been detected on both logical qubits, we
                # apply a correction on both qubits 3.
                if (1 in syn1) and (1 in syn2):
                    total_corr1 = ['Z' if i==3 else 'I' for i in range(7)]
                    total_corr2 = ['Z' if i==3 else 'I' for i in range(7)]
    
                    total_corr = total_corr1[:] + ['I' for i in range(7)] + total_corr2[:]
                    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                    self.stabs = corr_state[0][:]
                    self.destabs = corr_state[1][:]
                
                else:
                    syn1 += [self.run_one_circ(4).values()[0][0]]
                    syn2 += [self.run_one_circ(13).values()[0][0]]
                   
                    #print 'syn1 =', syn1
                    #print 'syn2 =', syn2

                    if (1 in syn1) and (1 in syn2):
                        total_corr1 = ['Z' if i==3 else 'I' for i in range(7)]
                        total_corr2 = ['Z' if i==3 else 'I' for i in range(7)]
    
                        total_corr = total_corr1[:] + ['I' for i in range(7)] + total_corr2[:]
                        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                        self.stabs = corr_state[0][:]
                        self.destabs = corr_state[1][:]

            
            # The next part is done to correct the possible hook errors that
            # occurred when we measured the first stabilizer with the flags

            total_corr1 = ['I' for i in range(7)]
            if outflags1==1:
                syn_flags1 = []
                syn_flags1 += [self.run_one_circ(6).values()[0][0]]
                syn_flags1 += [self.run_one_circ(8).values()[0][0]]
                
                #print 'synflags1 =', syn_flags1

                if sum(syn_flags1) > 0:
                    if syn_flags1[0] == 0:
                        total_corr1[3] = 'X'
                        total_corr1[4] = 'X'
                    else:
                        if syn_flags1[1] == 0:
                            total_corr1[3] = 'X'
                        else:
                            total_corr1[6] = 'X'
                    
                    total_corr = total_corr1 + ['I' for i in range(2*7)]
                    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                    self.stabs = corr_state[0][:]
                    self.destabs = corr_state[1][:]
                    

            #total_corr2 = ['I' for i in range(7)]
            #if outflags2==1:
            #    syn_flags2 = []
            #    syn_flags2 += [self.run_one_circ(16).values()[0][0]]
            #    syn_flags2 += [self.run_one_circ(17).values()[0][0]]
            #    if sum(syn_flags2) > 0:
            #        if syn_flags2[0] == 0:
            #            total_corr2[1] = 'Z'
            #            total_corr2[2] = 'Z'
            #        else:
            #            if syn_flags2[1] == 0:
            #                total_corr2[1] = 'Z'
            #            else:
            #                total_corr2[6] = 'Z'
           
            #if (outflags1==1) or (outflags2==1):
            #    if (total_corr1.count('Z') > 0) or (total_corr2.count('Z') > 0):
            #        total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
            #        #print 'total corr =', total_corr
            #        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
            #        self.stabs = corr_state[0][:]
            #        self.destabs = corr_state[1][:]
                
            return error_det
                
            ##########################################################################


        
        #outflags1 += [0 for i in range(3-len(outflags1))]
        #outflags2 += [0 for i in range(3-len(outflags2))]
            
        #print 'error det 1 =', error_det1
        #print 'error det 2 =', error_det2

        #print 'error_det =', error_det
        #print 'syns =', syn1, syn2
        #print 'outflags =', outflags1, outflags2
    
        base_corr1 = ['I' for i in range(7)]
        base_corr2 = ['I' for i in range(7)]
        
        if sum(syn1)>0 and sum(syn2)>0:
            # if both syndromes are non-trivial, we always apply a correction
            # on qubit 1 to account for the fact that the measurement outcomes
            # on stabilizers (1,2,5,6) are random.
            base_corr1[error_index] = Pauli_error
            base_corr2[error_index] = Pauli_error
    
            syn1[0] = (syn1[0]+1)%2
            syn2[0] = (syn2[0]+1)%2

            #print 'new syn1 =', syn1
            #print 'new syn2 =', syn2

        # if the target flag was triggered in the previous step
        if sum(inflags1)>0:
            # if the target syndrome doesn't show an error, but the ancilla syndrome does,
            # then flip the random stabilizers.
            if sum(syn1) < sum(syn2):
                syn1[0] = 1
                syn2[0] = 0
        
        elif sum(inflags2)>0:
            if sum(syn1) > sum(syn2):
                syn1[0] = 0
                syn2[1] = 1
            
        new_corr1 = st.Code.total_lookup_table_one_flag[tuple(inflags1)][tuple(syn1)]
        new_corr1 = [Pauli_error if oper=='E' else oper for oper in new_corr1]
        total_corr1 = ['I' if base_corr1[i]==new_corr1[i] else Pauli_error 
                       for i in range(7)]
        new_corr2 = st.Code.total_lookup_table_one_flag[tuple(inflags2)][tuple(syn2)]
        new_corr2 = [Pauli_error if oper=='E' else oper for oper in new_corr2]
        total_corr2 = ['I' if base_corr2[i]==new_corr2[i] else Pauli_error 
                       for i in range(7)]
    
        #print 'total corr 1 =', total_corr1
        #print 'total corr 2 =', total_corr2

        total_corr = total_corr1[:] + ['I' for i in range(7)] + total_corr2[:]
        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
        self.stabs = corr_state[0][:]
        self.destabs = corr_state[1][:]

        #return error_det, outflags1, outflags2
        return error_det



    def run_jointQECX_ion(self, error_det, inflags1, inflags2):
        '''
        error_det:  True if an error has already been detected
                    False otherwise.
                    If error_det is True, then we just run the
                    stabilizers 1 time nonFT.
        flags1:  the previous flags for the first logical qubit
        flags2:  the previous flags for the second logical qubit
        '''

        error_index = 3
        Pauli_error = 'Z'
        n_subcirc = 14
        
        subcircs_run_QECnonanc = {}
        subcircs_run_QECanc = {}
        for i in range(n_subcirc):
            subcircs_run_QECnonanc[i] = 0
            subcircs_run_QECanc[i] = 0

        #print error_det, inflags1, inflags2
        
        # undefined stab is the index of the stabilizer that is undefined.
        # With the current numbering, it's the first stab (stab 0)
        undefined_stab = 0

        outflags1, outflags2 = [], []
        
        #print 'inflags1 =', inflags1
        #print 'inflags2 =', inflags2

        if error_det:
            # In this case, we just run the non-FT circuits
        
            # Because the order of the stabilizer measurements in the QECz during
            # the merging step is different from the standard order, we need to
            # change the flag orderings
            inflags1 = [inflags1[2], inflags1[0], inflags1[1]]
            inflags2 = [inflags2[2], inflags2[0], inflags2[1]]

            # First logical qubit
            syn1 = []
            for i in range(5,8):
                syn1 += [self.run_one_circ(i).values()[0][0]]
                subcircs_run_QECnonanc[i] = 1
            # Reorder 2 to 0
            self.run_one_circ(13)
            subcircs_run_QECnonanc[13] = 1

            # Second logical qubit
            syn2 = []
            for i in range(n_subcirc+5,n_subcirc+8):
                syn2 += [self.run_one_circ(i).values()[0][0]]
                subcircs_run_QECanc[i] = 1
            # Reorder 2 to 0
            self.run_one_circ(n_subcirc+8)
            subcircs_run_QECanc[8] = 1
            
            outflags1, outflags2 = [0,0,0], [0,0,0]

            #print 'syn1 =', syn1
            #print 'syn2 =', syn2
            

        else:
            # In this case, both inflags are [0,0,0]

            ##########################################################################
            # Alternative:  only measure the undefined stabilizers
            error_det1 = False
            syn1 = []
            out_run1 = self.run_one_circ(0).values()
            subcircs_run_QECnonanc[0] = 1
            # Because the sub-circuit has 6 MS gates, we flip the outcomes
            syn1 += [(out_run1[0][0]+1)%2]
            outflags1 = (out_run1[1][0]+1)%2
            if outflags1==1:  error_det1 = True

            #brow.from_circuit(self.circuits[1], True)

            #print 'syn1 =', syn1
            #print 'outflags1 =', outflags1

            error_det2 = False
            syn2 = []
            out_run2 = self.run_one_circ(n_subcirc).values()
            subcircs_run_QECanc[0] = 1
            # Because the sub-circuit has 6 MS gates, we flip the outcomes
            syn2 += [(out_run2[0][0]+1)%2]
            outflags2 = (out_run2[1][0]+1)%2
            if outflags2==1:  error_det2 = True
            #syn2 += [self.run_one_circ(n_subcirc+5).values()[0][0]

            #brow.from_circuit(self.circuits[7], True)
            #sys.exit(0)

            #print 'syn2 =', syn2
            #print 'outflags2 =', outflags2
            
            # If a flag was triggered in either one of the stabs, error_det = True
            error_det = error_det1 or error_det2
            
            # If both syndromes were 1, correct both stabilizers
            if syn1[0]==1 and syn2[0]==1:
                total_corr1 = ['Z' if i==3 else 'I' for i in range(7)]
                total_corr2 = ['Z' if i==3 else 'I' for i in range(7)]
    
                total_corr = total_corr1[:] + ['I' for i in range(7)] + total_corr2[:]
                corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                self.stabs = corr_state[0][:]
                self.destabs = corr_state[1][:]

                # if that flag was not triggered we reorder from stab 1 to stab 0
                # and we're done.
                if outflags1==0:
                    self.run_one_circ(12)
                    subcircs_run_QECnonanc[12] = 1
                if outflags2==0:
                    self.run_one_circ(n_subcirc+8)
                    subcircs_run_QECanc[8] = 1

            elif syn1[0]==0 and syn2[0]==0:
                # if that flag was not triggered we reorder from stab 1 to stab 0
                # and we're done.
                if outflags1==0:
                    self.run_one_circ(12)
                    subcircs_run_QECnonanc[12] = 1
                if outflags2==0:
                    self.run_one_circ(n_subcirc+8)
                    subcircs_run_QECanc[8] = 1

            # If the 2 syndromes are different, error_det = True and we measure
            # the other 2 stabilizer pairs with no flags.
            elif syn1[0] != syn2[0]:
                error_det = True
                syn1 += [self.run_one_circ(7).values()[0][0]]
                subcircs_run_QECnonanc[7] = 1
                #syn1 += [self.run_one_circ(5).values()[0][0]]
                syn2 += [self.run_one_circ(n_subcirc+7).values()[0][0]]
                subcircs_run_QECanc[7] = 1
                #syn2 += [self.run_one_circ(11).values()[0][0]]

                #print 'syn1 =', syn1
                #print 'syn2 =', syn2

                #brow.from_circuit(self.circuits[3], True)
                #brow.from_circuit(self.circuits[9], True)

                # If an error has been detected on both logical qubits, we
                # apply a correction on both qubits 3.
                if (1 in syn1) and (1 in syn2):
                    total_corr1 = ['Z' if i==3 else 'I' for i in range(7)]
                    total_corr2 = ['Z' if i==3 else 'I' for i in range(7)]
    
                    total_corr = total_corr1[:] + ['I' for i in range(7)] + total_corr2[:]
                    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                    self.stabs = corr_state[0][:]
                    self.destabs = corr_state[1][:]
                
                else:
                    syn1 += [self.run_one_circ(6).values()[0][0]]
                    subcircs_run_QECnonanc[6] = 1
                    syn2 += [self.run_one_circ(n_subcirc+6).values()[0][0]]
                    subcircs_run_QECanc[6] = 1
                   
                    #print 'syn1 =', syn1
                    #print 'syn2 =', syn2

                    if (1 in syn1) and (1 in syn2):
                        total_corr1 = ['Z' if i==3 else 'I' for i in range(7)]
                        total_corr2 = ['Z' if i==3 else 'I' for i in range(7)]
    
                        total_corr = total_corr1[:] + ['I' for i in range(7)] + total_corr2[:]
                        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                        self.stabs = corr_state[0][:]
                        self.destabs = corr_state[1][:]

            
            # The next part is done to correct the possible hook errors that
            # occurred when we measured the first stabilizer with the flags

            total_corr1 = ['I' for i in range(7)]
            if outflags1==1:
                syn_flags1 = []
                syn_flags1 += [self.run_one_circ(9).values()[0][0]]
                subcircs_run_QECnonanc[9] = 1
                syn_flags1 += [self.run_one_circ(11).values()[0][0]]
                subcircs_run_QECnonanc[11] = 1
                # Reorder from 2 to 0
                self.run_one_circ(13)
                subcircs_run_QECnonanc[13] = 1

                #print 'synflags1 =', syn_flags1

                if sum(syn_flags1) > 0:
                    if syn_flags1[0] == 0:
                        total_corr1[3] = 'X'
                        total_corr1[4] = 'X'
                    else:
                        if syn_flags1[1] == 0:
                            total_corr1[3] = 'X'
                        else:
                            total_corr1[6] = 'X'
                    
                    total_corr = total_corr1 + ['I' for i in range(2*7)]
                    corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
                    self.stabs = corr_state[0][:]
                    self.destabs = corr_state[1][:]
                    

            # We don't have to worry about flags on the ancilla logical qubit
            # because the hook errors are X and we are going to measure in the
            # X basis.
            #total_corr2 = ['I' for i in range(7)]
            #if outflags2==1:
            #    syn_flags2 = []
            #    syn_flags2 += [self.run_one_circ(n_subcirc+9).values()[0][0]]
            #    syn_flags2 += [self.run_one_circ(n_subcirc+11).values()[0][0]]
                # Reorder from 2 to 0
            #    self.run_one_circ(n_subcirc+13)
            #    if sum(syn_flags2) > 0:
            #        if syn_flags2[0] == 0:
            #            total_corr2[1] = 'Z'
            #            total_corr2[2] = 'Z'
            #        else:
            #            if syn_flags2[1] == 0:
            #                total_corr2[1] = 'Z'
            #            else:
            #                total_corr2[6] = 'Z'
           
            #if (outflags1==1) or (outflags2==1):
            #    if (total_corr1.count('Z') > 0) or (total_corr2.count('Z') > 0):
            #        total_corr = ['I' for i in range(7)] + total_corr1[:] + total_corr2[:]
                    #print 'total corr =', total_corr
            #        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
            #        self.stabs = corr_state[0][:]
            #        self.destabs = corr_state[1][:]
                
            return error_det, subcircs_run_QECnonanc, subcircs_run_QECanc
                
            ##########################################################################


        
        #outflags1 += [0 for i in range(3-len(outflags1))]
        #outflags2 += [0 for i in range(3-len(outflags2))]
            
        #print 'error det 1 =', error_det1
        #print 'error det 2 =', error_det2

        #print 'error_det =', error_det
        #print 'syns =', syn1, syn2
        #print 'outflags =', outflags1, outflags2
    
        base_corr1 = ['I' for i in range(7)]
        base_corr2 = ['I' for i in range(7)]
        
        if sum(syn1)>0 and sum(syn2)>0:
            # if both syndromes are non-trivial, we always apply a correction
            # on qubit 1 to account for the fact that the measurement outcomes
            # on stabilizers (1,2,5,6) are random.
            base_corr1[error_index] = Pauli_error
            base_corr2[error_index] = Pauli_error
    
            syn1[0] = (syn1[0]+1)%2
            syn2[0] = (syn2[0]+1)%2

            #print 'new syn1 =', syn1
            #print 'new syn2 =', syn2

        # if the target flag was triggered in the previous step
        if sum(inflags1)>0:
            # if the target syndrome doesn't show an error, but the ancilla syndrome does,
            # then flip the random stabilizers.
            if sum(syn1) < sum(syn2):
                syn1[0] = 1
                syn2[0] = 0
        
        elif sum(inflags2)>0:
            if sum(syn1) > sum(syn2):
                syn1[0] = 0
                syn2[1] = 1
            
        new_corr1 = st.Code.total_lookup_table_one_flag[tuple(inflags1)][tuple(syn1)]
        new_corr1 = [Pauli_error if oper=='E' else oper for oper in new_corr1]
        total_corr1 = ['I' if base_corr1[i]==new_corr1[i] else Pauli_error 
                       for i in range(7)]
        new_corr2 = st.Code.total_lookup_table_one_flag[tuple(inflags2)][tuple(syn2)]
        new_corr2 = [Pauli_error if oper=='E' else oper for oper in new_corr2]
        total_corr2 = ['I' if base_corr2[i]==new_corr2[i] else Pauli_error 
                       for i in range(7)]
    
        #print 'total corr 1 =', total_corr1
        #print 'total corr 2 =', total_corr2

        total_corr = total_corr1[:] + ['I' for i in range(7)] + total_corr2[:]
        corr_state = qfun.update_stabs(self.stabs, self.destabs, total_corr)
        self.stabs = corr_state[0][:]
        self.destabs = corr_state[1][:]

        #return error_det, outflags1, outflags2
        return error_det, subcircs_run_QECnonanc, subcircs_run_QECanc



class Supra_Circuit(object):
    '''
    a supra-circuit is composed of several quantum operations.
    '''

    def __init__(self, initial_state, circuit, code, chp_location,
                 bare_ancilla=False, ion=False):
        '''
        bare_ancilla refers to whether or not the ancillae are bare
        qubits.
        '''
   
        self.state = initial_state[:]
        self.quant_opers = circuit.gates[:]
        self.code = code
        self.chp_loc = chp_location
        self.bare = bare_ancilla
        self.error_det = False
        self.ion = ion
        self.total_subcircs_run = {}


    def run_one_oper(self, quant_gate):
        '''
        runs one quantum operation and returns the final state
        quant_gate should be a Gate object with subcircuits.
        '''
   
        sub_circ = quant_gate.circuit_list[0]
 
        if quant_gate.gate_name == 'JointQECZ_flags':

            quant_circs = [g.circuit_list[0] for g in sub_circ.gates]
            quant_circs0 = [g.circuit_list[0] for g in quant_circs[0].gates]
            quant_circs1 = [g.circuit_list[0] for g in quant_circs[1].gates]
            reordered_quant_circs = quant_circs0 + quant_circs1

            q_oper = QEC_with_flags(self.state[:], reordered_quant_circs, self.chp_loc) 
            #output = q_oper.run_jointQECZ(self.error_det, self.flags1[-1], self.flags2[-1])
            if self.ion:
                output = q_oper.run_jointQECZ_ion(self.error_det, self.flags1, self.flags2)
            else:
                output = q_oper.run_jointQECZ(self.error_det, self.flags1, self.flags2)
            
            self.state = [q_oper.stabs[:], q_oper.destabs[:]]
            
            return output

        elif quant_gate.gate_name == 'JointQECX_flags':

            quant_circs = [g.circuit_list[0] for g in sub_circ.gates]
            quant_circs0 = [g.circuit_list[0] for g in quant_circs[0].gates]
            quant_circs1 = [g.circuit_list[0] for g in quant_circs[1].gates]
            reordered_quant_circs = quant_circs0 + quant_circs1

            q_oper = QEC_with_flags(self.state[:], reordered_quant_circs, self.chp_loc) 
            #output = q_oper.run_jointQECZ(self.error_det, self.flags1[-1], self.flags2[-1])
            if self.ion:
                output = q_oper.run_jointQECX_ion(self.error_det, self.flags1, self.flags2)
            else:
                output = q_oper.run_jointQECX(self.error_det, self.flags1, self.flags2)
            
            self.state = [q_oper.stabs[:], q_oper.destabs[:]]
            
            return output


        elif quant_gate.gate_name[:8] == 'JointQEC':
           
            #print 'Running JointQEC'

            stab_kind = quant_gate.gate_name[-1]
        
            quant_circs = [g.circuit_list[0] for g in sub_circ.gates]
            quant_circs0 = [g.circuit_list[0] for g in quant_circs[0].gates]
            quant_circs1 = [g.circuit_list[0] for g in quant_circs[1].gates]
            reordered_quant_circs = quant_circs0 + quant_circs1

            q_oper = QEC_d3(self.state[:], reordered_quant_circs, self.chp_loc)
            n_rep1, n_rep2 = q_oper.run_jointQEC(stab_kind)
            
            self.state = [q_oper.stabs[:], q_oper.destabs[:]]
        
            return n_rep1, n_rep2 


        elif quant_gate.gate_name[-7:] == 'Correct':
            
            quant_circs = [g.circuit_list[0] for g in sub_circ.gates]
            q_oper = QEC_d3(self.state[:], quant_circs, self.chp_loc)
            
            if self.code == 'Steane':
                n_rep = q_oper.run_fullQEC_CSS_d3(self.code, self.bare)
            
            elif self.code=='Cross' or self.code=='5qubit':
                n_rep = q_oper.run_fullQEC_nonCSS(self.code, self.bare)
        
            self.state = [q_oper.stabs[:], q_oper.destabs[:]]

            return n_rep

        
        elif quant_gate.gate_name[:16] == 'Measure2logicals':
            
            quant_circs = [g.circuit_list[0] for g in sub_circ.gates]
            #faulty_gate = quant_circs[0].gates[2]
            #if faulty_gate.gate_name == 'CX':
                #faulty_qubit = faulty_gate.qubits[1]
                #err_g = quant_circs[0].insert_gate(faulty_gate, [faulty_qubit], '', 'Z', False)
                #brow.from_circuit(quant_circs[0], True)
            q_oper = Measure_2_logicals(self.state[:], quant_circs, self.chp_loc)
            #parity, n_rep1, n_rep2 = q_oper.run_all(quant_gate.gate_name[16])

            if len(quant_circs) < 10:
                parity, corr_type = q_oper.run_all(quant_gate.gate_name[-1])
                #print 'parity =', parity
                #print 'corr type =', corr_type

                self.state = [q_oper.stabs[:], q_oper.destabs[:]]
                return parity, corr_type
            
            else:
                parity, rep, corr_type, anc_corr = q_oper.run_all_long(quant_gate.gate_name[-1])
                #print 'parity =', parity
                #print 'corr type =', corr_type
                #print 'anc corr =', anc_corr

                self.state = [q_oper.stabs[:], q_oper.destabs[:]]

                #print self.state[0]
                #sys.exit(0)
                return parity, rep, corr_type, anc_corr


        elif quant_gate.gate_name=='MeasureXX_flags' or quant_gate.gate_name=='MeasureZZ_flags':

            #print 'Running Measure %s%s' %(quant_gate.gate_name[7],quant_gate.gate_name[8])

            quant_circs = [g.circuit_list[0] for g in sub_circ.gates]
            q_oper = Measure_2_logicals(self.state[:], quant_circs, self.chp_loc)
            #output = q_oper.run_boundary_oper_flags(self.error_det)
            if self.ion:
                output = q_oper.run_boundary_oper_flags_short_ion(self.error_det)
            else:
                output = q_oper.run_boundary_oper_flags_short(self.error_det)
            self.state = [q_oper.stabs[:], q_oper.destabs[:]]
            #print self.state[0] 

            return output


        else:
            # assume it's just transversal logical gate or something that
            # doesn't require feedback based on measurements.
            
            q_oper = Quantum_Operation(self.state[:], [sub_circ], self.chp_loc)
            if quant_gate.gate_name[:7] == 'Measure':
                sub_circ.to_ancilla([q.qubit_id for q in sub_circ.qubits()])
                q_oper.n_d_q = len(q_oper.stabs) - len(sub_circ.ancilla_qubits())
            
            #brow.from_circuit(sub_circ, True)

            output_dict = q_oper.run_one_circ(0)            

            self.state = [q_oper.stabs[:], q_oper.destabs[:]]
            
            return output_dict



class CNOT_latt_surg(Supra_Circuit):
    '''
    '''

    def run_all_gates(self):
        '''
        '''

        n_repEC = []
        
        gate_i = 0

        #print 'Initial state =', self.state[0]


        for q_oper in self.quant_opers:
            #print q_oper.gate_name
            output = self.run_one_oper(q_oper)
           
            if q_oper.gate_name == 'JointQECZ_flags':
                self.error_det = output[0]
                subcircs_run = {0:output[1], 1:output[2]}
                self.total_subcircs_run[1] = subcircs_run
                #for stab in self.state[0]:
                #    if 'Z' in stab:
                #        print stab

            elif q_oper.gate_name == 'JointQECX_flags':
                self.error_det = output[0]
                subcircs_run = {0:output[1], 1:output[2]}
                self.total_subcircs_run[3] = subcircs_run

            elif q_oper.gate_name[:8] == 'JointQEC':
                n_rep1, n_rep2 = output
                #print 'State after JoinQEC:'
                #print self.state[0]
                #sys.exit(0)


            elif q_oper.gate_name[-7:] == 'Correct':
                n_repEC += [output]
                #print 'State after EC =', self.state[0]
            
            elif q_oper.gate_name == 'Measure2logicalsX':
                parX = output[0]
                corr_type = output[1]
                #n_rep1X, n_rep2X = output[1]
                #print self.state[0]
                clause1 = corr_type == 'normal' and parX == 1
                clause2 = corr_type == 'alternative' and parX == 0
                if clause1 or clause2:
                    #print 'Z correction after M_xx? Yes'
                    # Z logical on control
                    Z_corr = ['Z' for i in range(7)] + ['I' for i in range(7*2)]
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   Z_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]
                    #print 'State after corr:'
                    #print self.state[0]

            
            elif q_oper.gate_name == 'MeasureXX_flags':
                #print output
                low_w, high_w = output[0], output[1]
                parM = (low_w[-1]+high_w[-1])%2
                self.flags1, self.flags2 = output[4], output[5]
                self.error_det = output[8]
                corr_targ, corr_anc = output[6], output[7]
                self.total_subcircs_run[0] = output[9]
                #print corr_targ, corr_anc
                #print 'parity =', parM
                #print output[8]
                clause1 = corr_targ == 'normal' and parM == 1
                clause2 = corr_targ == 'alternative' and parM == 0
               
                # if either clause is True, we apply logical Z on ctrl
                if (clause1 or clause2):
                    log_corr = ['Z' for i in range(7)] + ['I' for i in range(7*2)]
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   log_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]
                   
                # if the ancilla correction is alternative, there was an error
                # on the ancilla qubit before the measurement of XX and we apply
                # logical Z on the ancilla
                if corr_anc == 'alternative':
                    log_corr = ['I' for i in range(7*2)] + ['Z' for i in range(7)]
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   log_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]
                    
            
            elif q_oper.gate_name == 'MeasureZZ_flags':
                #print output
                low_w, high_w = output[0], output[1]
                parM = (low_w[-1]+high_w[-1])%2
                self.flags1, self.flags2 = output[4], output[5]
                self.error_det = output[8]
                corr_targ, corr_anc = output[6], output[7]
                self.total_subcircs_run[2] = output[9]
                #print corr_targ, corr_anc
                #print 'parity =', parM
                #print output[8]
                clause1 = corr_targ == 'normal' and parM == 1
                clause2 = corr_targ == 'alternative' and parM == 0

                #print 'corr_targ =', corr_targ
                #print 'parM =', parM

                # if either clause is True, we apply logical X on target
                if (clause1 or clause2):
                    log_corr = ['I' for i in range(7)]
                    log_corr += ['X' for i in range(7)]
                    log_corr += ['I' for i in range(7)]
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   log_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]
                  
                    #print 'Applied X correction on target'

                # if the ancilla correction is alternative, there was an error
                # on the ancilla qubit before the measurement of ZZ and we apply
                # logical X on the ancilla
                
                #print 'corr_anc =', corr_anc

                if corr_anc == 'alternative':
                    log_corr = ['I' for i in range(7)]
                    log_corr += ['X' for i in range(7)]
                    log_corr += ['I' for i in range(7)]
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   log_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]

                #print low_w, high_w
                #for stab in self.state[0]:
                #    if 'Z' in stab:
                #        print stab
                


            elif q_oper.gate_name[:20] == 'Measure2logicalslong':
                stab = q_oper.gate_name[-1]
                parM = output[0]
                rep_low, rep_high = output[1]
                corr_type = output[2]
                anc_corr = output[3]
                #print self.state[0]
                
                # decide whether or not to apply the logical correction
                clause1 = corr_type == 'normal' and parM == 1
                clause2 = corr_type == 'alternative' and parM == 0
                if (clause1 or clause2):
                    #print 'Correction after M_%s%s? Yes' %(stab, stab)
                    
                    # if stab is X, we apply Z logical on control
                    if stab == 'X':
                        log_corr = ['Z' for i in range(7)] + ['I' for i in range(7*2)]
                    # if stab is Z, we apply X logical on target
                    elif stab == 'Z':
                        log_corr = ['I' for i in range(7)]
                        log_corr += ['X' for i in range(7)]
                        log_corr += ['I' for i in range(7)]
                        
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   log_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]
                    #print 'State after corr:'
                    #print self.state[0]

                # decide whether or not to apply the logical operator on the ancilla
                if anc_corr:
                    #print 'Correction after M_%s%s on ancilla? Yes' %(stab, stab)
                    # if stab is X, Z logical on ancilla
                    if stab == 'X':
                        log_corr = ['I' for i in range(7*2)] + ['Z' for i in range(7)]
                    # if stab is Z, X logical on target
                    elif stab == 'Z':
                        log_corr = ['I' for i in range(7)]
                        log_corr += ['X' for i in range(7)]
                        log_corr += ['I' for i in range(7)]
                        
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   log_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]
                    #print 'State after corr:'
                    #print self.state[0]



            elif q_oper.gate_name == 'Measure2logicalsZ':
                parZ = output[0]
                corr_type = output[1]
                #n_rep1Z, n_rep2Z = output[1]
                #print self.state[0]
                clause1 = corr_type == 'normal' and parZ == 1
                clause2 = corr_type == 'alternative' and parZ == 0
                if clause1 or clause2:
                    #print 'X correction after M_zz? Yes'
                    # X logical on target
                    X_corr = ['I' for i in range(7)]
                    X_corr += ['X' for i in range(7)]
                    X_corr += ['I' for i in range(7)]
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   X_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]
                    #print 'State after corr:'
                    #print self.state[0]

            elif q_oper.gate_name == 'MeasureX':
                meas_dict = output
                meas_outcomes = [val[0] for val in meas_dict.values()]
                parXanc = st.Code.parity_meas_Steane_EC(meas_outcomes)
                #print self.state[0]
                if parXanc == 1:
                    #print 'Z correction after M_x?  Yes'
                    # Z logical on control
                    Z_corr = ['Z' for i in range(7)] + ['I' for i in range(7)]
                    corr_state = qfun.update_stabs(self.state[0][:],
                                                   self.state[1][:],
                                                   Z_corr)
                    self.state = [corr_state[0][:], corr_state[1][:]]
                    #print 'State after corr:'
                    #print self.state[0]

                #sys.exit(0)
        return None

