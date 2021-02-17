import sys
import os
from circuit import *
from visualizer import browser_vis as brow



class Surface_Code_ion_trap:
    '''
    Measure stabilizers with a bare ancilla for each stabilizer, no cat states.
    Class defined to implement surface-17 code.
    The gates are already compiled to the primitve gates for ion traps.
    '''


    @classmethod
    def insert_I_errors_to_circ(cls, whole_circ, Is_after_2q,
                                Is_after_1q, MS_heating,
                                Stark, qubit_assignment,
                                time_MS, time_1q):
        '''
        '''
        if Is_after_1q:
            for g in whole_circ.gates[::-1]:
                if (len(g.qubits) == 1) and (g.gate_name[0] != 'I'):
                    new_g = whole_circ.insert_gate(g, g.qubits, '', 'I', False)

        if Is_after_2q:
            for g in whole_circ.gates[::-1]:
                if (len(g.qubits) == 2) and (g.gate_name[:2] != 'II'):
                    new_g = whole_circ.insert_gate(g, [g.qubits[0]], '', 'I', False)
                    new_g = whole_circ.insert_gate(g, [g.qubits[1]], '', 'I', False)

        if MS_heating:
            for g in whole_circ.gates[::-1]:
                if g.gate_name[:2] == 'MS':
                    x_ion1 = qubit_assignment.index(g.qubits[0].qubit_id)  # ion1 position
                    x_ion2 = qubit_assignment.index(g.qubits[1].qubit_id)  # ion2 position
                    ion_distance = abs(x_ion1 - x_ion2)   # distance between ions
                    MS_duration = time_MS(ion_distance)   # duration of MS gate
                    new_g = whole_circ.insert_gate(g, g.qubits, '', 'II_heat', False)
                    new_g.duration = MS_duration 

        if Stark:
            for g in whole_circ.gates[::-1]:
                if g.gate_name[:2] == 'MS':
                    x_ion1 = qubit_assignment.index(g.qubits[0].qubit_id)  # ion1 position
                    x_ion2 = qubit_assignment.index(g.qubits[1].qubit_id)  # ion2 position
                    ion_distance = abs(x_ion1 - x_ion2)   # distance between ions
                    MS_duration = time_MS(ion_distance)   # duration of MS gate
                    new_g1 = whole_circ.insert_gate(g, [g.qubits[0]], '', 'I_stark', False)
                    new_g1.duration = MS_duration
                    new_g2 = whole_circ.insert_gate(g, [g.qubits[1]], '', 'I_stark', False)
                    new_g2.duration = MS_duration
                elif (len(g.qubits) == 1) and (g.gate_name[0] != 'I') and (g.gate_name[0] != 'P'):
                    new_g = whole_circ.insert_gate(g, [g.qubits[0]], '', 'I_stark', False)
                    new_g.duration = time_1q

        return whole_circ
                    
                    


    @classmethod
    def generate_one_stab_ion_trap(cls, stabilizer, n_total, 
                                   i_ancilla, meas_errors=True):
        '''
        generates the circuit corresponding to 1 stabilizer
        measurement with ion-trap gates.
        The only I-like errors are measurements.
        '''
        stab_kind = stabilizer[0][0]    # 'X' or 'Z'
        stab_weight = len(stabilizer)

        stab_circ = Circuit()
        if len(stabilizer) == 2:
            stab_circ.add_gate_at([i_ancilla], 'PrepareZMinus')
        elif len(stabilizer) == 4:
            stab_circ.add_gate_at([i_ancilla], 'PrepareZPlus')
           
        if stab_kind == 'Z':
            for oper in stabilizer:
                stab_circ.add_gate_at([oper[1]], 'RY +')
                stab_circ.add_gate_at([oper[1]], 'RX -')
                stab_circ.add_gate_at([oper[1], i_ancilla], 'MS')
                stab_circ.add_gate_at([oper[1]], 'RY -')

        elif stab_kind == 'X':
            for oper in stabilizer:
                stab_circ.add_gate_at([oper[1]], 'RX +')
                stab_circ.add_gate_at([oper[1], i_ancilla], 'MS')

        if meas_errors:
            stab_circ.add_gate_at([i_ancilla], 'ImZ')
        stab_circ.add_gate_at([i_ancilla], 'MeasureZ')
        
        return stab_circ
        


    @classmethod
    def generate_one_stab_abstract(cls, stabilizer, n_total,
                                   i_ancilla, meas_errors=True,
                                   Is_after_two=False):
        '''
        '''
        stab_kind = stabilizer[0][0]    # 'X' or 'Z'
        stab_weight = len(stabilizer)

        stab_circ = Circuit()
        stab_circ.add_gate_at([i_ancilla], 'PrepareZPlus')
        
        if stab_kind == 'Z':
            for oper in stabilizer:
                stab_circ.add_gate_at([oper[1], i_ancilla], 'CX')
                if Is_after_two:
                    stab_circ.add_gate_at([oper[1]], 'I')
                    stab_circ.add_gate_at([i_ancilla], 'I')
                 
        elif stab_kind == 'X':
            stab_circ.add_gate_at([i_ancilla], 'H')
            for oper in stabilizer:
                stab_circ.add_gate_at([i_ancilla, oper[1]], 'CX')
                if Is_after_two:
                    stab_circ.add_gate_at([oper[1]], 'I')
                    stab_circ.add_gate_at([i_ancilla], 'I')
            stab_circ.add_gate_at([i_ancilla], 'H')

        if meas_errors:
            stab_circ.add_gate_at([i_ancilla], 'ImZ')
        stab_circ.add_gate_at([i_ancilla], 'MeasureZ')

        return stab_circ

   
 
    @classmethod
    def generate_stabs_meas(cls, stabilizer_list, ancillae_indexes,
                            input_I_round=False, meas_errors=True,
                            ion_trap=False, add_errors=False,
                            Is_after_2q=False, Is_after_1q=False,
                            MS_heating=False, Stark=False,
                            qubit_assignment=None, time_MS=None,
                            time_1q=None): 
        '''
        add_errors should only be True for debugging purposes.
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
       
        n_data = 9
        n_stabs = len(stabilizer_list)
        n_steps = n_stabs/2
        if ion_trap:
            stab_func = Surface_Code_ion_trap.generate_one_stab_ion_trap
        else:
            stab_func = Surface_Code_ion_trap.generate_one_stab_abstract
            
        n_total = n_data + n_stabs
        
        stabs_circ = Circuit()
        if input_I_round:
            for i in range(n_data):
                stabs_circ.add_gate_at([i], 'I')
        
        #if add_errors:
            #stabs_circ.add_gate_at([3], 'Y')
            #new_g = stabs_circ.gates[0]
            #new_g.is_error = True      
            #stabs_circ.add_gate_at([1], 'X')
            #new_g = stabs_circ.gates[1]
            #new_g.is_error = True      
 
        gate_name = 'Two_parallel_stabilizers'
        for i in range(n_steps):
            circ1 = stab_func(stabilizer_list[2*i], n_total, 
                              ancillae_indexes[2*i],
                              meas_errors)
            #circ1.to_ancilla([ancillae_indexes[2*i]])
            circ2 = stab_func(stabilizer_list[2*i+1], n_total, 
                              ancillae_indexes[2*i+1],
                              meas_errors)
            #circ2.to_ancilla([ancillae_indexes[2*i+1]])
            circ1.join_circuit(circ2, True)
            circ1 = Surface_Code_ion_trap.insert_I_errors_to_circ(
                                            circ1, Is_after_2q, Is_after_1q, 
                                            MS_heating, Stark, qubit_assignment,
                                            time_MS, time_1q)
            #if i == 0:
            #    for g in circ1.gates:
            #        print g.gate_name, [q.qubit_id for q in g.qubits], g.duration

            two_stab_circ = Encoded_Gate(gate_name,[circ1]).circuit_wrap()
            stabs_circ.join_circuit(two_stab_circ, True)

            

        # if we set the last qubits to be ancillae at this point
        # we get an error.  Need to do it early.
        #stabs_circ.to_ancilla(range(n_data, n_total))

        #brow.from_circuit(stabs_circ, True)
        #sys.exit(0)        
        
        return stabs_circ
            


    @classmethod
    def generate_logical_state_surf17(cls, stabilizer_list, ancillae_indexes,
                                      logical_oper):
        '''
        Generates circuit to measure the stabilizers and the logical operator.
        '''

        stabs_circ = Surface_Code_ion_trap.generate_stabs_meas(stabilizer_list, 
                                                               ancillae_indexes)
        last_anc = max(ancillae_indexes) + 1
        log_circ = Circuit()
        log_circ.add_gate_at([last_anc], 'PrepareZPlus')
        log_circ.add_gate_at([last_anc], 'H')
        
        for oper in logical_oper:
            log_circ.add_gate_at([last_anc, oper[1]], 'C%s' %oper[0])
        
        log_circ.add_gate_at([last_anc], 'MeasureX')
        stabs_circ.join_circuit(log_circ, True)
        stabs_circ = stabs_circ.unpack()
        stabs_circ.to_ancilla(range(9,18))       
 
        return stabs_circ 
        
                                        

class Flag_Correct:
    '''
    Measure stabilizers with a bare ancilla and flags, not cat states.
    Class defined to implement the 4.8.8 color codes
    '''

    @classmethod
    def generate_Reichardt_d3_1_flag(cls, meas_errors, Is_after_two, stabs='Z',
                                     with_flag=True, initial_I=False):
        '''
        Specific to the circuit in Figure 8b in arXiv:1705.02329
        '''
        n_data = 7
        n_anc = 4
        CNOT_indices = [[6,7], [10,7], [7,9], [5,7], [7,8], [4,7], [3,7],
                        [1,8], [2,8], [4,9], [0,9], [2,9], [10,7], [10,8], [10,9]]
        Z_indices = [7,8,9]
        X_indices = [10]
        
        if not with_flag:
            n_anc = 3
            CNOT_indices = [[6,7], [7,9], [5,7], [7,8], [4,7], [3,7],
                            [1,8], [2,8], [4,9], [0,9], [2,9]]
            X_indices = []

        if stabs == 'X':
            CNOT_indices = [[CNOT_i[1], CNOT_i[0]] for CNOT_i in CNOT_indices]
            Z_indices = X_indices[:]
            X_indices = [7,8,9]

        Reich_circ = Circuit()
        if initial_I:
            for i in range(n_data):
                Reich_circ.add_gate_at([i], 'I')

        for i in Z_indices:
            Reich_circ.add_gate_at([i], 'PrepareZPlus')
        for i in X_indices:
            Reich_circ.add_gate_at([i], 'PrepareXPlus')

        for CNOT_i in CNOT_indices:
            Reich_circ.add_gate_at(CNOT_i, 'CX')
            if Is_after_two:
                Reich_circ.add_gate_at([CNOT_i[0]], 'I')
                Reich_circ.add_gate_at([CNOT_i[1]], 'I')
        
        for i in Z_indices:
            if meas_errors:
                Reich_circ.add_gate_at([i], 'ImZ')
            Reich_circ.add_gate_at([i], 'MeasureZ')
        for i in X_indices:
            if meas_errors:
                Reich_circ.add_gate_at([i], 'ImX')
            Reich_circ.add_gate_at([i], 'MeasureX')
        
        Reich_circ.to_ancilla(range(n_data, n_data+n_anc))
        
        return Reich_circ


    
    @classmethod
    def generate_whole_QEC_Reichardt(cls, meas_errors, Is_after_two, n_rep=3, group_reps=False,
                                     initial_I=False):
        '''
        '''
        
        complete_circ = Circuit()
        for rep_i in range(n_rep):
            if rep_i == 0 and initial_I:
                circ1 = Flag_Correct.generate_Reichardt_d3_1_flag(meas_errors, Is_after_two, 'X',
                                                                  True, True)  
            else:
                circ1 = Flag_Correct.generate_Reichardt_d3_1_flag(meas_errors, Is_after_two, 'X') 
            circ2 = Flag_Correct.generate_Reichardt_d3_1_flag(meas_errors, Is_after_two, 'Z')
            
            circ1 = Encoded_Gate('Stabilizers_X%i'%rep_i, [circ1]).circuit_wrap()
            circ2 = Encoded_Gate('Stabilizers_Z%i'%rep_i, [circ2]).circuit_wrap()
            circ1.join_circuit(circ2)
            
            if group_reps:
                circ1 = Encoded_Gate('Rep_%i'%rep_i, [circ1]).circuit_wrap()
            complete_circ.join_circuit(circ1)

        return complete_circ


    
    @classmethod
    def generate_whole_QEC_Reichardt_special(cls, meas_errors, Is_after_two, initial_I=False):
        '''
        special circuit to run the Reichardt circuit using only 1 flag.
        The circuit consists of Sx(f), Sz(f), ( Sx(b), Sz(b) )
        '''
        
        QEC_circ = Flag_Correct.generate_whole_QEC_Reichardt(meas_errors, Is_after_two, 1,
                                                             False, initial_I)

        circ1 = Flag_Correct.generate_Reichardt_d3_1_flag(meas_errors, Is_after_two, 'X', False) 
        circ2 = Flag_Correct.generate_Reichardt_d3_1_flag(meas_errors, Is_after_two, 'Z', False)
        circ1 = Encoded_Gate('Stabilizers_X_bare', [circ1]).circuit_wrap()
        circ2 = Encoded_Gate('Stabilizers_Z_bare', [circ2]).circuit_wrap()
        circ1.join_circuit(circ2)
        QEC_circ.join_circuit(circ1)

        return QEC_circ



    @classmethod
    def generate_one_flagged_stab(cls, i_first_anc, stabilizer, flags, 
                                  meas_errors, Is_after_two, to_ancilla,
                                  initial_I=False):
        '''
        '''

        n_flags = len(flags)
        Pauli_oper = stabilizer[0][0]
        coupling_gate = 'C'+Pauli_oper

        stab_circ = Circuit()
        if initial_I:
            for i in range(i_first_anc):
                stab_circ.add_gate_at([i], 'I')

        stab_circ.add_gate_at([i_first_anc], 'PrepareXPlus')
        for i in range(n_flags):
            stab_circ.add_gate_at([i_first_anc+i+1], 'PrepareZPlus')

        for i in range(len(stabilizer)):
            for flag in flags:
                if i in flag:
                    flag_i = flags.index(flag)
                    stab_circ.add_gate_at([i_first_anc, i_first_anc+flag_i+1], 'CX')
                    if Is_after_two:
                        stab_circ.add_gate_at([i_first_anc], 'I')
                        stab_circ.add_gate_at([i_first_anc+flag_i+1], 'I')

            stab_circ.add_gate_at([i_first_anc, stabilizer[i][1]], coupling_gate)
            if Is_after_two:
                stab_circ.add_gate_at([i_first_anc], 'I')
                stab_circ.add_gate_at([stabilizer[i][1]], 'I')

        if meas_errors:
            stab_circ.add_gate_at([i_first_anc], 'ImX')
        stab_circ.add_gate_at([i_first_anc], 'MeasureX')
        for i in range(n_flags):
            if meas_errors:
                stab_circ.add_gate_at([i_first_anc+i+1], 'ImZ')
            stab_circ.add_gate_at([i_first_anc+i+1], 'MeasureZ')

        if to_ancilla:
            stab_circ.to_ancilla(range(i_first_anc, i_first_anc+n_flags+1))

        return stab_circ



    @classmethod
    def generate_all_flagged_stabs(cls, stabilizer_list, flags_list,
                                   meas_errors, Is_after_two, n_data):
        '''
        '''

        stabs_circ = Circuit()
        for i in range(len(stabilizer_list)):

            stabs_circ1 = Flag_Correct.generate_one_flagged_stab(n_data,
                                                            stabilizer_list[i],
                                                            flags_list[i],
                                                            meas_errors,
                                                            Is_after_two,
                                                            True)
        
            stabs_circ.join_circuit(stabs_circ1)
       

        return stabs_circ



    @classmethod
    def generate_high_indeterminacy_circuit(cls, stabilizer_list, flags_list,
                                            meas_errors, Is_after_two, n_data,
                                            initial_I=False):
        '''
        '''

        stabs_circ = Circuit()
        for i in range(len(stabilizer_list)):
            local_initial_I = False
            if i==0 and initial_I:  local_initial_I = True

            circ_flag = Flag_Correct.generate_one_flagged_stab(n_data,
                                                            stabilizer_list[i],
                                                            flags_list[i],
                                                            meas_errors,
                                                            Is_after_two,
                                                            True,
                                                            local_initial_I)
            circ_flag = Encoded_Gate('S%i_flag'%i, [circ_flag]).circuit_wrap()
            
            circ_bare = Flag_Correct.generate_one_flagged_stab(n_data,
                                                            stabilizer_list[i],
                                                            [],
                                                            meas_errors,
                                                            Is_after_two,
                                                            True)
            circ_bare = Encoded_Gate('S%i_bare'%i, [circ_bare]).circuit_wrap()
            
            stabs_circ.join_circuit(circ_flag)
            stabs_circ.join_circuit(circ_bare)

        return stabs_circ



    @classmethod
    def generate_whole_QEC_circ(cls, n_rep, stabilizer_list, flags_list,
                                meas_errors, Is_after_two, n_data,
                                initial_I, group_reps=False, high_indet=False):
        '''
        '''

        # We assume that the first half of the stabilizers are one kind (X)
        # and the second half are the other kind.

        n_stabs = len(stabilizer_list)
        stabs1, flags1 = stabilizer_list[:n_stabs/2], flags_list[:n_stabs/2]
        stabs2, flags2 = stabilizer_list[n_stabs/2:], flags_list[n_stabs/2:]
    
        if high_indet:
            circ_function = Flag_Correct.generate_high_indeterminacy_circuit
        else:
            circ_function = Flag_Correct.generate_all_flagged_stabs


        complete_circ = Circuit()
        for rep_i in range(n_rep):
            local_initial_I = False
            if rep_i == 0 and initial_I:  local_initial_I = True

            circ1 = circ_function(stabs1, flags1, meas_errors, Is_after_two, n_data,
                                  local_initial_I)
            circ2 = circ_function(stabs2, flags2, meas_errors, Is_after_two, n_data)
            circ1 = Encoded_Gate('Stabilizers_X%i'%rep_i, [circ1]).circuit_wrap()
            circ2 = Encoded_Gate('Stabilizers_Z%i'%rep_i, [circ2]).circuit_wrap()
            circ1.join_circuit(circ2)
            
            if group_reps:
                circ1 = Encoded_Gate('Rep_%i'%rep_i, [circ1]).circuit_wrap()
            complete_circ.join_circuit(circ1)

        return complete_circ

    
    
    @classmethod
    def generate_one_flagged_stab_ion(cls, i_first_anc, stabilizer, 
                                      meas_errors, initial_I=False,
                                      dephasing_during_MS=True,
                                      reordering_after=0):
        '''
        Specific to the d3 color code
        '''

        #n_flags = len(flags)
        n_data = 7
        n_flags = 1
        n_total = n_data + 1 + n_flags
        Pauli_oper = stabilizer[0][0]
        #coupling_gate = 'C'+Pauli_oper
        
        # The Is are just place-holders; no errors will be inserted after them.
        if Pauli_oper == 'X':  extra_gates = ['I','I']
        elif Pauli_oper == 'Z':  extra_gates = ['RY +', 'RY -']

        stab_circ = Circuit()
        if initial_I:
            for i in range(n_data):
                stab_circ.add_gate_at([i], 'Ism')


        stab_circ.add_gate_at([i_first_anc], 'PrepareZPlus')
        for i in range(n_flags):
            stab_circ.add_gate_at([i_first_anc+i+1], 'PrepareZPlus')


        # shuttling d1 to the operating zone and cooling 
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Icool')


        # couple d1 to s
        d1 = stabilizer[0][1]
        stab_circ.add_gate_at([d1], extra_gates[0])
        stab_circ.add_gate_at([d1], 'RX -')
        inter_qs = [d1, i_first_anc]
        stab_circ.add_gate_at(inter_qs, 'MS')
        stab_circ.add_gate_at([d1], extra_gates[1])
        idle_is = [i for i in range(n_total) if i not in inter_qs]
        if dephasing_during_MS:
            for i in idle_is:
                stab_circ.add_gate_at([i], 'Ism')

        # shuttling d1 back to the storage zone S2
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Icool')
        
        # couple s to f
        stab_circ.add_gate_at([i_first_anc,i_first_anc+1], 'MS')
        if dephasing_during_MS:
            for i in range(i_first_anc):
                stab_circ.add_gate_at([i], 'Ism')
        
        # shuttling d2 to the operating zone and cooling 
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Icool')

        # couple d2 to s
        d2 = stabilizer[1][1]
        stab_circ.add_gate_at([d2], extra_gates[0])
        stab_circ.add_gate_at([d2], 'RX -')
        inter_qs = [d2, i_first_anc]
        stab_circ.add_gate_at(inter_qs, 'MS')
        stab_circ.add_gate_at([d2], extra_gates[1])
        idle_is = [i for i in range(n_total) if i not in inter_qs]

        if dephasing_during_MS:
            for i in idle_is:
                stab_circ.add_gate_at([i], 'Ism')
        
        # shuttling d3 to the operating zone and cooling 
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Icool')

        # couple d3 to s
        d3 = stabilizer[2][1]
        stab_circ.add_gate_at([d3], extra_gates[0])
        stab_circ.add_gate_at([d3], 'RX -')
        inter_qs = [d3, i_first_anc]
        stab_circ.add_gate_at(inter_qs, 'MS')
        stab_circ.add_gate_at([d3], extra_gates[1])
        idle_is = [i for i in range(n_total) if i not in inter_qs]
        if dephasing_during_MS:
            for i in idle_is:
                stab_circ.add_gate_at([i], 'Ism')

        # shuttling d3 back to the storage zone S2
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Icool')

        # couple s to f
        stab_circ.add_gate_at([i_first_anc,i_first_anc+1], 'MS')
        if dephasing_during_MS:
            for i in range(i_first_anc):
                stab_circ.add_gate_at([i], 'Ism')
        
        # shuttling d4 to the operating zone and cooling 
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Icool')

        # couple d4 to s
        d4 = stabilizer[3][1]
        stab_circ.add_gate_at([d4], extra_gates[0])
        stab_circ.add_gate_at([d4], 'RX -')
        inter_qs = [d4, i_first_anc]
        stab_circ.add_gate_at(inter_qs, 'MS')
        stab_circ.add_gate_at([d4], extra_gates[1])
        idle_is = [i for i in range(n_total) if i not in inter_qs]
        if dephasing_during_MS:
            for i in idle_is:
                stab_circ.add_gate_at([i], 'Ism')
        
        # shuttling d4 back to the storage zone S2
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')

        if meas_errors:
            stab_circ.add_gate_at([i_first_anc], 'ImZ')
            stab_circ.add_gate_at([i_first_anc+1], 'ImZ')

        stab_circ.add_gate_at([i_first_anc], 'MeasureZ')
        stab_circ.add_gate_at([i_first_anc+1], 'MeasureZ')
        for i in range(i_first_anc):
            stab_circ.add_gate_at([i], 'Icool')
        
        # reordering after
        if reordering_after == 0:  delay = 0
        elif reordering_after == 1:  delay = 4
        elif reordering_after == 2:  delay = 5
        elif reordering_after == 3:  delay = 6

        for i in range(i_first_anc):
            for n_Is in range(delay):
                stab_circ.add_gate_at([i], 'Ism')

        # convert the last two qubits to ancilla
        stab_circ.to_ancilla([i_first_anc, i_first_anc+1])

        return stab_circ


    
    @classmethod
    def generate_one_nonFT_stab_ion(cls, i_first_anc, stabilizer, 
                                    meas_errors, initial_I=False,
                                    dephasing_during_MS=True,
                                    reordering_after=0):
        '''
        Uses a 5-qubit MS gate
        '''
        
        #n_flags = len(flags)
        n_data = 7
        n_total = n_data + 1
        Pauli_oper = stabilizer[0][0]
        coupling_gate = 'C'+Pauli_oper

        stab_circ = Circuit()
        if initial_I:
            for i in range(n_data):
                stab_circ.add_gate_at([i], 'Ism')
        
        stab_circ.add_gate_at([i_first_anc], 'PrepareXMinus')

        # shuttling d1, d2, d3, d4 to the operating zone and cooling 
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')
            stab_circ.add_gate_at([i], 'Icool')

        # performing the 5-qubit MS gate
        # Since I don't know how to compile the 5-qubit MS gate in
        # terms of C, H, and P, we just perform 4 CXs or CZs and
        # then add a gate name IMS5 to add the errors.
        stab_qs = [oper[1] for oper in stabilizer]
        for i in stab_qs:
            stab_circ.add_gate_at([i_first_anc, i], coupling_gate)
       
        stab_qs += [i_first_anc]
        idle_is = [i for i in range(n_total) if i not in stab_qs]
        if dephasing_during_MS:
            for i in idle_is:
                stab_circ.add_gate_at([i], 'Ism')
    
        stab_circ.add_gate_at(stab_qs, 'IMS5')

        # shuttling d1, d2, d3, d4 to the operating zone and cooling
        for i in range(n_total):
            stab_circ.add_gate_at([i], 'Ism')

        if meas_errors:
            stab_circ.add_gate_at([i_first_anc], 'ImX')
        stab_circ.add_gate_at([i_first_anc], 'MeasureX')

        # reordering after
        if reordering_after == 0:  delay = 0
        elif reordering_after == 1:  delay = 4
        elif reordering_after == 2:  delay = 5
        elif reordering_after == 3:  delay = 6

        for i in range(i_first_anc):
            for n_Is in range(delay):
                stab_circ.add_gate_at([i], 'Ism')

        # convert the last qubit to ancilla
        stab_circ.to_ancilla([i_first_anc])

        return stab_circ
        
        
    
    @classmethod
    def generate_whole_QEC_d3_ion(cls, stabilizers, meas_errors,
                                  initial_I=True,
                                  dephasing_during_MS=True,
                                  decoding='old'):
        '''
        decoding:  'new' refers to the new idea by Markus:
                   as soon as we detect 1 error, we measure
                   all the stabilizers 1 time non-FTly and
                   keep the last syndromes.

        '''
        # this is the case if the stabilizers are going to be
        # alternating: Sx1, Sz1, Sx2, Sz2, Sx3, Sz3
        if stabilizers[0][0][0] != stabilizers[1][0][0]:
            reordering_values = [0,1,0,2,0,3]
        # this is the case if the stabilizers are going to be
        # non-alternating: Sx1, Sx2, Sx3, Sz1, Sz2, Sz3
        else:
            reordering_values = [1,2,3,1,2,3]
        
        n_rep = 2
        FT_func = Flag_Correct.generate_one_flagged_stab_ion 
        nonFT_func = Flag_Correct.generate_one_nonFT_stab_ion 
        complete_circ = Circuit()
       
        #print 'Im in correction.py'
        #print reordering_values

        if decoding == 'new':
            
            # First round of stabilizers: only FT
            # only works for now if order of stabs is alternating
            for i_stab in range(len(stabilizers)):
                local_initial_I = False
                if i_stab == 0 and initial_I:  local_initial_I = True
            
                # the FT construction with 1 flag
                circ1 = FT_func(7, stabilizers[i_stab],
                                meas_errors,
                                local_initial_I,
                                dephasing_during_MS,
                                0)
                circ1 = Encoded_Gate('Stab_FT%i'%i_stab, [circ1]).circuit_wrap()
                
                if i_stab%2 == 1:
                    if i_stab == 1:  delay = 4
                    elif i_stab == 3:  delay = 5
                    elif i_stab == 5:  delay = 6
                    circ_shut = Circuit()
                    for i in range(7):
                        for n_Is in range(delay):
                            circ_shut.add_gate_at([i], 'Ism')
                    gate_name = 'Shuttling_%i_to_%i'%(i_stab/2, i_stab/2+1)
                    circ_shut = Encoded_Gate(gate_name, [circ_shut]).circuit_wrap()
                    circ1.join_circuit(circ_shut)

                complete_circ.join_circuit(circ1)

            # Second round of stabilizers: only non-FT
            # only works for now if order of stabs is alternating
            for i_stab in range(len(stabilizers)):
                    
                # the non-FT construction with no flags
                # and 2 5-qubit MS gates
                circ2 = nonFT_func(7, stabilizers[i_stab],
                                   meas_errors,
                                   False,
                                   dephasing_during_MS,
                                   0)
                circ2 = Encoded_Gate('Stab_nonFT%i'%i_stab, [circ2]).circuit_wrap()
                
                if i_stab%2 == 1:
                    if i_stab == 1:  delay = 4
                    elif i_stab == 3:  delay = 5
                    elif i_stab == 5:  delay = 6
                    circ_shut = Circuit()
                    for i in range(7):
                        for n_Is in range(delay):
                            circ_shut.add_gate_at([i], 'Ism')
                    gate_name = 'Shuttling_%i_to_%i'%(i_stab/2, i_stab/2+1)
                    circ_shut = Encoded_Gate(gate_name, [circ_shut]).circuit_wrap()
                    circ2.join_circuit(circ_shut)
                
                complete_circ.join_circuit(circ2)

            return complete_circ
                


        # First round of stabilizers: FT and non-FT
        for i_stab in range(len(stabilizers)):
            local_initial_I = False
            if i_stab == 0 and initial_I:  local_initial_I = True
                    
            # the FT construction with 1 flag
            circ1 = FT_func(7, stabilizers[i_stab],
                            meas_errors,
                            local_initial_I,
                            dephasing_during_MS,
                            reordering_values[i_stab])
            circ1 = Encoded_Gate('Stab_FT%i'%i_stab, [circ1]).circuit_wrap()
                    
            # the non-FT construction with no flags
            # and 2 5-qubit MS gates
            circ2 = nonFT_func(7, stabilizers[i_stab],
                               meas_errors,
                               False,
                               dephasing_during_MS,
                               reordering_values[i_stab])
            circ2 = Encoded_Gate('Stab_nonFT%i'%i_stab, [circ2]).circuit_wrap()
            circ1.join_circuit(circ2)
            complete_circ.join_circuit(circ1)

        
        # Second round of stabilizers: only non-FT
        for i_stab in range(len(stabilizers)):
                    
            # the non-FT construction with no flags
            # and 2 5-qubit MS gates
            circ2 = nonFT_func(7, stabilizers[i_stab],
                               meas_errors,
                               False,
                               dephasing_during_MS,
                               reordering_values[i_stab])
            circ2 = Encoded_Gate('Stab_nonFT%i'%i_stab, [circ2]).circuit_wrap()
            complete_circ.join_circuit(circ2)

        return complete_circ


    @classmethod
    def QEC_split_qubits(cls, Pauli_type='X', log_qubit=0, last_round=True):
        '''
        Measurement of stabilizers on a given logical qubit to perform
        splitting after merging during lattice surgery.
        We do one round of 3 stabilizers with flags + one round 
        without flags (+ one round of the other Pauli type to 
        catch possible w-2 hook errors)
        log_qubit: 0 (ctrl); 1 (targ); 2 (anc)
        '''
        
        n_code = 7
        n_total = 7*3
        stabs = [[3,4,5,6], [1,2,5,6], [0,2,4,6]]

        if Pauli_type == 'X':
            ent_gate = 'CX'
            other_gate = 'CZ'
        elif Pauli_type == 'Z':
            ent_gate = 'CZ'
            other_gate = 'CX'

        QEC_circ = Circuit()
        for i_stab in range(len(stabs)):
            stab_circ = Circuit()
            stab_circ.add_gate_at([n_total], 'PrepareXPlus')
            stab_circ.add_gate_at([n_total+1], 'PrepareZPlus')
            for i in range(len(stabs[i_stab])):
                q_index = stabs[i_stab][i]
                data_index = log_qubit*n_code + q_index
                stab_circ.add_gate_at([n_total, data_index], ent_gate)
                if i==0 or i==2:
                    stab_circ.add_gate_at([n_total, n_total+1], 'CX')
            stab_circ.add_gate_at([n_total], 'ImX')
            stab_circ.add_gate_at([n_total], 'MeasureX')
            stab_circ.add_gate_at([n_total+1], 'ImZ')
            stab_circ.add_gate_at([n_total+1], 'MeasureZ')
        
            stab_circ.to_ancilla(range(n_total, n_total+2))
            stab_circ = Encoded_Gate('Stab_FT%i'%i_stab, [stab_circ]).circuit_wrap()
            QEC_circ.join_circuit(stab_circ)

        for i_stab in range(len(stabs)):
            stab_circ = Circuit()
            stab_circ.add_gate_at([n_total], 'PrepareXPlus')
            for q_index in stabs[i_stab]:
                data_index = log_qubit*n_code + q_index
                stab_circ.add_gate_at([n_total, data_index], ent_gate)
            stab_circ.add_gate_at([n_total], 'ImX')
            stab_circ.add_gate_at([n_total], 'MeasureX')
            
            stab_circ.to_ancilla([n_total])
            stab_circ = Encoded_Gate('Stab_nonFT%i'%i_stab, [stab_circ]).circuit_wrap()
            QEC_circ.join_circuit(stab_circ)

        # One round of the other stabilizers.  These are only run if a flag was triggered 
        # during the measurement of the first stabilizers
        if last_round:
            for i_stab in range(len(stabs)):
                stab_circ = Circuit()
                stab_circ.add_gate_at([n_total], 'PrepareXPlus')
                for q_index in stabs[i_stab]:
                    data_index = log_qubit*n_code + q_index
                    stab_circ.add_gate_at([n_total, data_index], other_gate)
                stab_circ.add_gate_at([n_total], 'ImX')
                stab_circ.add_gate_at([n_total], 'MeasureX')
            
                stab_circ.to_ancilla([n_total])
                stab_circ = Encoded_Gate('Stab_nonFT%i'%i_stab, [stab_circ]).circuit_wrap()
                QEC_circ.join_circuit(stab_circ)
        

        QEC_circ = Encoded_Gate('EC%s_total'%Pauli_type, [QEC_circ]).circuit_wrap()
       
        return QEC_circ



    @classmethod
    def joint_QEC_split_qubits(cls, Pauli_type='X', log_qubits=[0,1]):
        '''
        '''
        QEC_circ = Flag_Correct.QEC_split_qubits(Pauli_type, log_qubits[0])        
        # When we measure the X stabilizers after the second merging, we don't
        # need to worry about w-2 hook X errors on the ancilla logical qubit
        # because we will measure it in the X basis.
        last_round = True
        if Pauli_type == 'X':
            last_round = False
        QEC_circ2 = Flag_Correct.QEC_split_qubits(Pauli_type, log_qubits[1], last_round)
        QEC_circ.join_circuit(QEC_circ2)
        QEC_circ = Encoded_Gate('JointQEC%s_flags'%Pauli_type, [QEC_circ]).circuit_wrap()

        return QEC_circ



    @classmethod
    def QEC_FT_lattsurg(cls, Pauli_type='X', log_qubit=0):
        '''
        This is the first QEC performed during the measurement of the XX operator.
        log_qubit: 0 (ctrl); 1 (targ); 2 (anc)
        '''
        
        n_code = 7
        n_total = 7*3

        # we always measure the non-boundary stabilizer first
        if Pauli_type == 'Z':
            stabs = [[1,2,5,6], [0,2,4,6], [3,4,5,6]]

        elif Pauli_type == 'X':
            #stabs = [[0,1,2,3], [1,2,4,5], [2,6,5,3]]   # Alejandro
            stabs = [[0,2,4,6], [3,4,5,6], [1,2,5,6]]    # Mauricio

        QEC_circ = Circuit()

        # First the FT version with flags
        for i_stab in range(len(stabs)):
            stab_circ = Circuit()
            stab_circ.add_gate_at([n_total], 'PrepareXPlus')
            stab_circ.add_gate_at([n_total+1], 'PrepareZPlus')
            for i in range(len(stabs[i_stab])):
                q_index = stabs[i_stab][i]
                data_index = log_qubit*n_code + q_index
                stab_circ.add_gate_at([n_total, data_index], 'C%s'%Pauli_type)
                if i==0 or i==2:
                    stab_circ.add_gate_at([n_total, n_total+1], 'CX')
            stab_circ.add_gate_at([n_total], 'ImX')
            stab_circ.add_gate_at([n_total], 'MeasureX')
            stab_circ.add_gate_at([n_total+1], 'ImZ')
            stab_circ.add_gate_at([n_total+1], 'MeasureZ')
        
            stab_circ.to_ancilla(range(n_total, n_total+2))
            stab_circ = Encoded_Gate('Stab_FT%i'%i_stab, [stab_circ]).circuit_wrap()
            QEC_circ.join_circuit(stab_circ)

        # then the non-FT version
        for i_stab in range(len(stabs)):
            stab_circ = Circuit()
            stab_circ.add_gate_at([n_total], 'PrepareXPlus')
            for q_index in stabs[i_stab]:
                data_index = log_qubit*n_code + q_index
                stab_circ.add_gate_at([n_total, data_index], 'C%s'%Pauli_type)
            stab_circ.add_gate_at([n_total], 'ImX')
            stab_circ.add_gate_at([n_total], 'MeasureX')
            
            stab_circ.to_ancilla([n_total])
            stab_circ = Encoded_Gate('Stab_nonFT%i'%i_stab, [stab_circ]).circuit_wrap()
            QEC_circ.join_circuit(stab_circ)

        QEC_circ = Encoded_Gate('QEC%s_FT'%Pauli_type, [QEC_circ]).circuit_wrap()

        return QEC_circ



    @classmethod
    def QECX_nonFT_lattsurg(cls, log_qubit=0):
        '''
        This is the second QEC performed during the measurement of the XX operator.
        '''
        
        n_code = 7
        n_total = 7*3
        #stabs = [[0,1,2,3], [1,2,4,5], [2,6,5,3]]   # Alejandro
        stabs = [[0,2,4,6], [3,4,5,6], [1,2,5,6]]    # Mauricio

        QEC_circ = Circuit()
        for i_stab in range(1,len(stabs)):
            stab_circ = Circuit()
            stab_circ.add_gate_at([n_total], 'PrepareXPlus')
            for q_index in stabs[i_stab]:
                data_index = log_qubit*n_code + q_index
                stab_circ.add_gate_at([n_total, data_index], 'CX')
            stab_circ.add_gate_at([n_total], 'ImX')
            stab_circ.add_gate_at([n_total], 'MeasureX')
            
            stab_circ.to_ancilla([n_total])
            stab_circ = Encoded_Gate('Stab_nonFT%i'%i_stab, [stab_circ]).circuit_wrap()
            QEC_circ.join_circuit(stab_circ)

        QEC_circ = Encoded_Gate('QECX_nonFT', [QEC_circ]).circuit_wrap()
        
        return QEC_circ


    
    @classmethod
    def measure_XXlogical(QEC_before=False):
        

        n_code = 7
        n_total = 7*3
        XX_circ = Circuit()

        if QEC_before:
            pass
        
        # measure the w-2 operator
        #qubits = [[1,6], [2,6]]   # Alejandro
        qubits = [[1,1], [2,1]]   # Mauricio
        X2_circ = Circuit()
        X2_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in qubits:
            X2_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], 'CX')
        X2_circ.add_gate_at([n_total], 'ImX')
        X2_circ.add_gate_at([n_total], 'MeasureX')
        X2_circ.to_ancilla([n_total])  
        X2_circ = Encoded_Gate('X2', [X2_circ]).circuit_wrap()
        XX_circ.join_circuit(X2_circ)

        # measure the w-4 operator
        #qubits = [[1,4], [2,4], [1,5], [2,5]]   # Alejandro
        qubits = [[1,3], [2,3], [1,5], [2,5]]   # Mauricio
        X4_circ = Circuit()
        X4_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in qubits:
            X4_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], 'CX')
        X4_circ.add_gate_at([n_total], 'ImX')
        X4_circ.add_gate_at([n_total], 'MeasureX')
        X4_circ.to_ancilla([n_total])  
        X4_circ = Encoded_Gate('X4', [X4_circ]).circuit_wrap()
        XX_circ.join_circuit(X4_circ)
        
        # measure the w-2 operator a second time
        #qubits = [[1,6], [2,6]]   # Alejandro
        qubits = [[1,1], [2,1]]   # Mauricio
        X2_circ = Circuit()
        X2_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in qubits:
            X2_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], 'CX')
        X2_circ.add_gate_at([n_total], 'ImX')
        X2_circ.add_gate_at([n_total], 'MeasureX')
        X2_circ.to_ancilla([n_total])  
        X2_circ = Encoded_Gate('X2', [X2_circ]).circuit_wrap()
        XX_circ.join_circuit(X2_circ)

        # measure the w-4 operator a second time
        #qubits = [[1,4], [2,4], [1,5], [2,5]]   # Alejandro
        qubits = [[1,3], [2,3], [1,5], [2,5]]   # Mauricio
        X4_circ = Circuit()
        X4_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in qubits:
            X4_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], 'CX')
        X4_circ.add_gate_at([n_total], 'ImX')
        X4_circ.add_gate_at([n_total], 'MeasureX')
        X4_circ.to_ancilla([n_total])  
        X4_circ = Encoded_Gate('X4', [X4_circ]).circuit_wrap()
        XX_circ.join_circuit(X4_circ)
        
        # FT measurement of X stabs on target (1) and ancilla (2)
        XX_circ.join_circuit(Flag_Correct.QEC_FT_lattsurg(Pauli_type, 1))
        XX_circ.join_circuit(Flag_Correct.QEC_FT_lattsurg(Pauli_type, 2))

        # nonFT measurement of X stabs on target (1) and ancilla (2)
        #XX_circ.join_circuit(Flag_Correct.QECX_nonFT_lattsurg(1))
        #XX_circ.join_circuit(Flag_Correct.QECX_nonFT_lattsurg(2))

        # measure the w-2 operator
        #qubits = [[1,6], [2,6]]   # Alejandro
        qubits = [[1,1], [2,1]]   # Mauricio
        X2_circ = Circuit()
        X2_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in qubits:
            X2_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], 'CX')
        X2_circ.add_gate_at([n_total], 'ImX')
        X2_circ.add_gate_at([n_total], 'MeasureX')
        X2_circ.to_ancilla([n_total])  
        X2_circ = Encoded_Gate('X2', [X2_circ]).circuit_wrap()
        XX_circ.join_circuit(X2_circ)

        # measure the w-4 operator
        #qubits = [[1,4], [2,4], [1,5], [2,5]]   # Alejandro
        qubits = [[1,3], [2,3], [1,5], [2,5]]   # Mauricio
        X4_circ = Circuit()
        X4_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in qubits:
            X4_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], 'CX')
        X4_circ.add_gate_at([n_total], 'ImX')
        X4_circ.add_gate_at([n_total], 'MeasureX')
        X4_circ.to_ancilla([n_total])  
        X4_circ = Encoded_Gate('X4', [X4_circ]).circuit_wrap()
        XX_circ.join_circuit(X4_circ)
   
        XX_circ = Encoded_Gate('MeasureXX_flags', [XX_circ]).circuit_wrap()

        return XX_circ
    
    
    
    @classmethod
    def measure_logical_boundary(cls, Pauli_meas='X', QEC_before=False):
        '''
        Circuit to FTly measure a logical XX or ZZ on the boundary of two
        logical qubits encoded on the d-3 color code.
        '''

        n_code = 7
        n_total = 7*3
        total_circ = Circuit()

        if QEC_before:
            pass
        
        if Pauli_meas == 'Z':
            w2_qubits = [[0,3],[2,3]]
            w4_qubits = [[0,0], [2,0], [0,4], [2,4]]
            ent_gate = 'CZ'
            nonanc_qubit = 0

        elif Pauli_meas == 'X':
            #qubits = [[1,6], [2,6]]   # Alejandro
            w2_qubits = [[1,1], [2,1]]   # Mauricio
            #qubits = [[1,4], [2,4], [1,5], [2,5]]   # Alejandro
            w4_qubits = [[1,3], [2,3], [1,5], [2,5]]   # Mauricio
            ent_gate = 'CX'
            nonanc_qubit = 1

        # measure the w-2 operator
        w2_circ = Circuit()
        w2_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in w2_qubits:
            w2_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], ent_gate)
        w2_circ.add_gate_at([n_total], 'ImX')
        w2_circ.add_gate_at([n_total], 'MeasureX')
        w2_circ.to_ancilla([n_total])  
        w2_circ = Encoded_Gate('%s2'%Pauli_meas, [w2_circ]).circuit_wrap()
        total_circ.join_circuit(w2_circ)

        # measure the w-4 operator
        w4_circ = Circuit()
        w4_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in w4_qubits:
            w4_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], ent_gate)
        w4_circ.add_gate_at([n_total], 'ImX')
        w4_circ.add_gate_at([n_total], 'MeasureX')
        w4_circ.to_ancilla([n_total])  
        w4_circ = Encoded_Gate('%s4'%Pauli_meas, [w4_circ]).circuit_wrap()
        total_circ.join_circuit(w4_circ)
        
        # measure the w-2 operator a second time
        w2_circ = Circuit()
        w2_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in w2_qubits:
            w2_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], ent_gate)
        w2_circ.add_gate_at([n_total], 'ImX')
        w2_circ.add_gate_at([n_total], 'MeasureX')
        w2_circ.to_ancilla([n_total])  
        w2_circ = Encoded_Gate('%s2'%Pauli_meas, [w2_circ]).circuit_wrap()
        total_circ.join_circuit(w2_circ)

        # measure the w-4 operator a second time
        w4_circ = Circuit()
        w4_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in w4_qubits:
            w4_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], ent_gate)
        w4_circ.add_gate_at([n_total], 'ImX')
        w4_circ.add_gate_at([n_total], 'MeasureX')
        w4_circ.to_ancilla([n_total])  
        w4_circ = Encoded_Gate('%s4'%Pauli_meas, [w4_circ]).circuit_wrap()
        total_circ.join_circuit(w4_circ)
        
        # FT measurement of X stabs on target (1) and ancilla (2)
        total_circ.join_circuit(Flag_Correct.QEC_FT_lattsurg(Pauli_meas, nonanc_qubit))
        total_circ.join_circuit(Flag_Correct.QEC_FT_lattsurg(Pauli_meas, 2))

        # nonFT measurement of X stabs on target (1) and ancilla (2)
        #XX_circ.join_circuit(Flag_Correct.QECX_nonFT_lattsurg(1))
        #XX_circ.join_circuit(Flag_Correct.QECX_nonFT_lattsurg(2))

        # measure the w-2 operator a third time
        w2_circ = Circuit()
        w2_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in w2_qubits:
            w2_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], ent_gate)
        w2_circ.add_gate_at([n_total], 'ImX')
        w2_circ.add_gate_at([n_total], 'MeasureX')
        w2_circ.to_ancilla([n_total])  
        w2_circ = Encoded_Gate('%s2'%Pauli_meas, [w2_circ]).circuit_wrap()
        total_circ.join_circuit(w2_circ)

        # measure the w-4 operator a third time
        w4_circ = Circuit()
        w4_circ.add_gate_at([n_total], 'PrepareXPlus')
        for pair in w4_qubits:
            w4_circ.add_gate_at([n_total,n_code*pair[0]+pair[1]], ent_gate)
        w4_circ.add_gate_at([n_total], 'ImX')
        w4_circ.add_gate_at([n_total], 'MeasureX')
        w4_circ.to_ancilla([n_total])  
        w4_circ = Encoded_Gate('%s4'%Pauli_meas, [w4_circ]).circuit_wrap()
        total_circ.join_circuit(w4_circ)
 
        total_circ = Encoded_Gate('Measure%s%s_flags'%(Pauli_meas,Pauli_meas), 
                                  [total_circ]).circuit_wrap()

        return total_circ



    @classmethod
    def latt_surg_CNOT(cls, initial_I):
        '''
        '''
        n_code = 7
        CNOT_circ = Circuit()
        
        if initial_I:
            I_circuit = Circuit()
            for i in range(3*n_code):
                I_circuit.add_gate_at([i], 'I')
            I_circuit = Encoded_Gate('Logical_I', [I_circuit]).circuit_wrap()
            #CNOT_circ.join_circuit_at(range(3*n_code), I_circuit)
            CNOT_circ.join_circuit(I_circuit)

        measureXX_circ = Flag_Correct.measure_XXlogical()
        CNOT_circ.join_circuit(measureXX_circ)

        return CNOT_circ



class Bare_Correct:
    '''
    Measure stabilizers with a bare ancilla for each stabilizer, no cat states.
    Class defined to implement Cross's [[7,1,3]] code.
    '''
    
    @classmethod
    def generate_bare_meas(cls, n_data, stabilizer_list, input_I_round=False, 
                           meas_errors=True, Is_after_two=False,
                           parallel=True, inputcirc=False):
        '''
        Generates stabilizer measurements with bare ancilla, no cat states.
        parallel = True means that we don't recicle the ancilla

        stabilizer_list is assumed to have the following format:
            [
              [('X',0), ('X',4)], ... 

            ]

        '''
        n_ancilla = len(stabilizer_list)
        bare_meas_circ = Circuit()
        if input_I_round:
            for i in range(n_data):
                bare_meas_circ.add_gate_at([i], 'I')

        if parallel:
            for i in range(n_data, n_data+n_ancilla):
	            bare_meas_circ.add_gate_at([i], 'PrepareXPlus')
        
            for i in range(len(stabilizer_list)):
                gateprefix = 'C'
                for gate in stabilizer_list[i]:
                    bare_meas_circ.add_gate_at([n_data+i,gate[1]], gateprefix+gate[0])
                    if Is_after_two:
                        bare_meas_circ.add_gate_at([gate[1]], 'I')
                        bare_meas_circ.add_gate_at([n_data+i], 'I')
        

            for i in range(n_data, n_data+n_ancilla):
                if meas_errors:
                    bare_meas_circ.add_gate_at([i], 'ImX')
                bare_meas_circ.add_gate_at([i], 'MeasureX')
                
            bare_meas_circ.to_ancilla(range(n_data, n_data+n_ancilla))

        else:
            pass

        return bare_meas_circ


    @classmethod
    def generate_rep_bare_meas(cls, n_data, stabilizer_list, n_rounds, input_I_round, 
                               meas_errors, Is_after_two_qubit, initial_trans=False,
                               ancilla_parallel=False, CSS=False):
        '''
        Is_after_two_qubits:  whether or not we want to add Is after 2-qubit gates,
                              to add errors on the Cross code.
        CSS:  if True, we divide the stabilizers into two groups which we measure
              separately.
        '''
        n = n_data
        n_ancilla = len(stabilizer_list)
        s_l = stabilizer_list[:]
        i_I = input_I_round
        m_e = meas_errors
        i_t = Is_after_two_qubit
        bare_meas_circ = Bare_Correct.generate_bare_meas
        rep_meas_circ = Circuit()


        # New part I added to run the surface-17 code MGA: 07/30/17.
        if CSS:
            n_stab = len(s_l)/2
            for i in range(n_rounds):
                if i == 0: i_I_local = i_I
                else:      i_I_local = False
                gate_nameX = 'SX_%i'%i
                stab_circX = bare_meas_circ(n, s_l[:n_stab], i_I_local, m_e, i_t)
                stab_circX = Encoded_Gate(gate_nameX, [stab_circX]).circuit_wrap()
                gate_nameZ = 'SZ_%i'%i
                stab_circZ = bare_meas_circ(n, s_l[n_stab:], False, m_e, i_t)
                stab_circZ = Encoded_Gate(gate_nameZ, [stab_circZ]).circuit_wrap()
                rep_meas_circ.join_circuit(stab_circX, ancilla_parallel)
                rep_meas_circ.join_circuit(stab_circZ, ancilla_parallel)

            return rep_meas_circ
        #######################################################


        if initial_trans != False:
            for i in range(n_data):
                rep_meas_circ.add_gate_at([i], initial_trans)

        for i in range(n_rounds):
            gate_name = 'Bare_Measurement_'+str(i+1)
            if i==0:
                stab_circ = Encoded_Gate(gate_name,[bare_meas_circ(n,s_l,i_I,m_e,i_t)]).circuit_wrap()
            else:
                stab_circ = Encoded_Gate(gate_name,[bare_meas_circ(n,s_l,False,m_e,i_t)]).circuit_wrap()
                
            rep_meas_circ.join_circuit(stab_circ, ancilla_parallel)
        
        rep_meas_circ = Encoded_Gate('EC_BareCorrect',
                        [rep_meas_circ]).circuit_wrap()

        #brow.from_circuit(rep_meas_circ, True)
        #sys.exit(0)

        return rep_meas_circ



    @classmethod
    def generate_rep_bare_meas_CSS(cls, n_data, stabilizer_list, n_rounds, initial_I, meas_errors,
                                   Is_after_twoq=False, initial_trans=False, anc_parallel=False):
        '''
        '''
        



class Cat_Correct:
    '''
    Methods that employ Shor's ancilla to perform EC.
    ''' 

    @classmethod
    def cat_syndrome_4_old(cls, stabilizer_list, redundancy=1, 
            verify=False, ancilla_parallel=True):
        '''
        This method is not used anymore.
        '''
        
        #set method syndrome extraction
        if verify==True:
            cat_circ=Cat_Correct.create_4_cat_verify
            corr_circ=Cat_Correct.create_4_cat_verify_correction
        else:
            cat_circ=Cat_Correct.create_4_cat_no_verify
            corr_circ=Cat_Correct.create_4_cat_no_verify_correction
            
        #important constants
        n_data = len(stabilizer_list[0])
        n_stab = len(stabilizer_list)

        cat_synd_circ=Circuit()
        cat_corr_circ=Circuit()
        
        for i_stab,stabilizer in enumerate(stabilizer_list):
        
            # Redundancy is needed for fault tolearance of syndrome
            # measurments.
            
            redundancy_circ=Circuit()

            non_identity_index=[]
            for i_pauli,pauli in enumerate(stabilizer):
                if (pauli!='I'):
                    non_identity_index.append(i_pauli)
            
            for i_redund in range(redundancy):
                stab_circ=cat_circ(stabilizer)
                redundancy_circ.join_circuit_at(non_identity_index, 
                                stab_circ)
                
            cat_synd_circ.join_circuit(redundancy_circ, ancilla_parallel)

            cat_synd_circ = Encoded_Gate('EC_Cat4Syndrome',
                        [cat_synd_circ]).circuit_wrap()
            
            #cat_corr_circ = corr_circ(stabilizer)
            
            #cat_synd_circ.join_circuit_at(non_identity_index, cat_corr_circ)

        return cat_synd_circ
    


    @classmethod
    def cat_syndrome_4_test(cls, stabilizer_list, redundancy=1, 
            verify=False, ancilla_parallel=True,
            diVincenzo=True, initial_I=False):
        '''
        Method to be employed temporarily while we change all the 
        functions that used to call the cat_syndrome_4_test.
        '''

        return Cat_Correct.cat_syndrome_4(stabilizer_list, redundancy,
                                        verify, ancilla_parallel,
                                        diVincenzo, initial_I)    


    @classmethod
    def cat_syndrome_4(cls, stabilizer_list, redundancy=1, 
                        verify=False, ancilla_parallel=True,
                        diVincenzo=True, initial_I=False,
                        initial_trans=False, code='Steane',
                        meas_errors=True, Is_after_two=False):
        '''
        Method to do Shor's ancilla EC.
        The method is designed for a code with weight-4
        stabilizers.
        Right now, we have only implemented the Steane
        and 5-qubit codes.   
        It can be generalized to other codes.
        MGA 6/16/2016.
        '''     

        #set method syndrome extraction
        if verify==True:
            cat_circ=Cat_Correct.create_4_cat_verify
            corr_circ=Cat_Correct.create_4_cat_verify_correction
        else:
            if diVincenzo == True:
                cat_circ = Cat_Correct.create_4_cat_diVincenzo
            else:
                cat_circ=Cat_Correct.create_4_cat_no_verify
            corr_circ=Cat_Correct.create_4_cat_no_verify_correction
            
        #important constants
        n_data = len(stabilizer_list[0])
        n_stab = len(stabilizer_list)
        
        cat_synd_circ=Circuit()
        if initial_trans != False:
            for i in range(n_data):
                cat_synd_circ.add_gate_at([i], initial_trans)

        if code == 'Steane':

            if n_stab < 6:
                stabs = stabilizer_list[:]
                for i in range(redundancy):
                    stab_circ = Cat_Correct.create_4_cat_stabilizer(n_data, initial_I, 
                                                                    stabs, cat_circ, 
                                                                    'stabs_Steane_', 
                                                                    i+1, meas_errors,
                                                                    Is_after_two)
                    cat_synd_circ.join_circuit(stab_circ, ancilla_parallel)
                

            else:
                X_stabs = stabilizer_list[:3]
                Z_stabs = stabilizer_list[3:]
                for i in range(redundancy):
                    stab_circX = Cat_Correct.create_4_cat_stabilizer(n_data, initial_I, 
                                                                     X_stabs, cat_circ, 
                                                                    'X_stabs_Steane_', 
                                                                     i+1, meas_errors,
                                                                     Is_after_two)
                    cat_synd_circ.join_circuit(stab_circX, ancilla_parallel)
                    initial_I = False
                    stab_circZ = Cat_Correct.create_4_cat_stabilizer(n_data, initial_I, 
                                                                     Z_stabs, cat_circ, 
                                                                    'Z_stabs_Steane_', 
                                                                     i+1, meas_errors,
                                                                     Is_after_two)
                    cat_synd_circ.join_circuit(stab_circZ, ancilla_parallel)


        elif code == '5qubit':
            for i in range(redundancy):
                stab_circ = Cat_Correct.create_4_cat_stabilizer(n_data, initial_I, 
                                                                stabilizer_list[:], 
                                                                cat_circ, 
                                                               'stabs_5qubit_', 
                                                                i+1, meas_errors,
                                                                Is_after_two)
                cat_synd_circ.join_circuit(stab_circ, ancilla_parallel)
                initial_I = False
       
 
        cat_synd_circ = Encoded_Gate('EC_CatCorrect',
                        [cat_synd_circ]).circuit_wrap()

        return cat_synd_circ



    @classmethod
    def create_4_cat_stabilizer(cls, n_q, initial_I, stabs, cat_circ_method, 
				                gate_name, rep, meas_errors, Is_after_two):
        '''
        
        '''
        total_circ = Circuit()
        if initial_I:
            for i in range(n_q):
                total_circ.add_gate_at([i], 'I')
        for i_stab, stab in enumerate(stabs):
            non_iden_index = []
            for i_pauli, pauli in enumerate(stab):
                if pauli != 'I':
                    non_iden_index.append(i_pauli)
            stab_circ = cat_circ_method(stab, meas_errors, Is_after_two)
            total_circ.join_circuit_at(non_iden_index, stab_circ)            

        full_gate_name = gate_name + str(rep)
        total_circ = Encoded_Gate(full_gate_name, [total_circ]).circuit_wrap()
        
        return total_circ


    
    @classmethod
    def create_stabs_d5_color_code(cls, initial_I, stabs, meas_errors, Is_after_two,
                                   rep_verif_weight8, rep):

        '''
        stabs are of the form [('X',8), ('X',3), ...]
        The weight-8 stab must always come first
        '''
        
        if len(stabs[0]) != 8:
            raise TypeError('The first stabilizer must always be the weight-8 one.')

        n_q = 17
        total_circ = Circuit()
        if initial_I:
            for i in range(17):
                total_circ.add_gate_at([i], 'I')
        #    for i in range(17,17+8):
        #        total_circ.add_gate_at([i], 'I_noerror')
        
        stab_kind = stabs[0][0][0]  # X or Z
        stabs_indices = [[oper[1] for oper in stabilizer] for stabilizer in stabs]

        # First we add the octagon stabilizer
        octagon_circ = Cat_Correct.measure_weight8_stab(rep_verif_weight8,
                                                        meas_errors, 
                                                        Is_after_two,
                                                        stab_kind)
        
        octagon_range = stabs_indices[0][:] + range(n_q,n_q+8)
        total_circ.join_circuit_at(octagon_range, octagon_circ)

        # Next we add the other seven
        stabs_orig = []
        for stabilizer in stabs_indices[1:]:
            stabs_orig += [[stab_kind if i in stabilizer else 'I' for i in range(n_q)]]
    
        gate_name = 'Weight4_stabs' + stab_kind
        cat_circ_method = Cat_Correct.create_4_cat_diVincenzo
        rest_circ = Cat_Correct.create_4_cat_stabilizer(n_q, False, stabs_orig, 
                                                        cat_circ_method, gate_name, rep, 
                                                        meas_errors, Is_after_two)

        total_circ.join_circuit(rest_circ)
        total_circ = Encoded_Gate('Stabs%s_%i'%(stab_kind,rep), [total_circ]).circuit_wrap()

        return total_circ



    @classmethod
    def EC_d5_color_code(cls, initial_I, stabs, meas_errors, Is_after_two, 
                         rep_verif_w8, redun):
        '''
        '''
        
        total_circ = Circuit()

        stabs_X = stabs[:8]
        stabs_Z = stabs[8:]

        for rep in range(redun):
            stab_circ_X = Cat_Correct.create_stabs_d5_color_code(initial_I, 
                                                                 stabs_X,
                                                                 meas_errors,
                                                                 Is_after_two,
                                                                 rep_verif_w8,
                                                                 rep+1)
            total_circ.join_circuit(stab_circ_X, False)
           
            # after first round we set initial_I to False
            initial_I = False
            stab_circ_Z = Cat_Correct.create_stabs_d5_color_code(initial_I, 
                                                                 stabs_Z,
                                                                 meas_errors,
                                                                 Is_after_two,
                                                                 rep_verif_w8,
                                                                 rep+1)
            total_circ.join_circuit(stab_circ_Z, False)

        total_circ = Encoded_Gate('QEC_d5', [total_circ]).circuit_wrap()

        return total_circ



    @classmethod
    def create_4_cat_verify(cls,stabilizer=['X','X','X','X']):
        """A helper method for creating a 4 qubit cat state and verify it.
        """
        cat_circuit = Circuit()
        for index in range(4,9):
            if index!=5:
                cat_circuit.add_gate_at([index],'PrepareZPlus')
            else:
                cat_circuit.add_gate_at([index],'PrepareXPlus')
        cat_circuit.add_gate_at([5,6],'CX')
        cat_circuit.add_gate_at([5,4],'CX')
        cat_circuit.add_gate_at([6,7],'CX')
        cat_circuit.add_gate_at([4,8],'CX')
        cat_circuit.add_gate_at([7,8],'CX')
        verify_gate = cat_circuit.add_gate_at([8],'MeasureZDestroy')
        
        final_gate = Verify_Gate(gate_name='4_cat_with_verify',
                                circuit_list=[cat_circuit])

        #create circuit with this final gate
        final_circ = Circuit()
        final_circ.add_gate_at(range(4,9),final_gate)
        
        non_identity_iter = 0
                
        for i_pauli,pauli in enumerate(stabilizer):
            if (pauli!='I'):
                gate_name = 'C' + pauli
                i1 = non_identity_iter+4
                final_circ.add_gate_at([i1+4, i1],gate_name)        
                non_identity_iter += 1

        if non_identity_iter>4:
            print ('Error!! stabilizer doesn\'t have four non-identity'+
                ' pauli operators.') 

        m_gates_cat=[]
        for i in range(4,8):
            m_gates_cat+=[final_circ.add_gate_at([i],'MeasureXDestroy')]
        
        final_circ.to_ancilla(range(4,9))

        return final_circ

    
    @classmethod    
    def create_4_cat_verify_correction(cls,stabilizer=['X','X','X','X']):
        #specify correction circuit
        corr_circ = Circuit()
        for i_pauli,pauli in enumerate(stabilizer):
            if (pauli!='I'):
                corr_circ.add_gate_at([i_pauli],pauli)
            
        #make correction gate
            
        corr_gate = Correction_Gate(gate_name='Cat_correction',
                circuit_list=[corr_circ])

        return corr_gate.circuit_wrap()
        
            
    @classmethod
    def create_4_cat_no_verify(cls,stabilizer):
        """A helper method for creating a 4 qubit cat state without a
        verification step. 
        
        This is specific for the Steane code, where the stabilizers contain
        either 4 X's or 4 Z's.  This could easily be generalized to other CSS
        codes.
        """     
        cat_circuit = Circuit()
        for index in range(4,8):
            if index==5:
                cat_circuit.add_gate_at([index],'PrepareXPlus')
            else:
                cat_circuit.add_gate_at([index],'PrepareZPlus')
        cat_circuit.add_gate_at([5,6],'CX')
        cat_circuit.add_gate_at([5,4],'CX')
        cat_circuit.add_gate_at([6,7],'CX')
                
        if stabilizer.count('X')==0:
            """Z stabilizers"""
            for index in range(4):
                cat_circuit.add_gate_at([index+4], 'H')
                cat_circuit.add_gate_at([index,index+4], 'CX')
        
        
        elif stabilizer.count('Z')==0:
            """X stabilizers"""
            for index in range(4):
                cat_circuit.add_gate_at([index+4,index], 'CX')

        m_gates_cat=[]
        

        if stabilizer.count('X')==0:
            """Z stabilizers"""
            for index in range(4,8):
                m_gates_cat += [cat_circuit.add_gate_at([index], 'MeasureZ')]

        elif stabilizer.count('Z')==0:
            """X stabilizers"""
            for index in range(4,8):
                m_gates_cat += [cat_circuit.add_gate_at([index], 'MeasureX')]


        cat_circuit.to_ancilla(range(4,8))

        return cat_circuit

    

    @classmethod
    def create_4_cat_diVincenzo(cls,stabilizer, meas_errors=True,
                                Is_after_two=False):
        """A helper method for creating a 4 qubit cat state without a
        verification step. 
        This is specific for a code where the stabilizers have weight 4.
        So far it's used for the Steane and the 5-qubit codes, but it
        is easily generalizable.
        if meas_errors == True: the measurements are faulty.
        if Is_after_two == True: we add I gates after 2-qubit gates.
                                 this is used on the ion_trap_simple
                                 error model.
        """ 
        cat_circuit = Circuit()
        for index in range(4,8):
            if index==5:
                cat_circuit.add_gate_at([index],'PrepareXPlus')
            else:
                cat_circuit.add_gate_at([index],'PrepareZPlus')
        cat_circuit.add_gate_at([5,6],'CX')
        if Is_after_two:
            cat_circuit.add_gate_at([5], 'I')
            cat_circuit.add_gate_at([6], 'I')
        cat_circuit.add_gate_at([5,4],'CX')
        if Is_after_two:
            cat_circuit.add_gate_at([5], 'I')
            cat_circuit.add_gate_at([4], 'I')
        cat_circuit.add_gate_at([6,7],'CX')
        if Is_after_two:
            cat_circuit.add_gate_at([6], 'I')
            cat_circuit.add_gate_at([7], 'I')
               
        #if stabilizer.count('X')==0:
        #    """Z stabilizers"""
        #    for index in range(4):
                #cat_circuit.add_gate_at([index+4], 'H')
                #cat_circuit.add_gate_at([index,index+4], 'CX')
        #        cat_circuit.add_gate_at([index+4,index], 'CZ')
            
        #elif stabilizer.count('Z')==0:
        #    """X stabilizers"""
        #    for index in range(4):
        #        cat_circuit.add_gate_at([index+4,index], 'CX')

        # More concise way to do the same thing: MGA 6/16/2016.
        index = 0
        for pauli in stabilizer:
            if pauli != 'I':
                cat_circuit.add_gate_at([index+4,index], 'C'+pauli)
                if Is_after_two:
                    cat_circuit.add_gate_at([index+4], 'I')
                    cat_circuit.add_gate_at([index], 'I')
                index += 1

        cat_circuit.add_gate_at([5,7], 'CX')
        if Is_after_two:
            cat_circuit.add_gate_at([5], 'I')
            cat_circuit.add_gate_at([7], 'I')
        cat_circuit.add_gate_at([6,4], 'CX')
        if Is_after_two:
            cat_circuit.add_gate_at([6], 'I')
            cat_circuit.add_gate_at([4], 'I')
        cat_circuit.add_gate_at([5,6], 'CX')    
        if Is_after_two:
            cat_circuit.add_gate_at([5], 'I')
            cat_circuit.add_gate_at([6], 'I')
        
        m_gates_cat=[]
        for index in range(4,8):
            if index == 5:
                if meas_errors:
                    cat_circuit.add_gate_at([index], 'ImX')
                m_gates_cat += [cat_circuit.add_gate_at([index], 'MeasureX')]
            else:
                if meas_errors:
                    cat_circuit.add_gate_at([index], 'ImZ')
                m_gates_cat += [cat_circuit.add_gate_at([index], 'MeasureZ')]   
            
        cat_circuit.to_ancilla(range(4,8))

        return cat_circuit

    
    
    @classmethod
    def create_8_cat_verify(cls, meas_errors=True,
                            Is_after_two=False):
        '''
        creates an 8-qubit cat state and verifies it twice.
        The verification involves measuring the parity of the first
        and eighth qubits.  
        This is used to measure the weight-8 stabilizer in the 4.8.8
        color code.  
        We have to verify twice to detect w-2 errors that propagate
        to form uncorrectable w-3 errors.
        if meas_errors == True: the measurements are faulty.
        if Is_after_two == True: we add I gates after 2-qubit gates.
                                 this is used on the ion_trap_simple
                                 error model.
        ''' 
        cat_circuit = Circuit()
        for index in range(10):
            if index==3:
                cat_circuit.add_gate_at([index],'PrepareXPlus')
            else:
                cat_circuit.add_gate_at([index],'PrepareZPlus')
        list_CNOTs = [[3,4], [4,5], [5,6], [6,7], [3,2], [2,1], [1,0],
                      [0,8], [7,8], [0,9], [7,9]]
        
        for CNOT in list_CNOTs:
            cat_circuit.add_gate_at(CNOT,'CX')
            if Is_after_two:
                cat_circuit.add_gate_at([CNOT[0]], 'I')
                cat_circuit.add_gate_at([CNOT[1]], 'I')
        
        if meas_errors:
            cat_circuit.add_gate_at([8], 'ImZ')
            cat_circuit.add_gate_at([9], 'ImZ')
        cat_circuit.add_gate_at([8], 'MeasureZ')
        cat_circuit.add_gate_at([9], 'MeasureZ')

        cat_circuit.replace_qubit_ids(range(8,18))

        #cat_circuit.to_ancilla([8,9])
        cat_circuit = Encoded_Gate('prep_cat_8', [cat_circuit]).circuit_wrap()

        return cat_circuit
    
    

    @classmethod
    def repeat_prep_8_cat(cls, rep, meas_errors, Is_after_two):
        '''
        Repeats the preparation and verification of an 8-qubit cat state
        "rep" times.  
        This is a somewhat artificial way for the fast sampler to handle 
        the unpredictibility of the number of times we'll have to prepare
        the 8-qubit cat state.
        '''
    
        rep_cat_circ = Circuit()
        for i in range(rep):
            cat_circuit = Cat_Correct.create_8_cat_verify(meas_errors, Is_after_two)
            rep_cat_circ.join_circuit(cat_circuit)        
        
        rep_cat_circ.to_ancilla([16,17])
        rep_cat_circ = Encoded_Gate('rep_cat_8_%i'%rep, [rep_cat_circ]).circuit_wrap()

        return rep_cat_circ
    


    @classmethod
    def measure_weight8_stab(cls, rep, meas_errors, Is_after_two, stab='X'):
        '''
        '''

        weight8_circ = Cat_Correct.repeat_prep_8_cat(rep, meas_errors, Is_after_two)
        
        coupling_circ = Circuit()
        for i in range(8):
            coupling_circ.add_gate_at([i+8, i], 'C'+stab)
            if Is_after_two:
                coupling_circ.add_gate_at([i+8], 'I')
                coupling_circ.add_gate_at([i], 'I')

        for i in range(8,16):
            if meas_errors:
                coupling_circ.add_gate_at([i], 'ImX')
            coupling_circ.add_gate_at([i], 'MeasureX')

        coupling_circ = Encoded_Gate('coupling_to_data', [coupling_circ]).circuit_wrap()
        weight8_circ.join_circuit(coupling_circ)
        weight8_circ = Encoded_Gate('weight8_S'+stab, [weight8_circ]).circuit_wrap()
        
        return weight8_circ

    

    @classmethod    
    def create_4_cat_no_verify_correction(cls,stabilizer=['X','X','X','X']):
        #specify correction circuit
        corr_circ = Circuit()
        for i_pauli,pauli in enumerate(stabilizer):
            if (pauli!='I'):
                corr_circ.add_gate_at([i_pauli],pauli)
            
        #make correction gate
            
        corr_gate = Correction_Gate(gate_name='Cat_correction',
                            circuit_list=[corr_circ])
        
        return corr_gate.circuit_wrap()
        
    @classmethod
    def create_cat_3(cls):
        """The great thing about three qubit cat states is that you don't have 
        to  verify them!--Dave Bacon"""
        cat_circuit = Circuit()
        for index in [0,2]:
            cat_circuit.add_gate_at([index],'PrepareZPlus')
        cat_circuit.add_gate_at([1],'PrepareXPlus')
        cat_circuit.add_gate_at([1,0],'CX')
        cat_circuit.add_gate_at([1,2],'CX')
        return cat_circuit
    
    @classmethod    
    def create_cat_n(cls, n):
        """Create a cat state on n qubits where n>3.
        """
        
        if n<=3:
            return 'Use a CNOT gate for n=2 or create_cat_3 for n=3.'
        
        circ = Circuit()
        for i in range(n):
            if not i == int(n/2):
                circ.add_gate_at([i],'PrepareZPlus')
            else:
                circ.add_gate_at([i],'PrepareXPlus')
                
        for i in range(n):
            if i < int(n/2):
                circ.add_gate_at([i+1,i],'CX')
            if i >= int(n/2):
                circ.add_gate_at([i,i+1],'CX')
                
        circ.add_gate_at([0,n],'CX')
        circ.add_gate_at([n-1,n],'CX')
        
        circ.add_gate_at([n],'MeasureZDestroy')

        #if MeasureZ gives -1, prepare again.

        gate = Prepare_Gate(gate_name=str(n)+' cat with verify', 
                            circuit_list = [circ])
        
        return gate.circuit_wrap()
        
#---------------#
#---------------#
        
class Steane_Correct:
    
    @classmethod
    def steane_syndrome(cls,ecc, redundancy =1, FT=True, kind='GFI'):
        """Steane syndrome measurement requires the creation of encoded zero and plus states in that code."""
        
        n_data = ecc.Code().block_size #number of data qubits in the code  
        n_ancilla = n_data #this code uses the same number of ancillas as data qubits
        # total data + ancilla qubits
        
        ###
        #Z syndrome
        ###
        
        #Create FT Encoded Zero
        synd_z_circ = Circuit()
        z_plus_circ, nothing = Steane_Correct.FT_encoded_zero_Steane()
        #brow.from_circuit(z_plus_circ, True)
        synd_z_circ.join_circuit_start_id(n_data,z_plus_circ)
        
        #kick phase down and measure    
        for i in range(n_data): 
            synd_z_circ.add_gate_at([n_data+i,i],'CX')
            synd_z_circ.add_gate_at([n_data+i],'MeasureX')
        
        synd_z_circ.to_ancilla(range(n_data, 2*n_data))
        synd_z_circ = Encoded_Gate('Steane_Syndrome_Z', [synd_z_circ]).circuit_wrap()
        
        #brow.from_circuit(synd_z_circ, True)
                
        ###
        #X syndrome
        ###
        
        #Create FT Encoded Plus
        synd_x_circ = Circuit()
        x_plus_circ, nothing = Steane_Correct.FT_encoded_plus_Steane()
        synd_x_circ.join_circuit_start_id(n_data,x_plus_circ)
        
        #Kick phase down and measure    
        for i in range(n_data): 
            synd_x_circ.add_gate_at([i,n_data+i],'CX')
            synd_x_circ.add_gate_at([n_data+i],'MeasureZ')
        
        synd_x_circ.to_ancilla(range(n_data,2*n_data))
        synd_x_circ = Encoded_Gate('Steane_Syndrome_X', [synd_x_circ]).circuit_wrap()
        
        #brow.from_circuit(synd_x_circ, True)
        
        main_circ = Circuit()
        main_circ.join_circuit(synd_z_circ)
        main_circ.join_circuit(synd_x_circ)
        
        return Container_Gate('EC_Steane_Correct',[main_circ]).circuit_wrap()

    @classmethod
    def encoded_zero_Steane(cls):
        """
        Creates an encoded zero for exclusive use in the 
        Steane[7,1,3] code
        """
        enc_zero = Circuit()
        for i in [0,1,3]:
            enc_zero.add_gate_at([i], 'PrepareXPlus')
        for i in [2,4,5,6]:
            enc_zero.add_gate_at([i], 'PrepareZPlus')
        cnots = [[0,2],[3,5],[1,6],[0,4],[3,6],[1,5],[0,6],[1,2],[3,4]]
        for i in cnots:
            enc_zero.add_gate_at(i, 'CX')
        return enc_zero



    @classmethod
    def encoded_plus_Steane(cls):
        """
        Creates an encoded plus for exclusive use in the 
        Steane[7,1,3] code.
        The Hadamard gate is transversal in the Steane code, so |+> is  
        almost identical to |0>; the only difference is the 
        CNOT's are flipped, |0>s are |+>s, and |+>s are |0>s.
        """
        enc_plus = Circuit()
        for i in [0,1,3]:
            enc_plus.add_gate_at([i], 'PrepareZPlus')
        for i in [2,4,5,6]:
            enc_plus.add_gate_at([i], 'PrepareXPlus')
        cnots = [[2,0],[5,3],[6,1],[4,0],[6,3],[5,1],[6,0],[2,1],[4,3]]
        for i in cnots:
            enc_plus.add_gate_at(i, 'CX')
        return enc_plus


    
    @classmethod
    def encoded_plus_i_Steane(cls):
        """
        Creates an encoded +i for exclusive use in the 
        Steane[7,1,3] code.
        It first creates a plus state and then applies a logical
        S gate.  S_L = S^t on each qubit.
        """
        enc_plus_i = cls.encoded_plus_Steane()
        for i in range(7):
            enc_plus_i.add_gate_at([i], 'Z')
            enc_plus_i.add_gate_at([i], 'P')
        return enc_plus_i

    

    @classmethod
    def FT_encoded_zero_Steane(cls, set_ancilla=True):
        circ1 = cls.encoded_zero_Steane()
        circ2 = cls.encoded_zero_Steane()
        circ1.join_circuit_at(range(7,14), circ2)

        meas_gates = []
        for i in range(7):
            circ1.add_gate_at([i,7+i], 'CX')
        for i in range(7):
            meas_gates += [circ1.add_gate_at([7+i], 'MeasureZDestroy')]
        
        if set_ancilla:  circ1.to_ancilla(range(7,14))
        
        return circ1, meas_gates

    
    
    @classmethod
    def FT_encoded_plus_Steane(cls, set_ancilla=True):
        circ1 = cls.encoded_plus_Steane()
        circ2 = cls.encoded_plus_Steane()
        circ1.join_circuit_at(range(7,14), circ2)

        #circ1.add_gate_at([0], 'Z')
        #circ1.add_gate_at([1], 'Z')
        meas_gates = []
        for i in range(7):
            circ1.add_gate_at([7+i,i], 'CX')
        for i in range(7):
            meas_gates += [circ1.add_gate_at([7+i], 'MeasureXDestroy')]
        
        if set_ancilla:  circ1.to_ancilla(range(7,14))
        
        return circ1, meas_gates
        


    @classmethod
    def detect_errors(cls, error_to_correct, n_data, initial_trans=False):
        '''
        very simple circuit to detect Z or X errors.
        It assumes the ancilla input is a FT logical |0> or |+> state
        '''
        circ = Circuit()
        m_gates = []
        # let's check this works:
        #if error_to_correct == 'X':
        #   circ.add_gate_at([0], 'X')
        #   circ.add_gate_at([1], 'X')

        if initial_trans != False:
            for i in range(n_data):
                circ.add_gate_at([i], initial_trans)

        if error_to_correct == 'X':
            for i in range(n_data):
                circ.add_gate_at([i, i+n_data], 'CX')
                m_gates += [circ.add_gate_at([i+n_data], 'MeasureZ')]
        elif error_to_correct == 'Z':
            for i in range(n_data):
                circ.add_gate_at([i+n_data, i], 'CX')
                m_gates += [circ.add_gate_at([i+n_data], 'MeasureX')]
        
        circ.to_ancilla(range(n_data, 2*n_data))

        return circ, m_gates


    @classmethod
    def steane_syndrome(cls, ancilla_parallel=False):
        """ 
        Steane syndrome measurement requires the creation of encoded zero and plus 
        states in that code.
        """
        
        #n_data = ecc.Code().block_size # number of data qubits in the code  
        
        n_data = 7                      # specific for the Steane code
        n_ancilla = n_data              # this code uses the same number of 
                                        # ancillas as data qubits

        ###
        #Z syndrome
        ###

        # create FT Encoded Zero
        synd_z_circ = Circuit()
        z_plus = Steane_Correct.FT_encoded_zero_Steane(False)[0]
        z_plus_circ = Encoded_Gate('Prepare_logical_Zero', [z_plus]).circuit_wrap()
        synd_z_circ.join_circuit_start_id(n_data, z_plus_circ)

        # kick phase down and measure    
        for i in range(n_data):
            synd_z_circ.add_gate_at([n_data+i,i],'CX')
            synd_z_circ.add_gate_at([n_data+i],'MeasureX')

        #synd_z_circ.to_ancilla(range(n_data, 2*n_data))
        synd_z_circ = Encoded_Gate('Steane_Syndrome_Z', [synd_z_circ]).circuit_wrap()

        ###
        #X syndrome
        ###

        # Create FT Encoded Plus
        synd_x_circ = Circuit()
        x_plus = Steane_Correct.FT_encoded_plus_Steane(False)[0]
        if ancilla_parallel:
            n_initial = n_data + 2*n_ancilla
        else:
            n_initial = n_data
        
        x_plus_circ = Encoded_Gate('Prepare_logical_Plus', [x_plus]).circuit_wrap()
        synd_x_circ.join_circuit_start_id(n_initial, x_plus_circ)

        # kick phase down and measure    
        for i in range(n_data):
            synd_x_circ.add_gate_at([i,n_initial+i],'CX')
            synd_x_circ.add_gate_at([n_initial+i],'MeasureZ')

        #synd_x_circ.to_ancilla(range(n_data, 2*n_data))
        synd_x_circ = Encoded_Gate('Steane_Syndrome_X', [synd_x_circ]).circuit_wrap()

        main_circ = Circuit()
        main_circ.join_circuit(synd_z_circ)
        main_circ.join_circuit(synd_x_circ)

        main_circ.to_ancilla(range(n_data, n_initial + 2*n_ancilla))

        steane_synd_circ = Encoded_Gate('EC_SteaneCorrect',[main_circ]).circuit_wrap()

        return steane_synd_circ


    @classmethod
    def steane_syndrome_old(cls, ecc, redundancy=1):
        """
        Steane syndrome measurement requires the creation of encoded zero
        and plus states. Currently we just assume that encoding is the Steane
        code. This can be expanded later.
        
        """
        
        n_data = ecc.Code().block_size  #number of data qubits in the code  
        n_anc = n_data                  # this code uses the same number of ancillas as  
                                        # total data + ancilla qubits will be 7+4*7=35
        
        main_circ=Circuit()
        
        
        ###
        #Z syndrome
        ###
        for i_redund in range(redundancy):

            # Create Encoded Zero and Verify
            enc_zero_circ, m_gates_enc_zero = cls.FT_encoded_zero_Steane()
            #gate = Verify_Gate(gate_name='Encoded Zero With Verify',
            #                   circuit_list=[enc_zero_circ])
            gate = Encoded_Gate(gate_name='FT_Encoded_Zero',
                                circuit_list=[enc_zero_circ])            

            synd_circ = gate.circuit_wrap()
            
            #(2) kick phase down and (3) measure Z
            
            synd_circ.replace_qubit_ids(range(n_data, n_data+2*n_anc))
            
            err_det_circ, m_gates_synd = cls.detect_errors('Z', n_data)
            synd_circ.join_circuit_at(range(n_data+n_anc), err_det_circ)
            synd_circ.to_ancilla(range(n_data,n_data+2*n_anc))
            main_circ.join_circuit_at(range(n_data), synd_circ)
            
        
        main_circ = Encoded_Gate('EC_SteaneSyndrome', [main_circ])
        
        corr_circ=Circuit()
        for i in range(7):
            corr_circ.add_gate_at([i],'Z')
        
        gate = Correction_Gate(gate_name='SteaneZ', circuit_list=[corr_circ])
        
        main_circ.join_circuit(gate.circuit_wrap())
                        
        ###
        #X syndrome
        ###
        for i_redund in range(redundancy):
            
            #Create Encoded Plus and Verify
            enc_plus_circ, m_gates_enc_plus = cls.FT_encoded_plus_Steane()
            gate = Verify_Gate(gate_name='Encoded Plus With Verify',
                    circuit_list=[enc_plus_circ])
            
            synd_circ = gate.circuit_wrap()
            
            #(2) kick phase down and (3) measure Z
            
            synd_circ.replace_qubit_ids(range(n_data, n_data+2*n_anc))
            
            err_det_circ, m_gates_synd = cls.detect_errors('X', n_data)
            synd_circ.join_circuit_at(range(n_data+n_anc), err_det_circ)
            synd_circ.to_ancilla(range(n_data,n_data+2*n_anc))
            main_circ.join_circuit_at(range(n_data), synd_circ)
            
            
        corr_circ=Circuit()
        for i in range(7):
            corr_circ.add_gate_at([i],'X')
                
        gate = Correction_Gate(gate_name='SteaneX', circuit_list=[corr_circ])

        main_circ.join_circuit(gate.circuit_wrap()) 
        
        return main_circ
        


class Knill_Correct:
    
    @classmethod
    def detect_errors(cls, n_data, order='normal', initial_trans=False):
        '''
        very simple circuit to detect errors.
        It assumes the input on the ancilla is a FT logical |+>
        and a FT logical |0>.
        It then prepares a logical Bell pair, applies the logical
        CNOT, and performs the logical X and Z measurements.
        
        if order == 'normal':     the circuit has data qubits on top and
                        logical Bell state at the bottom.
        elif order == 'reverse':  the other way around
        '''

        circ = Circuit()
        m_gates = []
        # let's check this works:
        #if error_to_correct == 'X':
        #   circ.add_gate_at([0], 'X')
        #   circ.add_gate_at([1], 'X')

        if order == 'normal':       
            # prepare the logical Bell pair
            for i in range(n_data):
                if initial_trans != False:
                    circ.add_gate_at([i], initial_trans)
                circ.add_gate_at([i+n_data, i+2*n_data], 'CX')
        
            # perform the logical CNOT and the measurements
            for i in range(n_data):
                circ.add_gate_at([i, i+n_data], 'CX')
                m_gates += [circ.add_gate_at([i], 'MeasureX')]
                m_gates += [circ.add_gate_at([i+n_data], 'MeasureZ')]

        elif order == 'reverse':
            # prepare the logical Bell pair
            for i in range(n_data):
                if initial_trans != False:
                    circ.add_gate_at([i+2*n_data], initial_trans)
                circ.add_gate_at([i+n_data, i], 'CX')
        
            # perform the logical CNOT and the measurements
            for i in range(n_data):
                circ.add_gate_at([i+2*n_data, i+n_data], 'CX')
                m_gates += [circ.add_gate_at([i+n_data], 'MeasureZ')]
                m_gates += [circ.add_gate_at([i+2*n_data], 'MeasureX')]

        circ.to_ancilla(range(n_data, 3*n_data))

        return circ, m_gates



    @classmethod
    def knill_syndrome(cls, ecc):
        '''
        '''
        n_data = 7
        n_ancilla = n_data
        bell_pair_circ = Circuit()
        
        # create FT encoded plus
        x_plus = Steane_Correct.FT_encoded_plus_Steane(False)[0]
        x_plus_circ = Encoded_Gate('Prepare_logical_Plus', [x_plus]).circuit_wrap()
        bell_pair_circ.join_circuit_start_id(n_data, x_plus_circ)
    
        # create FT encoded zero
        z_plus = Steane_Correct.FT_encoded_zero_Steane(False)[0]
        z_plus_circ = Encoded_Gate('Prepare_logical_Zero', [z_plus]).circuit_wrap()
        bell_pair_circ.join_circuit_start_id(n_data + 2*n_ancilla, z_plus_circ)
        
        # create Bell pair
        cnot = ecc.Generator.create_encoded_circuit('CX')
        cnot_circ = Encoded_Gate('logical_CX', [cnot]).circuit_wrap()
        cnot_qubits = range(n_data, n_data+n_ancilla) + \
                    range(n_data+2*n_ancilla, n_data+3*n_ancilla)
        bell_pair_circ.join_circuit_at(cnot_qubits, cnot_circ)
        
        bell_pair_circ = Encoded_Gate('Prepare_logical_Bell_pair', 
                                    [bell_pair_circ]).circuit_wrap()
        
    
        # couple to data and measure
        cnot = ecc.Generator.create_encoded_circuit('CX')
        cnot_circ = Encoded_Gate('logical_CX', [cnot]).circuit_wrap()
        cnot_qubits = range(n_data + n_ancilla)
        bell_pair_circ.join_circuit_at(cnot_qubits, cnot_circ)

        meas_X = ecc.Generator.create_encoded_circuit('MeasureXDestroy')
        meas_X_circ = Encoded_Gate('MeasureX', [meas_X]).circuit_wrap()
        meas_X_qubits = range(n_data)
        bell_pair_circ.join_circuit_at(meas_X_qubits, meas_X_circ)
        
        meas_Z = ecc.Generator.create_encoded_circuit('MeasureZDestroy')
        meas_Z_circ = Encoded_Gate('MeasureZ', [meas_X]).circuit_wrap()
        meas_Z_qubits = range(n_data, n_data + n_ancilla) 
        bell_pair_circ.join_circuit_at(meas_Z_qubits, meas_Z_circ)

        bell_pair_circ.to_ancilla(range(n_data, n_data + 4*n_ancilla))
        bell_pair_circ = Encoded_Gate('EC_KnillCorrect', [bell_pair_circ]).circuit_wrap()       

        return bell_pair_circ



class Aliferis_Cross_Correct:
    
    @classmethod
    def measure_T_operator(cls, data1, data2, ancilla, n=3, t='X'):
        """Measurement of two-qubit operator that results from the spliting of 
        a stabilizer. The label T comes from the original notation in Dave 
        Bacon's paper.
        """
        temp_anc = n**2
        T_meas = Circuit()
        if t == 'X':
            preparation = 'PrepareXPlus'
            control1, control2 = temp_anc, temp_anc
            target1, target2 = data1, data2
            measurement = 'MeasureXDestroy'
        else:
            preparation = 'PrepareZPlus'
            control1, control2 = data1, data2
            target1, target2 = temp_anc, temp_anc
            measurement = 'MeasureZDestroy'
        
        T_meas.add_gate_at([temp_anc], preparation)
        T_meas.add_gate_at([control1, target1], 'CX')
        T_meas.add_gate_at([control2, target2], 'CX')
        T_meas.add_gate_at([temp_anc], measurement)
        print T_meas.qubits()
        q = T_meas.qubits()[2]
        q.qubit_type = 'ancilla'
        q.qubit_id = ancilla
        T_meas.update_map()
                
        return T_meas


    @classmethod
    def measure_stabilizer(cls, first_data, second_data, first_ancilla, 
        n=3, t='X'):
        """
        Measurement of a stabilizer
        
        """
        Stabilizer_measurement = Circuit()

        for i in range(n):
            if t=='X':
                d1 = first_data*n + i
                d2 = second_data*n + i
            else:
                d1 = i*n + first_data
                d2 = i*n + second_data
            
            meas_op = Aliferis_Cross_Correct().measure_T_operator
            T_measure = meas_op(d1, d2, first_ancilla+i,n,t)
            T_gate = Encoded_Gate(gate_name='T_measure',
                                circuit_list=[T_measure])
            #T_gate = Correction_Gate(gate_name='T_measure', circuit_list=[T_measure])
            Stabilizer_measurement.add_gate(T_gate)
            
        return Stabilizer_measurement

    
    @classmethod
    def measure_all_stabilizers(cls, n=3):
        """
        Currently only implemented for the Bacon-Shor [9,1,3] code
        """
        
        All_stabilizers = Circuit()
        meas_stab = Aliferis_Cross_Correct().measure_stabilizer
        a = 0
        for t in ['X','Z']:
            for s in [[0,1],[0,2],[1,2]]:
                meas_stab = meas_stab(s[0],s[1],a,n,t)
                meas_stab_gate = Encoded_Gate(gate_name='stab_measurement',
                                            circuit_list=[meas_stab])
                All_stabilizers.add_gate(meas_stab_gate)
                a += n

        #for t in ['X','Z']:
        #   for i in range(n-1):
        #       meas_stab = Aliferis_Cross_Correct().measure_stabilizer(first=i, n=n, t=t)
        #       meas_stab_gate = Encoded_Gate(gate_name='stab_measurement', circuit_list=[meas_stab])
        #       All_stabilizers.add_gate(meas_stab_gate)
            
        #All_stabilizers.to_ancilla(range(n**2, n**2 + n))
        
        return All_stabilizers

        
        
