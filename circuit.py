'''
circuit.py
new trial version by Yu tomita (Sep 15 2011)
'''
import time
import copy
import operator

#from faultTolerant import decoder

class CircuitException(Exception): pass
class SizeMismatchError(CircuitException): pass 

class Qubit(object):
    '''
    '''
    
    def __init__(self, qubit_id=None, qubit_type='data', level=0):
        '''
        Qubits have qubit_id, qubit_type and level. 
        qubit_type can be 'data' or 'ancilla',
        level is its level of encoding (i.e. level 1 qubit is a logical qubit encoded once).
        
        This is semi-singleton (only one Qubit object per qubit_id is created).
        
        If your qubit_id is the same as an existing q.qubit_id, it returns q instead of a new Qubit object.
        '''
        self.qubit_id = qubit_id
        self.qubit_type = qubit_type
        self.level = level

    def __repr__(self):
        return ('<Qubit (qubit_id=' + repr(self.qubit_id) + ', qubit_type=' + repr(self.qubit_type) + 
                    ', level='+ repr(self.level)  +')>')

    def __eq__(self, other):
        return (self.qubit_id == other.qubit_id and
                self.qubit_type == other.qubit_type and
                self.level == other.level)
    
    def __hash__(self):
        return hash(self.qubit_id) ^ hash(self.qubit_type) ^ hash(self.level)

    def __cmp__(self, other):
        '''
        A < B (both Qubit obj) is true when 
            A.level < B.level, or,
            A.level == B.level and A.qubit_id < B.qubit_id
        When both qubit_id and level are the same between A and B, they are considered
        as the same object when compared (with different object id)
        '''
        if self.level < other.level: return -1
        if self.level > other.level: return 1
        if self.qubit_type != other.qubit_type:
            if self.qubit_type=='data':return -1
            if self.qubit_type=='ancilla':return 1
        if self.qubit_id < other.qubit_id:  return -1
        if self.qubit_id == other.qubit_id: return 0
        if self.qubit_id > other.qubit_id:  return 1

    def __hash__(self):
        '''
        Returns a hash value of this object. The value is xor of qubit_id and level.
        qubit_type is not considered as __cmp__ does not involve it.
        '''
        qubit_info = repr(self.qubit_id) + repr(self.qubit_type) +repr(self.level)
        return int(''.join([str(ord(s)) for s in qubit_info]))

class Gate(object):
    '''
    '''
    
    def __init__(self, gate_name=None, qubits=None, skip=None, duration=None):
        '''
        Gates have the name of the gate operation, the qubit being acted upon, and a boolean used in error correction.
        The gate can be a string of any Clifford operation with any modifiers.
        The qubits are a list of object instances from the Qubit class.
        '''
        self.gate_name = gate_name
        self.qubits = qubits
        self.skip = skip
        self.duration = duration
    
    #def __eq__(self, other):
    #   return (self.gate_name == other.gate_name and
    #           self.qubits == other.qubits and
    #           self.skip == other.skip)
        
    def set_sub_gates(self, gates):
        ''' Make sub_circuit of given gates. This removes the sub_circuit if already created.
        '''
        self.sub_circuit = Circuit()
        self.sub_circuit.gates = gates
        self.sub_circuit.update_map()

    def get_sub_gates(self):
        if not sub_circuit: return None
        return self.sub_circuit.gates

    def replace_qubit_ids(self, qubit_ids, qubit_type='data'):
        if qubit_type =='data':
            qubits = self.data_qubits()
        else:
            qubits = self.ancilla_qubits()
    
        for i,q in enumerate(qubits):
            q.qubit_id=qubit_ids[i] 
        
        return self
        
    def qubits_ordered(self, ascending=True):
        ''' Returns a list of qubits by ascending or descending order of their ids. 
        '''
        qubits = self.qubits
        qubits = sorted(qubits) if ascending else sorted(qubits)[::-1]
        return qubits
        
    def data_and_ancilla_qubits(self , ordered = True):
        '''Returns a tuple of list of data and ancilla qubits. i.e. ([data qubits], [ancilla qubits])
        '''
        all_qubits = self.qubits_ordered(ordered)
        ancilla = []
        data = []
        for q in all_qubits:
            if 'ancilla' in q.qubit_type:
                ancilla += [q]
            elif 'data' in q.qubit_type:
                data += [q]
            else:
                print 'qubit %r has qubit_type not data nor ancilla'%q
        return data, ancilla
        
    def data_qubits(self, ordered = True):
        return self.data_and_ancilla_qubits(ordered)[0]
        
    def ancilla_qubits(self, ordered = True):
        return self.data_and_ancilla_qubits(ordered)[1]
        
    def circuit_wrap(self):
        circ= Circuit()
        circ.add_gate(self)
        return circ
    
    def clean_qubits(self, qubits={}):
        
        if not qubits:
            qubits = {}
            for q in self.qubits:
                qubits[(q.qubit_id,q.qubit_type)] = q
        
        new_qubits=[]
        for i,q in enumerate(self.qubits):
            try:
                self.qubits[i] = qubits[(q.qubit_id,q.qubit_type)]
            except KeyError:
                qubits[(q.qubit_id,q.qubit_type)] = q
        if isinstance(self,Container_Gate):
            for circuit in self.circuit_list:
                qubits = circuit.clean_qubits(qubits)
            #self.update_qubits()
        if isinstance(self,Permute_Gate):

            for i,q in enumerate(self.output_qubits):
                try:
                    self.output_qubits[i] = qubits[(q.qubit_id,q.qubit_type)]
                except KeyError:
                    qubits[(q.qubit_id,q.qubit_type)] = q
    
    def unpack(self):
        ''' Returns all the gates in its sub-circuit, recursively all the way. 
        '''
        if isinstance(self, Encoded_Gate):
            if not self.circuit_list[0]:    return [self]
            gates = [x for g in self.circuit_list[0].gates for x in g.unpack()]
            if self.gate_name[0:2]=='EC':
                for g in gates:
                    g.skip = False
            return gates
        elif isinstance(self, Correction_Gate):
            self.skip= False
            return [self]
        elif isinstance(self, Oracle_Gate):     # Added this later on.  Mauricio
            return self.circuit_list[0].gates       
        else:
            return [self]

    def __repr__(self):
        return ('<Gate (gate_name=' + repr(self.gate_name) + ', qubits=' +repr(self.qubits)  +')>')
                
class Physical_Gate(Gate):
    def __init__(self, gate_name=None,qubits=None):
        self.start_time=None
        self.end_time=None
        self.is_error = False
        Gate.__init__(self, gate_name, qubits)
        
    def update_times(self, start_time, end_time):
        self.start_time = start_time
        self.end_time = end_time
        #print 'start time ',start_time, ' end time ',end_time
    
class Container_Gate(Gate):
    '''General category of gate that contain circuits. Used as an inhereting class for 
            encoded gate, verify gates, correction gates.'''
    def __init__(self, gate_name=None, circuit_list=None):
        #print gate_name
        #print circuit_list
        self.start_time=None  #should only be for physical gates
        self.end_time=None #should only be for physical gates
        qubits=[]
        for circ in circuit_list:
            qubits=list(set(qubits)|set(circ.qubits()))
            
        Gate.__init__(self, gate_name, qubits)
        self.circuit_list=[] if not circuit_list else circuit_list

    def update_times(self, start_time, end_time):
        self.start_time = start_time
        self.end_time = end_time
        #print 'start time ',start_time, ' end time ',end_time
        
    def update_qubits(self):
        qubits=[]
        for circ in self.circuit_list:
            qubits=list(set(qubits)|set(circ.qubits())) 
        self.qubits=qubits
        
    def qubits_ordered(self, ascending=True):
        self.update_qubits()
        return sorted(self.qubits)
    
    def __repr__(self):
        return ('<Gate (gate_name=' + repr(self.gate_name) + 
             ', qubits=' +repr(self.qubits)  +')>')
                
class Classically_Controlled_Gate(Container_Gate):
    def __init__(self, gate_name=None, classical_control=None, circuit_list=None):
        Container_Gate.__init__(self, gate_name, circuit_list)
        self.classical_control = classical_control
        self.target_and_syndrome_qubits = (self.classical_control.target_qubits,
                                           self.classical_control.syndrome_qubits)
        
    def decoder(self, function):
        '''
        Calls the decoder with the given function, returning the
        correction circuit.
        '''
        
        if function == 'parity':
            circ = Decoder.Decoder().parity(-1, self.target_and_syndrome_qubits[0], self.gate_name)
        elif function == 'level-by-level':
            print 'Not implemented yet!'
        else:
            raise ValueError('Error. Invalid function name')
            
        
class Encoded_Gate(Container_Gate):
    def __init__(self, gate_name=None, circuit_list=None):
                Container_Gate.__init__(self, gate_name, circuit_list)
        
        
class Verify_Gate(Classically_Controlled_Gate):
    def __init__(self, gate_name=None, classical_control=None, circuit_list=None):
        Classically_Controlled_Gate.__init__(self, gate_name, classical_control, circuit_list)
        
class Correction_Gate(Classically_Controlled_Gate):
    def __init__(self, gate_name=None, classical_control= None, circuit_list=None):
        Classically_Controlled_Gate.__init__(self, gate_name, classical_control, circuit_list)
        
class Permute_Gate(Gate):
    def __init__(self, gate_name=None, input_qubits=None, output_qubits=None):
        Gate.__init__(self, gate_name, input_qubits)
        self.output_qubits=output_qubits
        
class Oracle_Gate(Container_Gate):
    def __init__(self, gate_name=None, circuit_list=None):
        Container_Gate.__init__(self, gate_name, circuit_list)

        
class Classical_Control:
    '''
    
    '''
    def __init__(self, function, syndrome_qubits = None, target_qubits =None, controls = None):
        '''All these qubits are at a logical level. The physical level qubits need not be specified.
        Using the logical qubit tree list should be enough to relavant qubit syndrome values.
        
        measured_syndrome_qubits: qubits from which to extract the syndromes. We have to wait until
        the physical level qubits are measured.
        
        target_qubits: These are the qubits that are really being measured. So if we wish to correct
        errors on qubits, these are the qubits that are corrected. These are the highest level qubits.
        
        '''
        self.function=function
        self.syndrome_qubits = [] if not syndrome_qubits else syndrome_qubits
        self.target_qubits = [] if not target_qubits else target_qubits
        self.correction_delay_circuit = None
        self.controls = [] if not controls else controls

    def get_syndromes_for_stabilizers(target_qubit_list, qubit_pairing_list, stabilizer_list, syndrome_measurements):
        '''
        Given a list of target qubits, a list of (target_qubits, syndrome_qubits) pairs,
        a list of stabilizers, and a dictionary of syndrome measurements for the syndrome 
        qubits, produces a list of {-1, +1}^N, where N is the number of stabilizers, 
        that are the eigenvalues of each stabilizer (in the order that they were given).
        
        Argument List
        -------------
        target_qubit_list:
            List of target qubits (in order).
        qubit_pairing_list:
            List of tuples of tuples of the form ((q1, q2), (s1, s2)), where q1 and q2 are
            target qubits and s1 and s2 are syndrome qubits. The tuples give information
            about which syndrome qubits hold the measurements for which target qubits. It
            should be noted that each tuple of target qubits shares a tuple with several
            different tuples of syndrome qubits, and this is the result of a repetition
            measurement.
        stabilizer_list:
            List of lists of single character strings in {'I', 'X', 'Z'}. These are the
            stabilizers that we are finding the eigenvalues for (in order).
        syndrome_measurements:
            A dictionary of {Qubit : int} that tells us what the measurements for each
            syndrome qubit were. These values will be in {-1, 1}.
        '''
        
        stabilizer_eigenvalues = []
        
        # Go through each stabilizer, finding the corresponding eigenvalue.
        for stabilizer in stabilizer_list:
            
            # Find the target qubits associated with the stabilizer's
            # errors.
            target_qubits = []
            for i in range(0, len(stabilizer)):
                if stabilizer[i] != 'I':
                    target_qubits.append(target_qubits[i])
            target_qubits = tuple(target_qubits)
            
            # Find all tuples in qubit_pairing_list that contain
            # the tuple of target qubits that we're looking at,
            # and grab the corresponding syndrome qubit tuples.
            qubit_tuples = []
            for (t_tuple, s_tuple) in qubit_pairing_list:
                if t_tuple == target_qubits:
                    qubit_tuples.append(s_tuple)
            
            # For each index (i.e. corresponding to a target qubit), take
            # a majority vote over all of the syndrome tuples to find
            # the value for that target qubit.
            syndrome_qubit_values = []
            for i in range(0, len(qubit_tuples[0])):
                index_values = []
                for syndrome_tuple in qubit_tuples:
                    index_values.append(syndrome_tuple[i])

                # Take a majority vote over all of the values in index_values
                value_dict = {}
                value_dict[-1] = 0
                value_dict[1] = 0
                
                for value in index_values:
                    value_dict[value] += 1
                
                # Sort value_dict by its values, producing a list of tuples
                # of the form (key, value).
                sorted_value_dict = sorted(value_dict.iteritems(), key = operator.itemgetter(1))
                
                # The majority vote is the value of the last tuple in the 
                # sorted list.
                majority_vote = sorted_value_dict[-1][1]
                syndrome_qubit_values.append(majority_vote)
            
            # Now that we have the values for each target qubit that
            # has an error for this stabilizer, we can find the
            # eigenvalue for this stabilizer by combining these values
            # (i.e. finding the product of them).
            eigenvalue = reduce(operator.mul, syndrome_qubit_values, 1)
            stabilizer_eigenvalues.append(eigenvalue)
        
        return stabilizer_eigenvalues
                    
    
    def __repr__(self):
        return ('< Classical_Control (\n' + ')>\n')
             
            
class Circuit:
    '''
    '''
    def __init__(self, gates=None, qubits=None, level=0):
        '''
        Circuit contains a list of gates, and a dictionary of qubit:gates pair. 
        The list of gates needs to be in relative-time order. Note that this order is sufficient, but 
        not necessary.
        level is its level of encoding (i.e. level 1 qubit is a logical qubit encoded once).
        '''
        self.gates = [] if not gates else gates
        self.qubit_gates_map = {} if not qubits else dict.fromkeys(qubits)
        self.level = level

    #Functions to be regularly used by users
        
    def add_gate(self, gate):
        data_qubit_ids = [q.qubit_id for q in gate.data_qubits(False)]

        return self.add_end_gate(data_qubit_ids,gate)
        
    def add_gate_at(self, data_qubit_ids, gate):
        return self.add_end_gate(data_qubit_ids,gate)

    def join_circuit(self, circuit, ancilla_parallel=True):
        ''' Takes a circuit to extend and attaches it to self. data_qubits specify the where data 
        qubits get mapped. All other qubits
        '''
        self.end_join_circuit(circuit, ancilla_parallel)

    def join_circuit_at(self, data_qubit_ids, circuit):
        ''' Takes a circuit to extend and attaches it to self. data_qubits specify the where data 
        qubits get mapped. All other qubits
        '''
        circuit.replace_qubit_ids(data_qubit_ids)
        self.end_join_circuit(circuit)
        
    def join_circuit_start_id(self, qid, circuit):
        ''' Takes a circuit to extend and attaches it to self. data_qubits specify the where data 
        qubits get mapped. All other qubits
        '''
        n_qubits=len(circuit.data_qubits())
        circuit.replace_qubit_ids(range(qid,qid+n_qubits))
        self.end_join_circuit(circuit)

    def to_ancilla(self, data_qubit_ids):
        next_ancilla = self._next_ancilla()
        for i,q in enumerate(self._get_qubits(data_qubit_ids)):
            q.qubit_type='ancilla'
            q.qubit_id=next_ancilla+i
        self.update_map()
        return range(next_ancilla, next_ancilla+len(data_qubit_ids))
    
    #Main body of functions

    def update_map(self):
        ''' update qubit_gates_map from gates list. If unknown qubit is used in gates, these qubits
        get appended into this circuit.
        '''

        #self.qubit_gates_map = dict.fromkeys(self.qubit_gates_map.keys())
        self.qubit_gates_map = {}
        for g in self.gates:
            
            for q in g.qubits:
                try:
                    self.qubit_gates_map[q].append(g)
                except (KeyError, AttributeError):
                    self.qubit_gates_map[q] = [g]
                # if isinstance(g,Container_Gate):
                #   for circuit in g.circuit_list:
                #       circuit.update_map()

    def clean_qubits(self, qubits={}):
        '''
        Make sure there is no multiple qubit object with the same qubit ids.
        This can prevent errors when we change qubit ids manually.
        This is recursive down to the lowest sub-circuit in each gate.
        '''
        if not qubits:
            qubits = {}
            for q in self.qubits():
                qubits[(q.qubit_id,q.qubit_type)] = q
        for gate in self.gates:
            gate.clean_qubits(qubits)

        self.update_map()
        return qubits

    def create_empty_circuit(self, number_of_qubits=1):
        for i in range(number_of_qubits):
            self.qubit_gates_map[Qubit(i)] = []

    def add_new_qubits(self, number_of_qubits=0):
        max_id = max([int(q.qubit_id) for q in self.qubit_gates_map.keys()])
        qubits = [Qubit(i) for i in range(max_id+1, max_id+number_of_qubits+1)]
        for q in qubits:
            self.qubit_gates_map[q] = []
        return qubits
    
    def remove_unused_qubits(self):
        ''' Removes all the qubits which do not have any gates.
        '''
        for qubit in [q for q in self.qubit_gates_map.keys() if not self.qubit_gates_map[q]]:
            del self.qubit_gates_map[qubit]

    def add_end_gate(self, qubits_or_qubit_ids=[], gate_or_gate_name=''):
        
        data_qubits = self._get_qubits(qubits_or_qubit_ids, 'data', True)
        
        
        if isinstance(gate_or_gate_name, str):
            gate = Physical_Gate(gate_or_gate_name,data_qubits)
        
        elif isinstance(gate_or_gate_name, Gate):
            gate=gate_or_gate_name
            gate.replace_qubit_ids(qubits_or_qubit_ids)
            next_ancilla = self._next_ancilla()
            ancilla_ids= range(next_ancilla,next_ancilla+len(gate.ancilla_qubits()))
            gate.replace_qubit_ids(ancilla_ids,'ancilla')
            
        else:
            print 'Warning! did not pass a string or gate to add_end_gate()'
            raise Exception
        
        self.gates.append(gate)
        self.clean_qubits()

        return gate
    
    # Used in error.py file to insert error syndrome.
    #
    def insert_gate(self, after_gate='', qubits_or_qubit_ids=[], qubit_type='', gate_or_gate_name='',add_qubits_if_not_found=False):
        """
        inserts gates at user specified points in the circuit.

        @type after_gate: Gate object
        @param after_gate: the existing gate after which the new gate should be inserted.
        @type qubits_or_qubit_ids: list of qubit or qubits ids
        @param qubits_or_qubit_ids: qubits which are acted on by the new gate.
        @type qubit_type: string
        @param qubit_type: the type of qubits for the new gate.  This only allows on type to be specified for all qubits.
        @param gate_or_gate_name: the gate or gate name for the new gate.
        @param add_qubits_if_not_found: boolean expression which adds the qubits of the new gate to the circuit if they are not currently there.
        """
        qubits = self._get_qubits(qubits_or_qubit_ids, qubit_type, add_qubits_if_not_found)
        gate = self._get_gate(gate_or_gate_name, qubits)
        #gate.qubits = qubits
        if not after_gate:
            gate_index=-1
        else:
            try:
                gate_index=self.gates.index(after_gate)
            except ValueError:
                em = 'ERROR: the starting gate provided to insert_gate does not exist in this circuit.'
                raise SystemExit(em)
        self.gates.insert(gate_index+1, gate)
        self.update_map()
        return gate

    def end_join_circuit(self, circuit_to_append=None, ancilla_parallel=True):
        if self.level != circuit_to_append.level:
            print 'warning: joining circuits with different levels'

        circuit = circuit_to_append

        if ancilla_parallel:
            anc_qubits = circuit.ancilla_qubits()
            next_ancilla = self._next_ancilla()
            ancilla_ids= range(next_ancilla,next_ancilla+len(anc_qubits))   
            circuit.replace_qubit_ids(ancilla_ids,'ancilla')

        self.gates += circuit_to_append.gates
        start=time.time()
        self.clean_qubits()
        end = time.time()
        #print 'cleaning qubits: ', end-start

        #print 'cleaning qubits: ', time.time()-end
        
    #   Not used yet. If used, make sure it is easy to use.
    #
    # def insert_circuit(self, list_of_before_gates=[], circuit_to_insert=None):
    #   print list_of_before_gates
    #   ind = self._max_ind(list_of_before_gates, self.gates) + 1   
    #   all_before_gates = self.gates[:ind]
    #   self.gates = self.gates[:ind] + circuit_to_insert.gates + self.gates[ind:]
    #   self.update_map()

    def add_error_correction(self, correction_type='EC_CatCorrect', after_prepare=True):
        ''' Inserts specified error correction gates (upper level, without actual sub_circuit)
        after every gate.  No EC is added after measurements and we can choose whether or not
        to add EC after preparations.
        '''
    
        gates = []
        for g in self.gates:
            gates.append(g)
            if g.gate_name[0]!='M':
                if (g.gate_name[:7]!='Prepare'):
                    gates += [Gate(correction_type, [q]) for q in g.qubits]     
                else:
                    if after_prepare:
                        gates += [Gate(correction_type, [q]) for q in g.qubits]

        self.gates = gates
        self.update_map()
        
    def replace_qubit_ids(self, new_qubit_ids, qubit_type='data'):
        if qubit_type=='data':
            qubits = self.data_qubits()
        else:
            qubits = self.data_and_ancilla_qubits()[1]
            
        if len(new_qubit_ids) is not len(qubits):
            print 'too many/less qubit ids are given'
            raise SizeMismatchError
        for i,q in enumerate(qubits):
            q.qubit_id = new_qubit_ids[i]
        #self.clean_qubits()
        return self

    def unpack(self):
        ''' Returns a circuit with all plain gates (no container blocks suck as 'verify')
        '''
        unpacked_circuit = Circuit()
        gates = [g for gate in self.gates for g in gate.unpack()]
        unpacked_circuit.gates = gates
        unpacked_circuit.update_map()
        return unpacked_circuit

    def gates_on_qubit(self, qubit_or_qubit_id=None, qubit_type='data', recursive=False):
        ''' Returns a list of gates the qubit experiences. When recursive=True, unpackes the 
        container gates.
        For example, when there is a 'verify' gate, 
            recursive=False -> returns gates including 'verify' gate, without any gates within it.
            recursive=True  -> returns CNOT, M, and other gates for verification, without 'verify' gate.
        ''' 
        qubit = _get_qubits([qubit_or_qubit_id],'data')[-1] 
        if not recursive:
            return self.qubit_gates_map[qubit]
        unpacked_circuit = self.unpack()
        return unpacked_circuit.qubit_gates_map[qubit]

    def read_gates(self, without_t=True):
        ''' Iterate over gates in without_t or each t (parallel) order. This does not go to check the gates inside
        a gate. For example, this yields a 'verify' gate ignoring the actual gates in its sub-circuit.

         without_t=True ->  iterates over gates in relative order (sufficient ordering).
         without_t=False->  iterates over group of gates that can be performend in parallel order,
                            yielding a tuple of time and list of gates. (t, gates)

        To obtain gates for each qubit (like Circuit.read_gate(serial=True) in the previous circuit.py),
        use gates_on_qubit(qubit, re) instead.
        '''
        if without_t:
            for g in self.gates:
                yield g
            return
        qubits = sorted(self.qubit_gates_map.keys())
        indices = [0] * len(qubits)
        t = 0
        while True:
            
            next_gates = [self.qubit_gates_map[q][i] for q,i in zip(qubits, indices)
                                                if len(self.qubit_gates_map[q]) > i]
                                                
            # indices
            #print [(g.gate_name, [q.qubit_id for q in g.qubits]) for g in next_gates]
            
            if not next_gates:
                return
            ready = []
            for i, gate in enumerate(next_gates):
                # all of this gate.qubits are ready for this gate 
                #print i, gate.gate_name, 'evaluating'
                #print next_gates.count(gate), len(gate.qubits), '<- is the same for one qubit gate'
                #print ''

                if next_gates.count(gate) == len(gate.qubits) and not gate in ready:
                    
                    ready += [gate]
                    for q in gate.qubits:
                        indices[qubits.index(q)] += 1
                    continue
            yield (t, list(set(ready))) # remove duplicates
            t += 1
            
    def qubits(self, ascending=True):
        ''' Returns a list of qubits by ascending or descending order of their ids. 
        '''
        qubits = self.qubit_gates_map.keys()
        qubits = sorted(qubits) if ascending else sorted(qubits)[::-1]
        return qubits
        
    def data_and_ancilla_qubits(self, ascending=True):
        '''Returns a tuple of list of data and ancilla qubits. i.e. ([data qubits], [ancilla qubits])
        '''
        all_qubits = self.qubits(ascending)
        ancilla = []
        data = []
        for q in all_qubits:
            if q.qubit_type == 'ancilla':
                ancilla += [q]
            elif q.qubit_type == 'data':
                data += [q]
            else:
                print 'qubit %r has qubit_type not data nor ancilla'%q
        return data,ancilla
        
    def ancilla_qubits(self):   
        return self.data_and_ancilla_qubits()[1]
        
    def data_qubits(self):
        return self.data_and_ancilla_qubits()[0]
    
    def number_of_gates(self):
        return len(self.gates)

    # Helper funtions from here these are meant to be called only from another Circuit methods.

    def _max_ind(self, new, original):
        return max([original.index(n) for n in new if n in original])

    def _get_qubits(self, qubits_or_qubit_ids=[], qubit_type='data', add_if_not_found=False):
        ''' Returns a list of qubits.
        If qubit_or_qubit_ids is a list of Qubit objects, it returns itself.
        If not, it finds correspinding qubit from qubit_gate_map with the given qubit ids.
        If the given qubit id is not found, it adds a new qubit to this circuit.
        '''
        if not qubits_or_qubit_ids: return []
        if isinstance(qubits_or_qubit_ids[0], Qubit):
            return qubits_or_qubit_ids
        qubits = []
        # gate qubit from qubit ids, if not found, adds new qubits
        for qid in qubits_or_qubit_ids:
            q = [q for q in self.qubit_gates_map.keys() if (q.qubit_id==qid and q.qubit_type==qubit_type)]
            if q:
                qubits += q
                continue
            if add_if_not_found:
                q = Qubit(qid,qubit_type)
                qubits += [q]
                self.qubit_gates_map[q] = []
                continue
            raise KeyError('qubit_id:'+repr(qid)+' is not found in this circuit')
        return qubits


    def _get_gate(self, gate_or_gate_name='', qubits=None):
        ''' Returns a gate from the given gate or gate name. If gate is given, 
        this method simply returns it. Otherwise it creates a new gate and returns it.
        '''
        if not gate_or_gate_name:   return []
        if isinstance(gate_or_gate_name, Gate):
            return gate_or_gate_name
        return Physical_Gate(gate_or_gate_name, qubits)
        
    def _next_ancilla(self):
        largest_id= -1
        for q in self.qubits():
            if q.qubit_type == 'ancilla':
                if q.qubit_id > largest_id:
                    largest_id=q.qubit_id
        return largest_id+1
        
    def _next_data(self):
        largest_id=-1
        for q in self.qubits():
            if q.qubit_type == 'data':
                if q.qubit_id > largest_id:
                    largest_id=q.qubit_id
        return largest_id+1

    # Operator overloading

    def __add__(self, other):
        '''
        Circuit A + Circuit B => combined Circuit C
        C will be a new circuit. Contents in A and B do not change.
        C.gates = A.gates + B.gates.
        '''
        new_circ = Circuit(gates = self.gates[:]+other.gates[:])
        new_circ.update_map()
        return new_circ

    def __sub__(self, other):
        '''
        Cirrcuit A + Circuit B => smalle Circuit C
        C contains all the gates in A except the ones contained in B.
        In other words, B.gates + C.gates == A.gates + B.gates.
        '''
        new_gates = [g for g in self.gates if g not in other.gates]
        new_circ = Circuit(new_gates)
        new_circ.update_map()
        return new_circ

    def __len__(self):
        '''
        Returns a number of gates when len(circuit) is called.
        This is also used for bool(). For example circ become False when there is no gate.
        In this case, you can use:
        if not circ:
            print 'circuit is empty'
        '''
        return len(self.gates)
        
        
class LogicalQubit(object):
    '''Stores information about how qubits interrelated between encoding levels. 
    This is needed whenever we want to measure encoded information.
    '''
    def __init__(self, parent=None, qubit = None, children= None, encoder= None):
        self.parent = parent
        self.children = [] if not children else children
        self.encoding=encoder
        self.qubit = qubit 
        
    def increment_level(self):
        self.qubit.level += 1
        for child in self.children:
            if isinstance(child,LogicalQubit):
                child.increment_level()
        
    def get_qubits_at_level(self, level=1):
        qubits = []
        if self.qubit.level is level:
            qubits += [self.qubit]
        for child in self.children:
            qubits += child.get_qubits_at_level(level)
        return qubits
        
    def get_root(self):
        if self.parent is None:
            return self
        else:
            self.parent.get_root()
    
    def highest_level(self):
        return self.find_root().level
        
    def get_logical(self, qubit):
        if qubit is self.qubit:
            return self
        else:
            if not self.children == []:
                for child in self.children:
                    logical = child.get_logical(qubit)
                    if isinstance(logical,LogicalQubit):
                        return logical
                return None
            else:
                return None
                
    def get_logicals_at_level(self, level, ordered=True):
        logical_qubits = []
        if self.qubit.level == level:
            logical_qubits += [self]
        for child in self.children:
            logical_qubits += child.get_logicals_at_level(level,False)
        if ordered:
            logical_qubits = sorted(logical_qubits)
        return logical_qubits
    
    def add_child(self, qubit_id, qubit_type='data', level=0, encoder=None):
        children = self.make_children([Qubit(qubit_id,qubit_type,level)])
        return children
        
    def make_children(self, qubits, encoder=None):
        children = []
        for q in qubits:
            new_logical = LogicalQubit(self,q)
            children += [new_logical]
            self.children.append(new_logical)
        return children
                
    def _repr(self, highest_level):
        level=self.qubit.level
        return_string = ('\n'+' '*4*(highest_level-level)+str(self.qubit.qubit_id) + 
                        self.qubit.qubit_type[0] + '_'+ str(self.qubit.level))
            #print 'test', node.children
        if self.children !=None:
            for child in self.children:
                return_string += child._repr(highest_level)
        return return_string
        
        
    def __repr__(self):
        return ('< Logical Qubit( qubit= ' + repr(self.qubit.qubit_id) + 
                    self.qubit.qubit_type[0] + '_'+ repr(self.qubit.level) + 
                    ', Children= '+ repr(self.children)+'>')
                    
    def __cmp__(self, other):
        return self.qubit.__cmp__(other.qubit)

                    
class LogicalQubitList(object):
    def __init__(self, logical_qubits=None):
        self.logical_qubits = [] if not logical_qubits else logical_qubits
        
    def increment_level(self):
        for logical_qubit in self.logical_qubits:
            logical_qubit.increment_level()
    
    def get_qubits_at_level(self, level=1):
        qubits = []
        for logical_qubit in self.logical_qubits:
            root_node = logical_qubit.get_root()
            qubits += root_node.get_qubits_at_level(level)
        return qubits
            
    def highest_level(self):
        highest_level = 0
        for logical_qubit in self.logical_qubits:
            root_node = logical_qubit.get_root()
            if root_node.qubit.level > highest_level:
                highest_level = root_node.qubit.level
        return highest_level
    
    def add_root_qubit(self, qubit):
        root = LogicalQubit(qubit=qubit)
        self.logical_qubits += [root]
        return root
        
    def add_root_qubits(self, qubits):
        new_logical_qubits=[]
        for q in qubits:
            new_logical_qubits += [self.add_root_qubit(q)]
        return new_logical_qubits

    def get_logical(self,qubit):
        for logical_qubit in self.logical_qubits:
            root_node = logical_qubit.get_root()
            node = root_node.get_logical(qubit)
            if isinstance(node,LogicalQubit):
                return node
        return None
        
    def get_logicals_at_level(self, level, ordered=True):
        logical_qubits = []
        for logical_qubit in self.logical_qubits:
            root_node = logical_qubit.get_root()
            logical_qubits += root_node.get_logicals_at_level(level,False)
        #print logical_qubits
        #for logical_qubit in logical_qubits:
        #   print logical_qubit, id(logical_qubit)
        #print '\n'
        #print id(logical_qubits)
        if ordered:
            logical_qubits.sort()
        #print logical_qubits
        #for logical_qubit in logical_qubits:
        #   print logical_qubit, id(logical_qubit)
        #print '\n'
        #print id(logical_qubits)
        return logical_qubits
            
                
    def add_children_to_logical(self, qubit, children):
        logical = self.get_logical(qubit)
        logical.children.append(children)
        return logical
        
    def next_id_at_level(self,level, qubit_type='data'):
        next_id=-1
        for logical_qubit in self.logical_qubits:
            root_node = logical_qubit.get_root()
            qubits= root_node.get_qubits_at_level(level)
            for q in qubits:
                if q.qubit_type is qubit_type:
                    if q.qubit_id>next_id:
                        next_id=q.qubit_id
        return next_id +1
    
    def check_form(self):
        pass
        
    def __getitem__(self,i):
        return self.logical_qubits[i]
        
    def __repr__(self):
        return_string = 'Logical Qubit List:'
        highest_level= self.highest_level()
        for logical in self.logical_qubits:
            node = logical.get_root()
            return_string += node._repr(highest_level)
            return_string += '\n'+'----'*(highest_level+1)
        return return_string
        
        
class Encoded_Circuit(Circuit):
    def __init__(self, gates=None, qubits=None, level=0, logical_qubits=None):
        Circuit.__init__(self, gates, qubits, level)
        self.logical_qubits = [] if not logical_qubits else logical_qubits
        
    @classmethod
    def create_encoded_circuit(cls, circuit, logical_qubits = None):
        encoded_circuit = Encoded_Circuit(circuit.gates,circuit.qubits(),circuit.level,logical_qubits)
        encoded_circuit.update_map()
        return encoded_circuit
    
    def return_circuit_instance(self):
        print self.qubits()
        return Circuit( gates=self.gates, qubits=self.qubits(), level=self.level)
        
