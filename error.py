"""
error.py

Handles error of monte carlo simulation, including reading error file and adding errors with given probabilities.

Yu Tomita <yu.t@gatech.edu>
Oct. 19 2011

plus some additions by Mauricio Gutierrez <mga7@gatech.edu>
Dec. 3, 2012
"""
import sys
import json
import random
from math import exp

class ErrorRates(object):
	"""
	Data holder of error rates. 
	"""
	def __init__(self,rates,onequbit,twoqubit,errorkind):
		self.rates=rates
		self.one_qubit = onequbit
		self.two_qubit = twoqubit
		self.error_kind = errorkind

	def __repr__(self):
		return repr([self.rates,self.one_qubit,self.two_qubit,self.error_kind])



class ErrorRatesAlternative(object):
	"""
	Data holder of error rates with the alternative format. 
	"""
	def __init__(self, error_dict, error_kind=1):
		if 'error_kind' in error_dict.keys():
			self.error_kind = error_dict.pop('error_kind', None)
		else:
			self.error_kind = error_kind

		self.dic = error_dict



def read_error(filename):
	"""
	Reads error information from file and returns it as an L{ErrorRates} object.
	@type filename: string
	@param filename: filename of a json format file with error information. Example file is below:
		{
			"error_rates":{
			   "X": 0.5,
			   "C": 0.02, 
			   "I": 0.02,
			   "H": 0.2,
			   "M": 0.5,
			   "C_elu": 0.8
			},
			"one_qubit_ratio":{
			   "X": 3,
			   "Y": 1,
			   "Z": 1.5
			},
			"two_qubit_ratio":{
			   "XX": 1,
			   "XY": 0.1,
			   "ZZ": 2.5
			   "C" : 1
			}
			"error_kind": 3
		}	

	Here, error_rates contains the rate of having an error after each type of gate. One_ and two_qubit_ratio contains one 
	qubit error gates and the ratio of their happening.  In the above example having a "X" Pauli error is three times more 
	likely than having a "Y" error.
	Finally, error_kind refers to 1 of 3 possibilities:
		(1) Errors can occur after any gate.
		(2) Errors only occur after gates that are NOT involved in error correction.  
		In other words, the error detection and correction steps are perfect.
		(3) Errors can occur during error correction steps but only on data qubits.
	
	@rtype: L{ErrorRates}
	@return: ErrorRates object with the read error information
	"""
	try:
		es = json.load(open(filename))
		rates = es['error_rates']	
		# [('X',1),('Y',1.5)...] str() is needed to convert from unicode in json
		one_qubit = [(str(x),y) for x,y in es['one_qubit_ratio'].items()]
		two_qubit = [(str(x),y) for x,y in es['two_qubit_ratio'].items()]
		error_kind = es['error_kind']
	except IOError:
		em = 'ERROR: error info file name %s does not exist.'%filename
		raise SystemExit(em) 
	return ErrorRates(rates,one_qubit,two_qubit,error_kind)	
	


def pick_error(error_ratio):
	"""
	Weighted random number selector.

	@type error_ratio: list of tuples
	@param error_ratio: list of tuples of two values (key, value) where key is selected with the value. Values can be any non-negative number,
	and they do not have to be normalized.
	@return: randomly selected key with its weight
	"""
	key,value = zip(*error_ratio)
	r = random.random()*sum(value)
	for i,v in enumerate(value):
		r-=v
		if r<0:
			return key[i]


def add_error(circ, error_info_or_filename):
	if type(error_info_or_filename) == str:
		error_info = read_error(error_info_or_filename)
	else:
		error_info = error_info_or_filename

	if hasattr(error_info, 'rates'):
		#print '1'
		#sys.exit(0)
		add_error_traditional(circ, error_info)
	else:
		#print '2'
		#sys.exit(0)
		add_error_alternative(circ, error_info)

	return



def add_error_traditional(circ, error_info, sampling='old', gate_indices=[]):
	'''
	'''
	if sampling == 'old':
	    for g in circ.gates[::-1]:
		    # The "for" is done in reverse order to prevent an
		    # infinite loop as error gates are added to the gate list. Mauricio 
		    try:
			    rate = error_info.rates[g.gate_name]
		    except KeyError: # no error info found after this gate
			    continue 
		    r = random.random()
		    if r<rate:
			    kind = error_info.error_kind
			    if len(g.qubits) == 1:
				    if ((kind==2 or kind==3) and g.qubits[0].qubit_type=='ancilla'):
					    continue
				    error_ratio = error_info.one_qubit
				    new_g = circ.insert_gate(g, g.qubits,'', pick_error(error_ratio), False)
				    new_g.is_error = True
			    else:
			        # This "else" was added later to insert two 1-qubit gates instead of one 
                    		# 2-qubit gate, which the chp wrapper did not know how to interpret.  
			        # We might have to change this again if the two-qubit error cannot be 
                    		# expressed as two 1-qubit gates.   Mauricio
				    error_ratio = error_info.two_qubit
				    q1, q2 = g.qubits[0], g.qubits[1]
				    if (kind==2 and (q1.qubit_type=='ancilla' or q2.qubit_type=='ancilla')):
					    continue
				    elif kind==3:
					    error_gate = pick_error(error_ratio)
					    if q2.qubit_type=='data':
						    new_g = circ.insert_gate(g, [q2], '', error_gate[1], False)
						    new_g.is_error = True	
					    if q1.qubit_type=='data':
						    new_g = circ.insert_gate(g, [q1], '', error_gate[0], False)
						    new_g.is_error = True
				    else:
					    error_gate = pick_error(error_ratio)
					    for q in g.qubits[::-1]:
						    new_g = circ.insert_gate(g, [q], '', 
                                                            error_gate[g.qubits.index(q)], False)
						    new_g.is_error = True

	elif sampling == 'Muyalon':
	    gate_indices.sort()  # this step is critical
	    for i in gate_indices[::-1]:
			g = circ.gates[i]
			if len(g.qubits) == 1:
				error_ratio = error_info.one_qubit
				new_g = circ.insert_gate(g, g.qubits, '', pick_error(error_ratio), False)
				new_g.is_error = True
			else:
				error_ratio = error_info.two_qubit
				error_gate = pick_error(error_ratio)
				new_g = circ.insert_gate(g, [g.qubits[1]], '', error_gate[1], False)
				new_g.is_error = True
				new_g = circ.insert_gate(g, [g.qubits[0]], '', error_gate[0], False)
				new_g.is_error = True

	return


def add_error_alternative(circ, error_info, sampling='old', gate_indices=[]):
	'''
	'''
	if sampling == 'old':
		for g in circ.gates[::-1]:
			# The "for" is done in reverse order to prevent an
			# infinite loop as error gates are added to the gate list. Mauricio 
			try:
				rate = error_info.dic[g.gate_name]['error_rate']
				# if the gate has a duration, it means the error rate should
				# be multiplied by the duration.
				# This should be improved later on, because it assumes that
				# the error rate scales linearly with the gate's duration.
				# This is true for incoherent errors like heating, but for
				# coherent errors like the Stark shifts, it should be quadratic.
				#if type(g.duration) == type(0.1):
				#	rate = rate*g.duration
			except KeyError: # no error info found after this gate
				continue 
			r = random.random()
			if r<rate:
				error_ratio_dic = error_info.dic[g.gate_name]['error_ratio']
				error_ratio = [(str(x),y) for x,y in error_ratio_dic.items()]
				kind = error_info.error_kind
				if len(g.qubits) == 1:
					if ((kind==2 or kind==3) and g.qubits[0].qubit_type=='ancilla'):
						continue
					new_g = circ.insert_gate(g, g.qubits,'', pick_error(error_ratio), False)
					new_g.is_error = True
				else:
			        	# This "else" was added later to insert two 1-qubit gates instead of one 
                    			# 2-qubit gate, which the chp wrapper did not know how to interprete.  
			        	# We might have to change this again if the two-qubit error cannot be 
                    			# expressed as two 1-qubit gates.   Mauricio
				    	q1, q2 = g.qubits[0], g.qubits[1]
				    	if (kind==2 and (q1.qubit_type=='ancilla' or q2.qubit_type=='ancilla')):
					    	continue
				    	elif kind==3:
					    	error_gate = pick_error(error_ratio)
					    	if q2.qubit_type=='data':
							new_g = circ.insert_gate(g, [q2], '', error_gate[1], False)
							new_g.is_error = True	
					    	if q1.qubit_type=='data':
						    new_g = circ.insert_gate(g, [q1], '', error_gate[0], False)
						    new_g.is_error = True
			    		else:
						error_gate = pick_error(error_ratio)
						for q in g.qubits[::-1]:
							new_g = circ.insert_gate(g, [q], '', 
                                        		error_gate[g.qubits.index(q)], False)
						    	new_g.is_error = True
                return

	elif sampling == 'Muyalon':
	    
            errors_added = []
            gate_indices.sort()  # this step is critical
	    for i in gate_indices[::-1]:
			g = circ.gates[i]
			error_ratio_dic = error_info.dic[g.gate_name]['error_ratio']
			error_ratio = [(str(x),y) for x,y in error_ratio_dic.items()]
			if len(g.qubits) == 1:
				new_g = circ.insert_gate(g, g.qubits, '', pick_error(error_ratio), False)
				new_g.is_error = True
                                  
                                # MGA 12/23/19.  New addition for the NN decoder
                                errors_added += [[i, new_g.gate_name]]
			else:
				error_gate = pick_error(error_ratio)
                                # MGA 12/23/19.  New addition for the NN decoder
                                errors_added += [[i, error_gate]]
				for q_index in range(len(g.qubits))[::-1]:
					new_g = circ.insert_gate(g, [g.qubits[q_index]], '', error_gate[q_index], False)
					new_g.is_error = True
	return errors_added


def add_decoherence(circ, error_info_or_filename):
        """Add errors on scheduled circuits based on the time that qubits sit idle."""
        decoherence_rate=1e-3
        if type(error_info_or_filename) == str:
		error_info = read_error(error_info_or_filename)
	else:
		error_info = error_info_or_filename
		
        for qubit in circ.qubit_gates_map:
            for i, gate in enumerate(circ.qubit_gates_map[qubit]):
                if i>0:
                        #print gate.start_time, prev_gate.end_time
                        idle_time = gate.start_time-prev_gate.end_time
                        if type(idle_time)is int or type(idle_time)is float:
                                decoherence_prob = 1-exp(-decoherence_rate*idle_time)
                                r = random.random()
                                if r < decoherence_prob:
                                        error_ratio = error_info.one_qubit
                                        circ.insert_gate(gate, [qubit], '', pick_error(error_ratio), False)
                prev_gate = gate
	return
