"""
chper.py

A wrapper module for using CHP program by Scott Aaronson.

Yu Tomita <yu.t@gatech.edu>
Oct 19 2011

With new functionalities added by Mauricio Gutierrez <mga7@gatech.edu>
Sept 2012
The new functionalities are basically the capability to start a chp run with 
a list of stabilizers (stabs) and destabilizers (destabs) instead of a ket 
state and likewise return a list of stabs and destabs at the end of the run.
"""

import sys
import os
import time as t
import random as rd
import circuit
from visualizer import browser_vis as brow
from subprocess import Popen, PIPE

class Chper(object):
	def __init__(self, circ, num_d_qub, num_a_qub, stabs, destabs,
	 	     anc_stabs, anc_destabs, states='None', input_output_files=True,
		     final_stab_init_i=0):
		"""
		Creates a chp runner from given circuit. This saves a list of chp gates 
		translated from the given circuit. If the circuit does not change and 
		when we need to run chp multiple times, it is better to keep this object 
		and call chper.run() again and again.
		This new method requires the user to input some extra parameters:
		   - num_d_qub: the number of data qubits in the circuit
		   - num_a_qub: the number of ancilla qubits in the circuit
		   - states:    'None' or a string.
		   - stabs:     a list of the stabilizers of the initial state. 
				For example, if our initial state is |000>, 
				stabs = ['+ZII', '+IZI', '+IIZ']
		   - destabs:   a list of the destabilizers.
		   - anc_stabs and anc_destabs:  same thing
                   Notice that stabs and destabs can either cover only the data qubits
		   or the whole circuit.  		

		     If the circuit contains prepare gates at the start, then the user 
		     doesn't need to specify stabilizers or destabilizers.
		   - input_output_files:  True if we want CHP to get the input from (txt) 
					       files and write the output on a file as well.
					  False if we don't want to use any files for I/O.
		
		Right now, we're assuming that the input states will always be None.
		Change this later if needed.

		final_stab_init_i:  the initial index of the final stabilizers.
				    This option was included for the Knill correction,
				    where the data gets teleported to the last qubits
				    and we need to get rid of the initial ones.

		NOTICE:  I DECIDED TO NOT HAVE DEFAULT VALUES FOR ANY LIST VARIABLE GIVEN 
			 THE BUG I HAD WITH THE STABILIZERS.  THIS HAS TO DO WITH LISTS
			 BEING MUTABLE.  WE COULD CHANGE EVERY LIST TO A TUPLE, BUT I
			 THINK IT'S EASIER TO FORCE THE USER TO DEFINE EVERY VARIABLE
			 WHEN INITIALIZING AN OBJECT.  MGA MARCH 2015
		"""
		self.replacements = {'C':[('C', [0,1])],
			'CX':[('C', [0,1])],
			'CNOT':[('C', [0,1])],
			'CZ':[('H', [1]), ('C', [0,1]), ('H', [1])],
			'CY':[('P', [1]), ('C', [0,1]), ('P',[1]),('P',[1]),('P',[1])],
            		'H':[('H', [0])],
			'P':[('P', [0])],
			'X':[('H', [0]),('P', [0]),('P', [0]),('H', [0])],
			'Y':[('P', [0]),('H', [0]),('P', [0]),('P', [0]),('H', [0]),('P', [0]),('P',[0]),('P',[0])],
			'Z':[('P', [0]),('P', [0])],
			'Measure':[('M', [0])],
			'MeasureZ':[('M', [0])],
            		'MeasureZDestroy':[('M', [0])],
			'MeasureX':[('H', [0]),('M', [0])],
            		'MeasureXDestroy':[('H', [0]),('M', [0])],
			'M':[('M', [0])],
			'I':[],
            'ImZ': [],
            'ImX': [],
            'I_stark': [],
            'II_heat': [],
            'Ism': [],
            'Icool': [],
            'I_idle': [],
            'I_gate': [],
            'I_cross': [],
			'IMS5':[],
            'RZ +': [('P', [0])],
            'RZ -': [('P', [0]),('P', [0]),('P', [0])],
            'RX +': [('H', [0]), ('P', [0]), ('H', [0])], 
            'RX -': [('H', [0]), ('P', [0]), ('P', [0]), ('P', [0]), ('H', [0])], 
            'RY +': [('P', [0]), ('H', [0]), ('P', [0]), ('P', [0]), ('P', [0]), ('H', [0]),
                     ('P',[0]), ('P',[0]), ('P',[0])],
            'RY -': [('P', [0]), ('H', [0]), ('P', [0]), ('H', [0]), ('P',[0]), ('P',[0]), ('P',[0])], 
            'MS': [('H',[0]), ('C',[0,1]), ('P',[0]), ('H',[0]),
                   ('H',[1]), ('P',[1]), ('H',[1])]
            }
		self.preparations ={"PrepareXPlus":'x', "PrepareXMinus":'X', "PrepareYPlus":'y',
				            "PrepareYMinus":'Y', "PrepareZPlus":'z', "PrepareZMinus":'Z'}
		self.stab_dic = {'z':'+Z', 'Z':'-Z', 'x':'+X', 'X':'-X'}
		self.destab_dic = {'z':'+X', 'Z':'+X', 'x':'+Z', 'X':'+Z'}
		self.circ = circ
		#brow.from_circuit(circ, True)
		#t.sleep(15)
		self.num_d_qub = num_d_qub
		self.num_a_qub = num_a_qub
		#print num_d_qub, num_a_qub
		self.num_qub = num_d_qub + num_a_qub
		self.final_stab_init_i = final_stab_init_i 
		self.qubit_num_map = {}	# key = Qubit, value = chp number
		self.measurement_gates = {} # key is chp number
		self.results = {}
		#print stabs
        	#print stabs, destabs, anc_stabs, anc_destabs
		#t.sleep(5)
		self.stabs, self.destabs = stabs[:], destabs[:]
		self.anc_stabs, self.anc_destabs = anc_stabs[:], anc_destabs[:]
		if len(self.stabs) == 0:
			if len(self.anc_stabs) == 0:          self.case = 5
			else:                                 self.case = 3
		else:
			if len(self.stabs) == self.num_qub:   self.case = 1
			else:
				if len(self.anc_stabs) == 0:  self.case = 2
				else:                         self.case = 4
		self.input_states = self.set_input_states(states)[:]
		self.final_stabs = []
		self.final_destabs = []

		#print self.case
		#print self.input_states
		
		#print 'qubit_num_map', self.qubit_num_map

		self.input_output_files = input_output_files
		if self.input_output_files:
			self.input_file, self.n_gates = self.make_input_file()
			self.arrange_stabilizers()
			self.input_stab_file = self.make_input_stab_file()
		else:
			self.input_circuit_string, self.n_gates = self.make_input_circuit_string()
			self.arrange_stabilizers()
			self.input_stab_string = self.make_input_stab_string()
	
			#print 'stabs =', self.stabs
			#print 'destabs =', self.destabs
			#print 'case =', self.case
			#print 'input circ=', self.input_circuit_string
			#print 'input stab=', self.input_stab_string
		

		#print 'Initialization was successful :)'
		#t.sleep(20)



	def set_input_states(self, states=None):
		"""
		set input states with given states.
		States need to be a list of characters (strings) ['z','Z',...] 
		or a string 'zZxxy...'
		"""
		if ((type(states) is str) and (states != 'None')):
			return [c for c in states]
		else:
			if (self.case == 1) or (self.case == 4):
				in_states = []
			elif self.case == 2:
				in_states = ['z' for i in range(self.num_a_qub)]
			elif self.case == 3:
				in_states = ['z' for i in range(self.num_d_qub)]
			elif self.case == 5:
				in_states = ['z' for i in range(self.num_qub)]

		return in_states



	def circ_to_chp_gates(self):
		"""
		data holder of gates for chp input files. contains only gate string 
		and a list of qubit numbers.
		"""
		#for i,q in enumerate(self.circ.qubits()):
		#	self.qubit_num_map[q] = i
		
		max_data_q = max([q.qubit_id for q in self.circ.data_qubits()])
		min_data_q = min([q.qubit_id for q in self.circ.data_qubits()])
		if max_data_q+1 > self.num_d_qub:
			ref_data_q = min_data_q
		else:
			ref_data_q = 0 
		
		if self.num_a_qub > 0:	
			first_anc_id = self.circ.ancilla_qubits()[0].qubit_id
		for i,q in enumerate(self.circ.qubits()):
			if q.qubit_type == 'data':
				self.qubit_num_map[q] = q.qubit_id - ref_data_q
			elif q.qubit_type == 'ancilla':
				self.qubit_num_map[q] = q.qubit_id-first_anc_id \
							+self.num_d_qub		
		#print first_anc_id
		#print 'qubit_num_map', self.qubit_num_map
		gates = []

                if isinstance(self.circ.gates[0], circuit.Physical_Gate):
                        """Assumption : if first gate is physical, then the entire 
			circuit can be treated as physical."""
                        for gate in self.circ.gates:
				#print gate.gate_name, gate.qubits[0].qubit_id
                                if 'Measure' in gate.gate_name:
                                        q = gate.qubits[0]
                                        self.measurement_gates[self.qubit_num_map[q]] = gate
                                gates += self.translate_gate(gate)
				
                        return gates
                else:
                        """Trying to unpack the logical circuit to physical gates."""
                        unpacked_circuit = self.circ.unpack()
                        read_gates = unpacked_circuit.read_gates()
                        for i in range(0, unpacked_circuit.number_of_gates()):
                                next_gate = read_gates.next()
                                while type(next_gate) != list:
                                        next_gate = next_gate.unpack()
                                        gate = next_gate[0]
                                        if "orrect" in gate.gate_name or "Encoded" in gate.gate_name:
                                            """Skip the gate."""
                                        else:
                                                if 'Me' in gate.gate_name:
                                                        q = gate.qubits[0]
                                                        self.measurement_gates[self.qubit_num_map[q]] = gate
                                                gates += self.translate_gate(gate)
                        return gates

	def translate_gate(self, gate):
		"""
		Given an Gate, returns translated list of Gates with only H,C,P,M gates.
		@type chp_gate: L{Gate}
		@param chp_gate: an L{Gate} object with any gate
		@rtype: list of L{Gate}
		@return: a list of L{Gate} oinly containig C, H, P, and M gates.
		"""
		if gate.gate_name in self.preparations:
			#print gate.gate_name, gate.qubits[0].qubit_id
			ind = self.qubit_num_map[gate.qubits[0]]
			#print gate.qubits[0]
			#oprint ind
			#if ind < len(self.input_states):
			#	self.input_states[ind] = self.preparations[gate.gate_name]
			#else:
			letter = self.preparations[gate.gate_name]			
			
			if letter != 'z':
				if (self.case == 3) or (self.case == 5):
					self.input_states[ind] = letter
				elif self.case == 2:
					self.input_states[ind-self.num_d_qub] = letter

			#print gate.gate_name, gate.qubits[0].qubit_id
			#print "Input states: ",self.input_states			
			return []	

		try:			
			gates = self.replacements[gate.gate_name]
			
		except KeyError:
			em = ('Gate %s is not C,H,P or M, and we do not know how'
				'to translate it into those gates. Please check "replacements"'
				'dict in monte/chper.py'%gate)
			raise SystemExit(em)
		new_gates = [circuit.Gate(g[0], [gate.qubits[qid] for qid in g[1]]) for g in gates]
		return new_gates



	def make_input_file(self, filename='/tmp/messtempinput'):
		"""
		Saves a chp input file from given general circuit. Filename can be specified 
		or will be a temp file in /tmp/ folder. If a file already exist it overwrites.
		@type circuit: L{Circuit}
		@param circuit: a L{Circuit} object to make CHP input file from.
		@type filename: string
		@param filename: filename of an CHP input file to be saved. This file is 
		overwritten if already exists.
		@rtype: string
		@return: filename of CHP input file created.
		"""
		gates = self.circ_to_chp_gates()[:]
		infile = open(filename, 'w')
		infile.write('#\n')
		for gate in gates:
			s = '%s'%(gate.gate_name)
			for num in [self.qubit_num_map[q] for q in gate.qubits]:
				s += ' %d'%(num)
			s += '\n'
			infile.write(s)
		infile.close()

		return filename, len(gates)



	def make_input_circuit_string(self):
		'''
		One of the new functions added to deal with non-file I/O CHP.  MGA Feb 2015
		Does 2 things:
			(1) Decomposes a general circuit into gates that CHP can handle,
			    just like the old make_input_file
			(2) Creates a string appropriate for CHP to read.  The string has to
			    to have the following format.  For example, if we want to have a
			    circuit that prepares a 2-qubit cat state and measures both qubits
			    H 0
			    C 0 1
			    M 0
			    M 1	
			    the actual string would be 'H_0nC_0_1nM_0nM_1n'.  
			    Notice that each gate is separated by an 'n' and we have an 
			    underscore between the gate name and the qubit numbers.
		'''
		gates = self.circ_to_chp_gates()[:]
		s = ''
		for gate in gates:
			s += gate.gate_name
			for num in [self.qubit_num_map[q] for q in gate.qubits]:
				s += '_%d'%(num)
			s += 'n'

		return s, len(gates)
		


	def make_input_stab_file(self, filename='/tmp/messtempinputstatet'):
		"""
		Makes the file with the destabilizers and stabilizers
		"""
		destabs_plus_stabs = self.destabs[:] + self.stabs[:]
		infile = open(filename, 'w')
		for operator in destabs_plus_stabs:
			infile.write(operator)
			infile.write('\n')
		infile.close()
		
		return filename
	


	def make_input_stab_string(self):
		'''
		Another function added to deal with non-file I/O CHP.  MGA Feb 2015
		Prepares the string with the destabilizers and stabilizers.
		For the previous example of the 2-qubit cat state, we want to start with
		|00>.  Therefore, the string of destabs and stabs would be:
		'+XI+IX+ZI+IZ'.  Notice that the first two are the destabs and the last 
		two are the stabs.
		'''
		destabs_plus_stabs = self.destabs[:] + self.stabs[:]
		s = ''
		for operator in destabs_plus_stabs:  
			s += operator	
		
		return s
	


	def translate_states_to_stabilizers(self):
		"""
		Returns the destabilizers and stabilizers for a given input state.
		"""
		destabils, stabils = [], []
		for i,state in enumerate(self.input_states):
			destabil = ['I' for j in range(len(self.input_states)+1)]
			stabil = ['I' for j in range(len(self.input_states)+1)]	
			destabil[0] = self.destab_dic[state][0]
			stabil[0] = self.stab_dic[state][0]
			destabil[i+1] = self.destab_dic[state][1]
			stabil[i+1] = self.stab_dic[state][1]
			destabils += [''.join(destabil)]
			stabils += [''.join(stabil)]
	
		return destabils, stabils		
			


	def combine_stabilizers(self, data_stabs, data_destabs, anc_stabs, anc_destabs):
		'''
		combines data and ancillary stabilizers and destabilizers by
		adding Is at the end of the data and in front of the ancilla
		and them joining the two lists
		'''
		n_data, n_anc = len(data_stabs), len(anc_stabs)
		for i in range(n_data):
			data_stabs[i] += ''.join(['I' for j in range(n_anc)])
			data_destabs[i] += ''.join(['I' for j in range(n_anc)])
		extra_Is = ''.join(['I' for i in range(n_data)])
		for i in range(n_anc):
			anc_stabs[i] = anc_stabs[i][0]+extra_Is+anc_stabs[i][1:]	
			anc_destabs[i] = anc_destabs[i][0]+extra_Is+anc_destabs[i][1:]	

		comb_stabs = data_stabs[:] + anc_stabs[:]
		comb_destabs = data_destabs[:] + anc_destabs[:]

		return comb_stabs, comb_destabs 



	def arrange_stabilizers(self):
		"""
		Adds the stabilizers of the ancilla qubits to the stabilizers of the 
		data qubits.
		"""
		if self.case == 1:  return None

		trans_destabs, trans_stabs = self.translate_states_to_stabilizers()

		if self.case == 5:
			self.stabs += trans_stabs[:]
			self.destabs += trans_destabs[:]
			return None
		
		if self.case == 2:
			self.anc_stabs = trans_stabs[:]
			self.anc_destabs = trans_destabs[:]
		elif self.case == 3:
			self.stabs = trans_stabs[:]
			self.destabs = trans_destabs[:]			
		
		self.stabs, self.destabs = self.combine_stabilizers(
						self.stabs, self.destabs, 
						self.anc_stabs, self.anc_destabs)
	
		return None  



	def parse_output(self, filename_or_string, output_file=True):
		"""
		Returns a dictionary of qubit number and measurements from the given 
		chp output filename.
		@type filename_or_string: string
		@param filename_or_string: filename of an existing CHP output file
					   or string with the output itself
		@type output_file:   Boolean
		@param output_file:  True if the output is a file
				     False if it is a string
		@rtype:     tuple 
		@return:    (dictionary of key = qubit number (int) and measurement-random 
		            tuple (int, boolean), 
			    list of the final data stabilizers,
			    list of the final data destabilizers
			    )
		"""
		if output_file:
			f = open(filename_or_string, 'r')
			list_lines = f.readlines()
			f.close()
		else:
			whole_string = filename_or_string
			list_lines = filename_or_string.split('\n')
		
		destabs_stabs = []
		STAB = 0
		meas = {}
		for pre_line in list_lines:
			line = pre_line.split()
			#print 'line =', line
			if 'Outcome' in line:
				q = int(line[4].split(':')[0])
				#m = 1 if not int(line[5]) else -1
				m = int(line[5])
				r = '(random)' in line
				#gate = self.measurement_gates[q]
				#meas[gate] = (m, r)
				meas[q] = (m, r)
			
			if STAB:
				destabs_stabs += line
			
			if 'Final' in line:
				STAB = 1

		destabs = destabs_stabs[:((len(destabs_stabs)-1)/2)]
		stabs = destabs_stabs[((len(destabs_stabs)+1)/2):]
	
		self.final_destabs = destabs[:]
		self.final_stabs = stabs[:]

		initial_i = self.final_stab_init_i
		final_i = initial_i + self.num_d_qub		

		data_destabs = destabs[initial_i : final_i]
		data_stabs = stabs[initial_i : final_i]


		for i in range(len(data_stabs)):
			data_destabs[i] = data_destabs[i][0] + \
					  data_destabs[i][initial_i+1 : final_i+1]
			data_stabs[i] = data_stabs[i][0] + \
					data_stabs[i][initial_i+1 : final_i+1]
		

		#print 'parsing...'
		#print 'meas:', meas
		#print 'data stabs:', data_stabs
		#print 'data destabs:', data_destabs
		#t.sleep(5)

		self.set_values_back_to_default()

		#print 'len of stabs =', len(stabs)
		#print 'number of data qubits =', self.num_d_qub
		#print 'stabs at the end =', stabs
		#print 'self.stabs =', self.stabs
		#print 'self.destabs =', self.destabs
		#print 'self.anc_stabs =', self.anc_stabs
		#print 'self.anc_destabs =', self.anc_destabs

		return (meas, data_stabs, data_destabs)


	def set_values_back_to_default(self):
		'''
		set all the attributes' values back to default to avoid bugs related to
		mutability of lists
		'''
		self.states = 'None'
		self.stabs, self.destabs = [], []
		self.anc_stabs, self.anc_destabs = [], []
		
		return None



	def run(self, chp_location='./chp_extended', output_filename='/tmp/messtempoutput',
		verbose=False):
		"""
		Runs chp with the given input file name and the directory/filename of chp program.
		states is a string of initial qubit states. 

			- z for |0> 		(this is default)
			- Z for |1>
			- x for |0>+|1>
			- X for |0>-|1>
			- y for |0>+i|1>
			- Y for |0>-i|1>

		for example, states = 'zzzzz' = |00000> (this time it's same as states = '')

		@type input_filename: string
		@param input_filename: filename of existing CHP input file
		@type output_filename: string
		@param output_filename: filename of CHP output file to be saved
		@type chp_location: string
		@param chp_location: directory and program name of CHP program (default='./chp')
		@rtype: L{dict}
		@return: dictionary of key = qubit number (int) and measurement-random tuple (int, boolean).
		"""

		# New part: this is where we generate the random number
		# that will be used as seed in CHP.  MGA: 03/02/2017.
		r_num = rd.random()
		r_num_int = int(r_num*10**8)
		
		#print 'rand int =', r_num_int
		#destabs = ''.join(self.destabs)
		#stabs = ''.join(self.stabs)
		#instates = destabs + stabs
		num_ancilla = len(self.circ.ancilla_qubits())	
		num_total = self.num_d_qub + num_ancilla
		#command = ('%s -q %s %s %s'%(chp_location, self.input_file, 
		#				instates, num_ancilla)
	 	
		#print stabs
		#print num_ancilla
		#print num_total

		#print 'chp location =', chp_location
		#print 'input file =', self.input_file	
		#print 'input stab file =', self.input_stab_file
		#print 'output filename =', output_filename

		#print 'chp location =', chp_location
		#print self.input_circuit_string
		#print self.input_stab_string
		#print 'num ancilla =', str(num_ancilla)
		#print 'num total =', str(num_total)
		#print 'num gates =', str(self.n_gates)	
			
	
		if self.input_output_files:

			command = ('%s -q %s %s %s %s %s > %s'%(chp_location, self.input_file, 
							self.input_stab_file, str(r_num_int),
							str(num_ancilla), str(num_total), 
							str(self.n_gates), output_filename))

			b = os.system(command)
	
			if verbose:
				print 'input state is this', self.input_states
				print command 
			if b:
				em = 'Error running chp\n'
				if b==32512:
					em+= '  chp program not found at %s\n'%chp_location
				else:
					em+= '  input file %s is contaminated\n'%self.input_file
					a = open(self.input_file)
					print self.input_file
					for i in a:
						print i
				raise SystemExit(em)
			
			return self.parse_output(output_filename, output_file=True)
		

		else:
			#print 'circ =', self.input_circuit_string
			#print 'stab =', self.input_stab_string
			chp_inputs = [chp_location, '-lq', self.input_circuit_string,
				      self.input_stab_string, str(r_num_int),
				      str(num_ancilla), str(num_total), str(self.n_gates)]
			chp_proc = Popen(chp_inputs, stdin=PIPE, stdout=PIPE)
			output_string, error = chp_proc.communicate()
			#print output_string
			#chp_proc.kill

			if error != None:
				em = 'Error running chp'
				raise SystemExit(em)			

			return self.parse_output(output_string, output_file=False)




	
