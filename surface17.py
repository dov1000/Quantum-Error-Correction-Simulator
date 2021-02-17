'''
surface17.py
'''

import sys
import os
import json
from circuit import *
import correction



class Code:
    '''
    Defines constants and methods useful for dealing with the Surface-17 code.
    '''

    code_name = 'Surface17'

    stabilizers = [
                   [('X',1), ('X',2)],
                   [('X',0), ('X',1), ('X',3), ('X',4)],
                   [('X',4), ('X',5), ('X',7), ('X',8)],
                   [('X',6), ('X',7)],
                   [('Z',0), ('Z',3)],
                   [('Z',1), ('Z',4), ('Z',2), ('Z',5)],
                   [('Z',3), ('Z',6), ('Z',4), ('Z',7)],
                   [('Z',5), ('Z',8)]
                  ]

    ancillae_indexes = [11, 16, 14, 9, 12, 10, 13, 15]

    logical_opers = {'X': [('X',2), ('X',4), ('X', 6)],
                     'Y': [('X',2), ('X',6), ('Z',0), ('Z',8), ('Y',4)],
                     'Z': [('Z',0), ('Z',4), ('Z',8)]
                    }

    stabilizer_CHP = ['+XXIXXIIII',
                      '+IXXIIIIII',
                      '+IIIIXXIXX',
                      '+IIIIIIXXI',
                      '+ZIIIZZIII',
                      '+IZZIZZIII',
                      '+IIIZZZIII',
                      '+IIIIIZZZI',
                      '+IIIIIIZZZ']

    destabilizer_CHP = ['+ZIIIIIIII',
                        '+ZZIIIIIII',
                        '+ZIIIZIIII', 
                        '+IIIIIIZII',
                        '+IIXXIXIXX', 
                        '+IIXIIIIII', 
                        '+IIIXIIIII', 
                        '+IIIIIIIXX', 
                        '+IIIIIIIIX']
    stabilizers_CHP = {
                        'Z': 
                            ['+XXIXXIIII', 
                             '+IXXIIIIII', 
                             '+IIIIXXIXX', 
                             '+IIIIIIXXI', 
                             '+ZIIIZZIII', 
                             '+IZZIZZIII', 
                             '+IIIZZZIII', 
                             '+IIIIIZZZI', 
                             '+IIIIIIZZZ'],
                        'X': ['+XXIXXIIII', '+IXXIIIIII', '+IIXIXIXII', '+IIIIXXIXX', '+IIIIIIXXI', '+ZIIIZIZZI', '+IZZIZZIII', '+IIIZZIZZI', '+IIIIIZIIZ'],
                        'Y': ['+XXIXXIIII', '+IXXIIIIII', '+ZIXIYIXIZ', '+IIIIXXIXX', '+IIIIIIXXI', '+ZIIIZIZZI', '+IZZIZZIII', '+IIIZZIZZI', '+IIIIIZIIZ']
                      }
    
    destabilizers_CHP = {
                        'Z':
                        ['+ZIIIIIIII', 
                         '+ZZIIIIIII', '+ZIIIZIIII', '+IIIIIIZII', '+IIXXIXIXX', '+IIXIIIIII', '+IIIXIIIII', '+IIIIIIIXX', '+IIIIIIIIX'],
                        'X': ['+ZIIIIIIII', '+ZZIIIIIII', '+ZZZIIIIII', '+IZZIZIIII', '+ZZZIIIZII', '+IIIXIIIXI', '+IIIIIXIIX', '+IIIXIIIII', '+IIIIIIIIX'],
                        'Y': ['+ZIIIIIIII', '+ZZIIIIIII', '+ZZZIIIIII', '+IZZIZIIII', '+ZZZIIIZII', '+IIIXIIIXI', '+ZZZIIXIIX', '+IIIXIIIII', '+ZZZIIIIIX']
                      }


    # Xstabs means that we need to insert an Z correction
    # Zstabs means that we need to insert a X correction
    lookuptable = {
        'Xstabs': { '0000': [],
                    '0001': [6],
                    '0010': [5],
                    '0011': [7],
                    '0100': [0],
                    '0101': [3,6],
                    '0110': [4],
                    '0111': [4,6],
                    '1000': [2],
                    #'1001': [1,4,7],  # or [2,6]
                    '1001': [2,6],
                    '1010': [2,5],
                    '1011': [2,7],
                    '1100': [1],
                    '1101': [1,6],
                    '1110': [2,4],
                    '1111': [1,7]
                  },
        'Zstabs': { '0000': [],
                    '0001': [8],
                    '0010': [6],
                    '0011': [7,8],
                    '0100': [1],
                    '0101': [5],
                    '0110': [4],
                    '0111': [4,8],
                    '1000': [0],
                    #'1001': [3,4,5],  # or [0,8]
                    '1001': [0,8],
                    '1010': [3],
                    '1011': [3,8],
                    '1100': [0,1],
                    '1101': [0,5],
                    '1110': [0,4],
                    '1111': [3,5]
                  }
        }
