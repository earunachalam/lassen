#!/usr/bin/python3

# enum_edge.py
# count multiplicities of edges for a config space coarse-grained by the number of unbound monomers and number of segments in a given microstate. Results are printed edge-by-edge in column format: #seg in origin state, #mon in origin state, #seg in target state, #mon in target state, multiplicity of said edge
# results can be plotted using companion script plot_edge.py

import math
import numpy as np
import sys

# total number of monomer units in system
try:
    L = int(sys.argv[1])
except:
    print('Syntax: enum_edge L \nwhere L = total number of monomer units in system')
    sys.exit()

# list of dictionaries containing information about each transition
transitions = []

from_state = np.ones(L, dtype=int)
final_state = np.array([L], dtype=int)
from_state_queue = []
from_state_past = []

# from_ustate is the 'current' ustate
while not np.array_equal(from_state, final_state):
    
    from_nseg = from_state.size
    from_nlnk = L - from_nseg
    from_nmon = np.count_nonzero(from_state == 1)
    from_nbnd = L - from_nmon

    for elem_idx, elem in enumerate(from_state):
        for elem2_idx, elem2 in enumerate(from_state):
            if elem_idx >= elem2_idx:
                continue
            
            to_state = np.delete(from_state, [elem_idx, elem2_idx])
            to_state = np.append(to_state, elem+elem2)
            to_state = np.sort(to_state)

            appendto_from_state_queue = True
            for already in from_state_past:
                if np.array_equal(already, to_state):
                    appendto_from_state_queue = False
            for already in from_state_queue:
                if np.array_equal(already, to_state):
                    appendto_from_state_queue = False
            if appendto_from_state_queue:
                from_state_queue.append(to_state)
            
            to_nseg = to_state.size
            to_nlnk = L - to_nseg
            to_nmon = np.count_nonzero(to_state == 1)
            to_nbnd = L - to_nmon
            
            create_new = True
            for existing in transitions:
                if (from_nlnk == existing['from_nlnk']) and (from_nbnd == existing['from_nbnd']) and (to_nlnk == existing['to_nlnk']) and (to_nbnd == existing['to_nbnd']):
                    existing['mult'] += 1
                    create_new = False
            if create_new:
                transitions.append({'from_nlnk': from_nlnk, 'from_nbnd': from_nbnd, 'to_nlnk': to_nlnk, 'to_nbnd': to_nbnd, 'mult': 1})

    from_state_past.append(from_state)
    from_state = from_state_queue.pop(0)

print('from_nlnk from_nbnd to_nlnk to_nbnd mult')
for t in transitions:
    print(t['from_nlnk'],' ', t['from_nbnd'], ' ', t['to_nlnk'], ' ', t['to_nbnd'], ' ', t['mult'])

