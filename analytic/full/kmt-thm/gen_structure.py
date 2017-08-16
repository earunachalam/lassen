#!/usr/bin/python3

import numpy as np
import sys

try:
    ifname = sys.argv[1]
    ofname_adjmx = sys.argv[2]
    ofname_vtxs = sys.argv[3]
except:
    print("Syntax: gen_structure input_filename output_filename_adjacency_matrix output_filename_vertex_locations")

with open(ifname, 'r') as f:
    dat = [d.split() for d in f.readlines()]
    dat = dat[1:]

# convert list of transitions to list of vertices and edges

transitions = [[int(e) for e in d] for d in dat]

vertices    = []
edges       = []

for t in transitions:

    vtx_from    = t[0:2]
    vtx_to      = t[2:4]

    add_from = add_to = True
    idx_from = idx_to = -999

    for ct, v in enumerate(vertices):
        if v == vtx_from:
            add_from = False
            idx_from = ct + 1
        if v == vtx_to:
            add_to = False
            idx_to = ct + 1

    if add_from:
        vertices.append(vtx_from)
        idx_from = len(vertices)

    if add_to:
        vertices.append(vtx_to)
        idx_to = len(vertices)

    edges.append([idx_from, idx_to])


# construct adjacency matrix

nv  = len(vertices)                     # number of vertices
M   = np.zeros((nv, nv), dtype=int)     # adjacency matrix

for e in edges:
    M[e[1]-1,e[0]-1] = 1
    M[e[0]-1,e[1]-1] = 1

# adjust vertex locations
vertices = [[v[1],v[0]] for v in vertices]

# save adjacency matrix and vertex locations
np.savetxt(ofname_adjmx, M, fmt='%4.0f')
np.savetxt(ofname_vtxs, vertices, fmt='%4.0f')
