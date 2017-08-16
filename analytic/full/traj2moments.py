#!/usr/bin/python3

# traj2moments.py
# calculates moments of segment length distribution from state populations (loaded from trajectory file)

from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npm
import sys

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

try:
    data    = np.loadtxt(sys.argv[1])
except:
    print('Syntax: traj2moments.py input_filename')
    sys.exit()

vtx_loc     = data[0:2,1:]                      # vertex locations on graph
traj        = data[2:,1:]                       # trajectory (i.e. p(vertex) over time)
t           = data[2:,0]                        # time
ntstep      = np.size(t)                        # number of timesteps

L           = np.max(vtx_loc)                   # number of monomer units in the system

vtx_nbound  = vtx_loc[0,]                       # number of bound monomers in config at each vertex
vtx_nbonds  = vtx_loc[1,]                       # number of bonds in config at each vertex
vtx_nmon    = L - vtx_nbound                    # number of monomers in config at each vertex
vtx_nseg    = L - vtx_nbonds                    # number of segments in config at each vertex
vtx_lbar    = np.divide(L, vtx_nseg)            # average segment length in config at each vertex

moments = np.empty([ntstep,0])
moment0 = np.transpose(np.sum(np.multiply(traj,npm.repmat(np.ones(np.size(vtx_lbar)),ntstep,1)),axis=1))
moment1 = np.sum(np.multiply(traj,npm.repmat(np.power(vtx_lbar,1.0),ntstep,1)),axis=1)
moment2 = np.sum(np.multiply(traj,npm.repmat(np.power(vtx_lbar,2.0),ntstep,1)),axis=1)
moment3 = np.sum(np.multiply(traj,npm.repmat(np.power(vtx_lbar,3.0),ntstep,1)),axis=1)
moments = np.concatenate((t[:,None],moment0[:,None],moment1[:,None],moment2[:,None],moment3[:,None]),axis=1)
print('\n'.join(' '.join(str(cell) for cell in row) for row in moments))
