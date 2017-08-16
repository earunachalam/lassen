#!/usr/bin/python3

# enum_edge.py
# uses edge-by-edge multiplicity data generated by companion script enum_edge.py to draw a config space (coarse-grained by the number of unbound monomers and number of bonds in a given microstate) and annotate each edge with its multiplicity

from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import sys

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

try:
    L = int(sys.argv[1])
    f = open(sys.argv[2])
    imgname = sys.argv[3]
except:
    print('Syntax: plot_edge L input_filename output_filename')
    sys.exit()

print('Reading', sys.argv[1], '...', end=' ')

data = f.read()
print('Done.\nPlotting', imgname, '...', end=' ')
data = np.fromstring(data[data.find('\n'):], dtype=int, sep=' ')
nt = int(np.size(data)/5)   # number of transitions
data = np.reshape(data, (nt, 5))

vtx_radius = 0.1
axmin = np.min(data[:,0:4])
axmax = np.max(data[:,0:4]);
axrange = axmax - axmin
axmin -= 2*vtx_radius
axmax += 2*vtx_radius
scale = 1
dims = (scale*axrange, scale*axrange)

fig = plt.figure(figsize=dims)
ax = fig.add_subplot(111)

def intp(start, end, ratio):
    return start + ratio*(end-start)

def get_vtx_mult(L, x1, y1):
    
    # vertex multiplicities (number of partitions of n into k parts; k incremented along axis 0)
    seqs     = np.ones([5,10])
    seqs[2,] = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    seqs[3,] = [1, 1, 2, 3, 4, 5, 7, 8, 10, 12]
    seqs[4,] = [1, 1, 2, 3, 5, 6, 9, 11, 15, 18]

    k = x1 - y1             # branch
    i = 2*y1 - x1 + 1       # element of sequence within this branch
    g = seqs[k,i-1]         # vertex multiplicity

    return g, i


for it in range(0, nt):
    
    from_nlnk   = data[it,0]
    from_nbnd   = data[it,1]
    to_nlnk     = data[it,2]
    to_nbnd     = data[it,3]
    mult        = data[it,4]
    
    ax.plot([from_nbnd, to_nbnd], [from_nlnk, to_nlnk], 'k', linewidth=2)
    
    circ_from = plt.Circle((from_nbnd, from_nlnk), vtx_radius, color='g')
    circ_to = plt.Circle((to_nbnd, to_nlnk), vtx_radius, color='g')
    ax.add_artist(circ_from)
    ax.add_artist(circ_to)

    content = str(mult) + "\n"

    x1 = from_nbnd
    y1 = from_nlnk
    x2 = to_nbnd
    y2 = to_nlnk
    g, i = get_vtx_mult(x2, x1, y1)
    
    if x1 == x2:
        content += str(int(g*(x1-y1)*(x1-y1-1)/2))
    elif (x2-x1) == 1:
        content += str(int(g*(x1-y1)*(L-x1)))
    elif (x2-x1) == 2:
        content += str(int(g*( 1 + (L-x2)*(0.5*(L+x2-1)-x1) )))
    ax.text(x1-vtx_radius, y1+vtx_radius, str(g), bbox={'facecolor':'red', 'alpha':0.90, 'pad':1})
    ax.text(intp(from_nbnd,to_nbnd,0.7)-vtx_radius, intp(from_nlnk,to_nlnk,0.7), content, bbox={'facecolor':'orange', 'alpha':0.90, 'pad':1})

plt.axis('square')
plt.axis([axmin, axmax, axmax, axmin])
ax.invert_yaxis()

plt.xlabel(r'Number of bound monomers')
plt.ylabel(r'Number of bonds')

plt.grid()

fig.savefig(imgname, bbox_inches='tight')
print('Done.')