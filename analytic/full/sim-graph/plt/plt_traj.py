#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys

traj = np.loadtxt('traj_4.dat')

t = traj[:,0]
t1 = traj[:,1]
t2 = traj[:,2]
t3 = traj[:,3]
t4 = traj[:,4]
t5 = traj[:,5]

for var in range(1,6):
    leg = str(var)
    plt.plot(t,traj[:,var],label=leg)

print(traj[-1,:])
print(traj[-1,1:])

# without first vertex

ptot_no_t1 = np.sum(traj[:,2:6], axis=1)
ptot_no_t1_tile = np.transpose(np.tile(ptot_no_t1,(4,1)))
traj_no_t1 = np.divide(traj[:,2:6],ptot_no_t1_tile)

# print(traj_no_t1[-1,:])
# for var in range(0,4):
    # plt.plot(t,traj_no_t1[:,var])

plt.legend()
plt.show()
