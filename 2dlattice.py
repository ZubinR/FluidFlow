# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 17:23:50 2015

@author: zubin
"""

import numpy as np
import math
import matplotlib.pyplot as plt
#T=300


rows=6
i=0
n=rows**2 #assume 1 particle at each site? Lattice spacing of 1.




#vol=n/dens
#L=vol**(1/3)
R=np.zeros((n,2),dtype=float)
#Rf=np.zeros((n,3),dtype=float)
for x in range(0,rows):
    for y in range(0,rows):
         R[i,:]=np.array([[0+x,0+y]]) 
         i=i+1

#Rf=R*L/rows
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(R[:,0], R[:,1], c='b', marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')

plt.show()

#########Algorithm#########################
pi=(math.pi)
sqrt=(math.sqrt)

dens=np.zeros(9) #initial condition zero velocity everywhere??
Angles=np.zeros(9)
e=np.zeros((2,9),dtype=float)
m=1
dt=1

dens[0]=3E10 #Is this a good pick?

w1=np.array([4, 1, 1, 1, 1, 1, 1, 1, 1],dtype=float)
w2=np.array([9, 36, 9, 36, 9, 36, 9, 36, 9],dtype=float)
w=np.divide(w1,w2)

Angles[1:(len(dens))]=np.arange(0,2*pi,pi/4)

e[:,1:len(dens)]=np.array([[np.cos(Angles[1:(len(dens))]),np.sin(Angles[1:(len(dens))])]]) #unit vectors for each of the 9 directions.

V=np.array([0, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt])

AvV=np.dot(V*dens,np.transpose(e))
dens_eq=w*sum(dens)*(1+3*np.dot(AvV,e)+9/2*np.dot(AvV,e)**2-3/2*np.dot(AvV,np.transpose(AvV)))
