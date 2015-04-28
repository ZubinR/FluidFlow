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

m=1
dt=1
tau=1 # What should this be???
T=100




dens=np.zeros((n,9),dtype=float) #initial condition zero velocity everywhere??
Angles=np.zeros(9)
e=np.zeros((2,9),dtype=float)
Totdens=np.zeros((n,1))
AvV_sq=np.zeros((n,1))
Ext=np.ones((1,9))
Tdens=np.zeros((n,9,T))
denspos=np.zeros((n,2,9))



dens[:,0]=3E10 #Is this a good pick?

w1=np.array([4, 1, 1, 1, 1, 1, 1, 1, 1],dtype=float)
w2=np.array([9, 36, 9, 36, 9, 36, 9, 36, 9],dtype=float)
w=np.divide(w1,w2)

Angles[1:(dens.shape[1])]=np.arange(0,2*pi,pi/4)

e[:,1:dens.shape[1]]=np.array([[np.cos(Angles[1:(dens.shape[1])]),np.sin(Angles[1:(dens.shape[1])])]]) #unit vectors for each of the 9 directions.

V=np.array([0, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt])

for k in range(0,len(V)):
    
    denspos[:,:,k]=R #optimize this step.

Totdens[:,0]=dens.sum(axis=1)

for j in range(0,len(V)):
  for i in range(0,n):
    for t in range(0,T):
         Tdens[i+dt,j,t]=dens[i+dt,j]
         
         denspos[i,:,j]=denspos[i,:,j] + dt*np.transpose(e[:,j])
        
         AvV=np.dot(V*dens,np.transpose(e))/Totdens
         
         AvV_sq[:,0]=(AvV**2).sum(axis=1)
         
         dens_eq=w*Totdens*(1+3*np.dot(AvV,e)+9/2*np.dot(AvV,e)**2-3/2*np.dot(AvV_sq,Ext))
         
         dens[i+dt,j]=dens[i,j]-1/tau*(dens[i,j]-dens_eq[i,j])
         
         


