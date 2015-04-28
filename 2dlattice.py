# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 17:23:50 2015

@author: zubin
"""

import numpy as np
import math
import matplotlib.pyplot as plt



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
T=10




dens=np.zeros((n,9,T+dt),dtype=float) #initial condition zero velocity everywhere??
Angles=np.zeros(9)
e=np.zeros((2,9),dtype=float)
Totdens=np.zeros((n,1,T+1))
AvV_sq=np.zeros((n,T))
Ext=np.ones((1,9))
Tdens=np.zeros((n,9,T))
denspos=np.zeros((n,2,9))
dens_eq=np.zeros((n,9,T))
AvV=np.zeros((n,2,T))
H=np.zeros((n,1))



dens[:,0]=3 #Is this a good pick?

w1=np.array([4, 1, 1, 1, 1, 1, 1, 1, 1],dtype=float)
w2=np.array([9, 36, 9, 36, 9, 36, 9, 36, 9],dtype=float)
w=np.divide(w1,w2)

V=np.array([0, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt])

Angles[1:len(V)]=np.arange(0,2*pi,pi/4)

e[:,1:len(V)]=np.array([[np.cos(Angles[1:len(V)]),np.sin(Angles[1:len(V)])]]) #unit vectors for each of the 9 directions.



for k in range(0,len(V)):
    
    denspos[:,:,k]=R #optimize this step.


for j in range(0,len(V)):
  for i in range(0,n):
    for t in range(0,T):
#         Tdens[i+dt,j,t]=dens[i+dt,j]
         
          
          if denspos[i,1,j]==1:
              
              if j>=6:                  
               denspos[i,:,j]=denspos[i,:,j] + dt*np.transpose(e[:,j])
               denspos[i,:,(j-4)]=denspos[i,:,j]
             
          elif denspos[i,1,j]==(rows-1):
              
               if 1<j<5:
                denspos[i,:,j]=denspos[i,:,j] + dt*np.transpose(e[:,j])
                denspos[i,:,(j+4)]=denspos[i,:,j]
                
          else: denspos[i,:,j]=denspos[i,:,j] + dt*np.transpose(e[:,j])
                
            
        
          Totdens[:,0,:]=dens.sum(axis=1)

          AvV[:,:,t]=np.dot(V*dens[:,:,t],np.transpose(e))/Totdens[:,:,t]
          AvV_sq[:,t]=(AvV[:,:,t]**2).sum(axis=1)
          
          H[:,0]=AvV_sq[:,t]

          dens_eq[:,:,t]=w*Totdens[:,:,t]*(1+3*np.dot(AvV[:,:,t],e)+9/2*np.dot(AvV[:,:,t],e)**2-3/2*np.dot(H,Ext))
         
          dens[:,:,t+dt]=dens[:,:,t]-1/tau*(dens[:,:,t]-dens_eq[:,:,t])
         
         


