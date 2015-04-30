# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 17:23:50 2015

@author: zubin
"""

import numpy as np
import math
import matplotlib.pyplot as plt



rows=6
columns=6
i=0
n=rows**2 #assume 1 particle at each site? Lattice spacing of 1.




#vol=n/dens
#L=vol**(1/3)
def initLattice(n,rows):
    i=0
    R=np.zeros((n,2),dtype=float)
    for x in range(0,rows):
        for y in range(0,rows):
             R[i,:]=np.array([[0+x,0+y]]) 
             i+=1
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
q=9



dens=np.zeros((columns,rows,q),dtype=float) #initial condition zero velocity everywhere??
Angles=np.zeros(q)
e=np.zeros((2,q),dtype=float)
Totdens=np.zeros((columns,rows))
AvV_sq=np.zeros((columns,rows))
Ext=np.ones((1,q))
#Tdens=np.zeros((n,9,T))
#denspos=np.zeros((columns,rows,q))
dens_eq=np.zeros((columns,rows,q))
Av=np.zeros((columns,rows,2))
AvV=np.zeros((columns,rows,2))
H=np.zeros((columns,rows))



dens[:,:,3]=3 #Is this a good pick?



#V=np.array([0, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt])

Angles[1:q]=np.arange(0,2*pi,pi/4)

e = np.array([[0,1,1,0,-1,-1,-1,0,1], [0,0,1,1,1,0,-1,-1,-1]]) # unit vectors

w = np.array([[4./9] , [1./36] , [1./9] , [1./36] , [1./9] , [1./36] , [1./9] , [1./36] , [1./9] ])


R=initLattice(n,rows) 

#for k in range(0,q):
#    
#    denspos[:,:,k]=R #optimize this step.



for t in range(0,T):
    for j in range(0,q):
        for i in range(0,columns):
              for k in range(0,rows):
#         Tdens[i+dt,j,t]=dens[i+dt,j]
         
                dens[i+e[0,j],k+e[1,j],j] = dens[i,k,j]
          
                if (dens[:,1,j]==np.ones(columns)).all():
              
                   if j>=6:                  
                     dens[i+e[0,j],k+e[1,j],j] = dens[i,k,j]
                     dens[:,:,(j-4)]=dens[:,:,j]
             
                elif (dens[:,1,j]==(rows-1)*np.ones(columns)).all():
              
                  if 1<j<5:
                    dens[i+e[0,j],k+e[1,j],j] = dens[i,k,j]
                    dens[:,:,(j+4)]=dens[:,:,j]
                
          
                
            
        
                  Totdens=dens.sum(axis=2)
          
                  Av=dens[i,k,j]*e[:,j]/Totdens[i,k]
                  AvV[i,k,:]=Av+AvV[i,k,:]
#          AvV[:,:,0]=np.dot(dens,np.transpose(e))/Totdens
#          AvV[:,:,1]=np.dot(dens,np.transpose(e))/Totdens
                  AvV_sq[:,:]=(AvV[:,:,:]**2).sum(axis=2)
          
#          H[:,0]=AvV_sq[:,:]

#          dens_eq[:,:,:]=np.transpose(w)*Totdens[:,:]*(1+3*np.dot(AvV[:,:,:],e)+9/2*np.dot(AvV[:,:,:],e)**2-3/2*np.dot(AvV_sq,Ext))
#         
#          dens[:,:,:]=dens[:,:,t]-1/tau*(dens[:,:,t]-dens_eq[:,:,t])
         
         


