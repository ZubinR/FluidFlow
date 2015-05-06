# -*- coding: utf-8 -*-
"""
Created on Tue May  5 14:36:59 2015

@author: zubin
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

#### MODEL PARAMETERS##########################################################
nx = 50 ; ny = 20 ; lx=3; ly=5 ; q = 9 ; dt = 1 ; tau = 1.85 ; maxiter = 3000; du = 0.005
e = np.array([[0,1,1,0,-1,-1,-1,0,1], [0,0,1,1,1,0,-1,-1,-1]]) # unit vectors
weight = np.array([4./9, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36])
###############################################################################
def equilibrium (rho,u):
    denseq = np.zeros((q,nx,ny))
    eu = np.dot(e.T,u.transpose(1,0,2))
    u2 = u[0]**2+u[1]**2
    for i in range(q):                      
        denseq[i,:,:] = rho * weight[i] * (1. + 3. * eu[i] + 0.5*9*eu[i]**2 - 
                        1.5*u2)
    return denseq
   
   
def Velx(y,a):  
    #a=dp/(2*visc*nx)
    return a*y*(y-(ny-1))


def curvature(f):
    dx = 1
    fdiff = np.zeros(len(f)-2)
    for i in range (1,len(f)-1):    
        fdiff[i-1] = dx**(-2) *(f[i-1] + f[i+1] -2*f[i])
    return fdiff
    
visc=(2*tau-1)/6
u0 = np.zeros((2,nx,ny)) # Initial condition
rho0=np.ones((nx,ny))
denseq = equilibrium(rho0,u0) 
densin = denseq.copy()
mask = np.ones((nx,ny), dtype=bool)
mask2 = np.ones((nx,ny), dtype=bool)
#mask[25:25+lx,10:10+ly] =False
mask[:,[0,ny-1]]=False
mask2[25:25+lx,10:10+ly] =False
mask2[:,[0,ny-1]]=False

up=np.arange(2,5)
down=np.arange(6,9)
left =np.array([4,5,6])
right=np.array([1,2,8])

# Main time loop
for time in range(maxiter): 
    densupbound=densin[down,:,1]
    denslowbound=densin[up,:,ny-2]
#    densin[:,~mask2]=0

    #Define densities which point towards the obstacle before streaming step
    densleftboundobst=densin[right,24,10:10+ly] 
    densrightboundobst=densin[left,25+lx+1,10:10+ly]
    densupboundobst=densin[down,25:25+lx,10+ly+1]
    denslowboundobst=densin[up,25:25+lx,9]    
    densin[down,:,1]=0 # Forces densities with velocities away from the bulk to zero.
    densin[up,:,ny-2]=0
    #Forces densities with velocities which point towards the obstacle to zero.
    densin[up,25:25+lx,9]=0
    densin[down,25:25+lx,10+ly+1]=0
    densin[left,25+lx+1,10:10+ly]=0
    densin[right,24,10:10+ly]=0
        
    
    
    
    
    for j in range(q): #Streaming of the bulk.
        densbulk = densin[j,mask].reshape(nx,ny-2) 
        densbulk = np.roll(np.roll(densbulk,e[0,j],axis=0),e[1,j],
                                axis=1)       
        densin[j,:,1:-1]=densbulk

    densin[up,:,1]=densin[up,:,1]+densupbound # Assigns density opposite velocity at pre-boundary nodes.
    densin[down,:,ny-2]=densin[down,:,ny-2]+denslowbound      
    # Assigns density opposite velocity at pre-boundary nodes for the obstacle.
    densin[down,25:25+lx,9]=densin[down,25:25+lx,9]+denslowboundobst
    densin[up,25:25+lx,10+ly+1]=densin[up,25:25+lx,10+ly+1]+densupboundobst
    densin[left,24,10:10+ly]=densin[left,24,10:10+ly]+densleftboundobst
    densin[right,25+lx+1,10:10+ly]=densin[right,25+lx+1,10:10+ly]+densrightboundobst
    
    rho = np.sum(densin,axis=0)
    u = np.dot(e,densin.transpose(1,0,2))/rho  
    u[0,mask2]+=du    
#    print(rho[1,1]) # Checks density conservation
    denseq = equilibrium(rho,u)
 
    for j in range (q): #Relaxation
        densin[j,mask2] = (1-1/tau)*densin[j,mask2] + denseq[j,mask2]/tau
          
    
    
Re =np.sum(u)*ny/(nx*ny*visc)
##plotting velocity profile 
U=u[0,25,:]
y=np.arange(ny)  
a,cov=opt.curve_fit(Velx,y,U,0.0,None) 
#plt.plot(U,y)   
plt.plot(U,y, 'r+', Velx(y,a),y) 
plt.xlabel('Horizontal Velocity') ; plt.ylabel('y') ;plt.title(  'Poiseuille Flow, '"Re=%g"%Re)
curve=curvature(U)
