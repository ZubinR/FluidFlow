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
mask[25:25+lx,10:10+ly] =False
mask[:,[0,ny-1]]=False
mask2=~mask







up=np.arange(2,5)
dowt
n=np.arange(6,9)
mid =np.array([0,1,5])

# Main time loop
for time in range(maxiter):    
    for j in range(q): #Streaming of the bulk.
        densbulk = densin[j,mask].reshape(nx,ny-2) 
        densbulk = np.roll(np.roll(densbulk,e[0,j],axis=0),e[1,j],
                                axis=1)       
        densin[j,:,1:-1]=densbulk

    densin[up,:,1]=densin[up,:,1]+densin[down,:,1] # Assigns density opposite velocity at pre-boundary nodes.
    densin[down,:,ny-2]=densin[down,:,ny-2]+densin[up,:,ny-2]      
    
    densin[down,:,9]=densin[down,:,9]+densin[up,:,9]
    densin[up,:,10+ly+1]=densin[up,:,10+ly+1]+densin[down,:,10+ly+1]
    densin[mid[2],24,:]=densin[mid[2],24,:]+densin[mid[1],24,:]
    densin[mid[1],25+lx+1,:]=densin[mid[1],25+lx+1,:]+densin[mid[2],25+lx+1,:]
    
    rho = np.sum(densin,axis=0)
    u = np.dot(e,densin.transpose(1,0,2))/rho  
    u[0,mask]+=du    
#    print(rho[1,1]) # Checks density conservation
    denseq = equilibrium(rho,u)
 
    for j in range (q): #Relaxation
        densin[j,:,1:-1] = (1-1/tau)*densin[j,:,1:-1] + denseq[j,:,1:-1]/tau
          
    densin[down,:,1]=0 # Forces densities with velocities away from the bulk to zero.
    densin[up,:,ny-2]=0
    densin[up,:,9]=0
    densin[down,:,10+ly+1]=0
    densin[mid[2],25+lx+1,:]=0
    densin[mid[1],24,:]=0
    
Re =np.sum(u)*ny/(nx*ny*visc)
#plotting velocity profile 
U=u[0,25,:]
y=np.arange(ny)  
a,cov=opt.curve_fit(Velx,y,U,0.0,None)    
plt.plot(U,y, 'r+', Velx(y,a),y) 
plt.xlabel('Horizontal Velocity') ; plt.ylabel('y') ;plt.title( 'Poisseuile Flow')
curve=curvature(U)
