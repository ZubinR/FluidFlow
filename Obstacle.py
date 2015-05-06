# -*- coding: utf-8 -*-
"""
Created on Tue May  5 14:36:59 2015

@author: zubin
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from pylab import *

#### MODEL PARAMETERS##########################################################
nx = 50 ; ny = 20
obx= nx/2; oby = ny/2  #coordinates of bottom left rectangular obstacle point
lx=3; ly=5 
q = 9 
dt = 1 ; tau = 1.85
maxiter = 3000
du = 0.005
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

up=np.arange(2,5)
down=np.arange(6,9)
left =np.array([4,5,6])
right=np.array([1,2,8])
u = np.ndarray((nx,ny),dtype=float)
# Main time loop
for time in range(maxiter): 
 
    for j in range(q): #Streaming of the bulk.
        densbulk = densin[j,:,:]
        densbulk = np.roll(np.roll(densbulk,e[0,j],axis=0),e[1,j],
                                axis=1)       
        densin[j,:,:] = densbulk
    densin[up,:,0] = densin[down,:,0]
    densin[down, :, ny-1] = densin[up, :, ny-1]
    densin[left,obx,oby:oby+ly]=densin[right,obx,oby:oby+ly]
    densin[right,obx+lx,oby:oby+ly]=densin[left,obx+lx,oby:oby+ly]
    densin[up,obx:obx+lx,oby+ly]=densin[down,obx:obx+lx,oby+ly]
    densin[down,obx:obx+lx,oby]=densin[up,obx:obx+lx,oby]
    
    rho = np.sum(densin,axis=0)
    u = np.dot(e,densin.transpose(1,0,2))/rho  
#    print 'a', u[0,:,1:ny-2]
#    print 'b', u[1,:,1:ny-2]
#    raw_input()
    u[0,:,:]+=du    
#    print(rho[1,1]) # Checks density conservation
    denseq = equilibrium(rho,u)
 
    for j in range (q): #Relaxation
        densin[j,:,1:-1] = (1-1/tau)*densin[j,:,1:-1] + denseq[j,:,1:-1]/tau
                 
#tempvel = np.copy(u)
# tempvel[outside] = [0.0, 0.0]
v = np.transpose(u)
Q = quiver(v[1:-1,:,0], v[1:-1,:,1])
l,r,b,t = axis()
dx, dy = r-l, t-b
axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])
title('velocity profile rectangular obstacle')
show()    
Re =np.sum(u)*ny/(nx*ny*visc)
##plotting velocity profile 
U=u[0,obx,:]
y=np.arange(ny)  
a,cov=opt.curve_fit(Velx,y,U,0.0,None)   
plt.plot(U,y, 'r+', Velx(y,a),y) 
plt.xlabel('Horizontal Velocity') ; plt.ylabel('y') ;plt.title(  'Poiseuille Flow, '"Re=%g"%Re)
curve=curvature(U)
