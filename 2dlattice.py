import numpy as np
import math
import matplotlib.pyplot as plt

def initLattice(nx,ny):
    i=0
    R=np.zeros((nx*ny,2),dtype=float)
    for x in range(nx):
        for y in range(ny):
             R[i,:]=np.array([[0+x,0+y]]) 
             i+=1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(R[:,0], R[:,1], c='b', marker='o')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    plt.show()
    return R

#### MODEL PARAMETERS##########################################################
nx = 50 ; ny = 20 ; q = 9 ; dt = 1 ; tau = 1 ; m = 1 ; maxiter = 10 ; P0=1
e = np.array([[0,1,1,0,-1,-1,-1,0,1], [0,0,1,1,1,0,-1,-1,-1]]) # unit vectors
weight = np.array([[4./9] , [1./36] , [1./9] , [1./36] , [1./9] , [1./36] , [1./9] , [1./36] , [1./9] ])
##############################################################################
def equilibrium (rho,u):
    for i in range(q):
        eu = np.dot(e.T[i],u.transpose(1,0,2))
        u2 = u[0]**2+u[1]**2
        denseq = np.zeros((q,nx,ny))
        denseq[i,:,:] = rho * weight[i] * (1. + 3. * eu + 0.5*9*eu**2 - 1.5*u2)
        return denseq

u0 = np.zeros((2,nx,ny)) # Initial condition
denseq = equilibrium(1.0,u0) ; densin = denseq.copy()

mask=np.zeros((nx,ny))
mask[:,(0,ny-1)]=1, mask = mask==1
#index = np.transpose(np.nonzero(mask))

for time in range(maxiter):
    rho = np.sum(densin,axis=0)
    u = np.dot(e,densin.transpose(1,0,2))/rho    
    
    
    for j in range(0,q): #Streaming w/ Bounce-Back
           densin[j,:,:] = np.roll(np.roll(densin[j,:,:],e[0,j],axis=0),e[1,j],axis=1)
           
           densin[j,mask] = densin[j+4,mask]
#           u[j,index[0],index[1]]= -1* u[j,index[0],index[1]]
                  
#               if j>=6:                  
##                 densin[j,:,:] = np.roll(np.roll(densin[j,:,:],e[0,j],axis=0),e[1,j],axis=2)
#                 densin[(j-4),:,0]=densin[j,:,:]
#         
#           elif (densin[j,:,1]==(ny-1)*np.ones(nx)).all():
#           
#              if 1<j<5:
##                densin[j,:,:] = np.roll(np.roll(densin[j,:,:],e[0,j],axis=0),e[1,j],axis=1)
#                densin[(j+4),:,ny]=densin[j,:,:]
            
    denseq = equilibrium(rho,u)
    u[0]+=0.05 
    
    for j in range (q): #Relaxation
        densin[j,:,:] = (1-1/tau)*densin[j,:,:] +denseq[j,:,:]