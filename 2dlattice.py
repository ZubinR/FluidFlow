import numpy as np
import math
import matplotlib.pyplot as plt
rows=6
columns=6
i=0
n=rows**2 #assume 1 particle at each site? Lattice spacing of 1.
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

#### MODEL PARAMETERS##########################################################
nx = 50 ; ny = 20 ; q = 9 ; dt = 1 ; tau = 1 ; m = 1 ; maxiter = 10
e = np.array([[0,1,1,0,-1,-1,-1,0,1], [0,0,1,1,1,0,-1,-1,-1]]) # unit vectors
weight = np.array([[4./9] , [1./36] , [1./9] , [1./36] , [1./9] , [1./36] , [1./9] , [1./36] , [1./9] ])
##############################################################################
i1 = np.array([1,2,8]) #unknown at left wall
i2 = np.array([0,3,7]) #unknown vertical
i3 = np.array([4,5,6]) #unknown at right wall

def equilibrium (rho,u):
    for i in range(q):
        eu = np.dot(e.T[i],u.transpose(1,0,2))
        u2 = u[0]**2+u[1]**2
        denseq = np.zeros((q,nx,ny))
        denseq[i,:,:] = rho * weight[i] * (1. + 3. * eu + 0.5*9*eu**2 - 1.5*u2)
        return denseq

u0 = np.zeros((2,nx,ny)) # Initial condition
denseq = equilibrium(1.0,u0) ; densin = denseq.copy()

#for time in range(maxiter):
rho = np.sum(densin,axis=0)
u = np.dot(e,densin.transpose(1,0,2))/rho    


























##########Algorithm#########################
#pi=(math.pi)
#sqrt=(math.sqrt)
#
#
#
#dens=np.zeros((columns,rows,q),dtype=float) #initial condition zero velocity everywhere??
#Angles=np.zeros(q)
#e=np.zeros((2,q),dtype=float)
#Totdens=np.zeros((columns,rows))
#AvV_sq=np.zeros((columns,rows,T))
#Ext=np.ones((1,q))
##Tdens=np.zeros((n,9,T))
##denspos=np.zeros((columns,rows,q))
#dens_eq=np.zeros((columns,rows,q))
#AvV=np.zeros((columns,rows,T))
#H=np.zeros((n,1))
#
#
#
#dens[:,0]=3 #Is this a good pick?
#
#
#
##V=np.array([0, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt, 1/dt, sqrt(2)/dt])
#
#Angles[1:q]=np.arange(0,2*pi,pi/4)
# 
#
##for k in range(0,q):
##    
##    denspos[:,:,k]=R #optimize this step.
#
#
#for j in range(0,q):
#    for t in range(0,T):
##         Tdens[i+dt,j,t]=dens[i+dt,j]
#         
#          dens[:,:,j]=dens[e[0,j],e[1,j],j] 
#          
#          if (dens[:,1,j]==np.ones(columns)).all():
#              
#              if j>=6:                  
#               dens[:,:,j]=dens[e[0,j],e[1,j],j] 
#               dens[:,:,(j-4)]=dens[:,:,j]
#             
#          elif (dens[:,1,j]==(rows-1)*np.ones(columns)).all():
#              
#               if 1<j<5:
#                dens[:,:,j]=dens[e[0,j],e[1,j],j] 
#                dens[:,:,(j+4)]=dens[:,:,j]
#                
#          
#                
#            
#        
#          Totdens=dens.sum(axis=2)
#
#          AvV[:,:,t]=np.dot(dens[:,:,j],np.transpose(e))/Totdens
#          AvV_sq[:,t]=(AvV[:,:,t]**2).sum(axis=1)
#          
#          H[:,0]=AvV_sq[:,t]
#
#          dens_eq[:,:,t]=np.transpose(w)*Totdens[:,:,t]*(1+3*np.dot(AvV[:,:,t],e)+9/2*np.dot(AvV[:,:,t],e)**2-3/2*np.dot(H,Ext))
#         
#          dens[:,:,t+dt]=dens[:,:,t]-1/tau*(dens[:,:,t]-dens_eq[:,:,t])
#         
#         


