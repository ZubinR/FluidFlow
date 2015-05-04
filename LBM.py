import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

#def initLattice(nx,ny):
#    i=0
#    R=np.zeros((nx*ny,2),dtype=float)
#    for x in range(nx):
#        for y in range(ny):
#             R[i,:]=np.array([[0+x,0+y]]) 
#             i+=1
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.scatter(R[:,0], R[:,1], c='b', marker='o')
#    ax.set_xlabel('X Label')
#    ax.set_ylabel('Y Label')
#    plt.show()
#    return R

#### MODEL PARAMETERS##########################################################
nx = 50 ; ny = 20 ; q = 9 ; dt = 1 ; tau = 1.85 ; m = 1 ; maxiter = 15; du = 0.005
e = np.array([[0,1,1,0,-1,-1,-1,0,1], [0,0,1,1,1,0,-1,-1,-1]]) # unit vectors
weight = np.array([4./9, 1./36, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36, 1./9])
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
    return a*(((ny-1)/2)**2-(y-(ny-1)/2)**2)


u0 = np.zeros((2,nx,ny));ut0=u0[0] # Initial condition
denseq = equilibrium(1.0,u0) 
visc=(2*tau-1)/6
densin = denseq.copy()
mask=np.zeros((nx,ny))

mask[:,(0,ny-1)]=1
mask = mask==1



for time in range(maxiter):    
    for j in range(q): #Streaming w/ Bounce-Back     
        densin[j,:,:] = np.roll(np.roll(densin[j,:,:],e[0,j],axis=0),e[1,j],
                                axis=1)
        if 1<j<5:
            densin[j,mask] = densin[j+4,mask]
        elif j>5:
            densin[j,mask] = densin[j-4,mask]


    rho = np.sum(densin,axis=0)

    u = np.dot(e,densin.transpose(1,0,2))/rho   

    u[0]+=du
    
    print(u[0,1,1])       
    denseq = equilibrium(rho,u)

    
           
#=======    
    for j in range (q): #Relaxation
        densin[j,:,:] = (1-1/tau)*densin[j,:,:] + denseq[j,:,:]/tau
##Steady State condition        
#    u1=u[0]
#    if np.sum(ut1-ut0)/np.sum(ut1)<10:
#        print (maxiter)
#        break
#    else: 
#        ut0=ut1
y=np.arange(ny)  

a,cov=opt.curve_fit(Velx,y,u[0,25,:],0.0,None) 
#Velx=du/(2*visc*nx)*(((ny-1)/2)**2-(y-(ny-1)/2)**2)  #comparison with theory
#a=du/(2*visc*nx)   
plt.plot(y,u[0,25,:], 'r+', y, Velx(y,a)) #plotting velocity profile 
