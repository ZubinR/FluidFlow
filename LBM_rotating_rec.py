import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from pylab import *
from matplotlib import animation
#### MODEL PARAMETERS##########################################################
nx = 50 ; ny = 20 ; q = 9 ; tau = 1.85 ; maxiter = 4; du = 0.005
obx= nx/2; oby = ny/2  #coordinates of bottom left rectangular obstacle point
lx=3; ly=5 
e = np.array([[0,0], [1,0], [1,1], [0,1], [-1,1], [-1,0], [-1,-1], [0,-1], [1,-1]])
weight = np.array([4./9, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36, 1./9, 1./36])
###############################################################################
def equilibrium (rho,u):
    denseq = np.zeros((q,nx,ny))
    eu = np.dot(e,u.transpose(1,0,2))
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
    return fdiff [1]  

#Initializing parameters/matrices
visc=(2*tau-1)/6
u = np.zeros((2,nx,ny))
uold=np.zeros((2,nx,ny))
unew=np.zeros((2,nx,ny))
ub = np.zeros((2,nx,ny))
rho=np.ones((nx,ny))
denseq = equilibrium(rho,u) 
densold = denseq.copy()
densnew = np.zeros((q,nx,ny))


Ux=np.zeros((maxiter+1,nx,ny))
Uy=np.zeros((maxiter+1,nx,ny))


Fy=Fx=np.zeros((q,nx,ny))
Ftot=np.zeros((2))
Ftot0=np.zeros((2))
mask = np.ones((nx,ny),dtype=bool)
mask[:,[0,-1]]=False
mask3= np.zeros((nx,ny),dtype=bool)
mask3[obx-lx:obx+lx,oby-ly:oby+ly] = True
mask[obx-lx:obx+lx,oby-ly:oby+ly] = False
objmask = np.zeros((nx,ny),dtype=bool)
objmask[obx-lx:obx+lx,[oby-ly,oby+ly-1]]=True   ; objmask[[obx-lx,obx+lx-1],oby-ly:oby+ly]=True
notbulk= ~mask
wall = np.zeros((q,nx,ny),dtype=bool) # Will contain the crossed boundary points.

for j in range(q):
    wall[j,:,:] = np.logical_and(notbulk,
                                 np.roll(np.roll(mask,e[j,0],axis=0),e[j,1],axis=1))
    
qflip = np.mod((np.arange(q) +3),8)+1 ; qflip[0]=0


index = [] # Contains indices of boundary points adjacent to bulk points.
for j in range(q):
    [x,y]=wall[j,:,:].nonzero()
    index.append([x,y])
        
Rcom=np.array([obx,oby])
for time in range(maxiter):
#    print(time) 
    uobj[:,mask3]=u[:,mask3]
    uold[:,mask3]=uobj[:,mask3]
    ub[:,objmask]=u[:,objmask]
    print(ub[:,objmask])
    eub = np.dot(e,ub.transpose(1,0,2))
#    print(densnew[5,2,2]-densold[5,2,2])
    for j in range(q): 
        densnew[j,:,:]=np.roll(np.roll(densold[j,:,:],e[j,0],axis=0),e[j,1],axis=1)
#    print(densnew[5,4,2]-densold[5,4,2])
    for j in range(q):
        densnew[j,index[j][0],index[j][1]] = densnew[qflip[j],index[j][0],
                                                   index[j][1]] - 6 * weight[j]*rho[index[j][0],index[j][1]]*eub[j,index[j][0],index[j][1]]
#    
    rho = np.sum(densnew,axis=0)
    u = np.dot(e.T,densnew.transpose(1,0,2))/rho  
    u[0,mask]+=du
    
#
#
#   
#    Ux[time,:,:]=u[0,:,:]
#    Uy[time,:,:]=u[1,:,:]
##    print(rho[1,1]) # Checks density conservation
    denseq = equilibrium(rho,u)
    for j in range (q): #Relaxation
        densnew[j,mask] = (1-1/tau)*densnew[j,mask] + denseq[j,mask]/tau
#    
        Fx[j,objmask]=2*(densold[j,objmask]- densnew[j,objmask] - 2*(weight[j]*3*rho[objmask]*eub[j,objmask]))*e[j,0]
        Fy[j,objmask]=2*(densold[j,objmask]- densnew[j,objmask] - 2*(weight[j]*3*rho[objmask]*eub[j,objmask]))*e[j,1]
#    
    
    Ftot[0]=sum(Fx); Ftot[1]=sum(Fy)
    F = np.array([0.01,0])#0.5 * (Ftot+Ftot0)
    

    
    unew[0,mask3]=uold[0,mask3]+2*F[0]
    unew[1,mask3]=uold[1,mask3]+2*F[1]
    u[:,mask3]=unew[:,mask3]
#    Ftot0=Ftot
#    Ux[time+1,objmask]=2*Ftot[0] #+Ux[time-1,objmask]
#    Uy[time+1,objmask]=2*Ftot[1] #Uy[time-1,objmask]+
##    
#    u[0,objmask]=Ux[time+1,objmask]
#    u[1,objmask]=Uy[time+1,objmask]
#    
    Rcom+=u[:,obx,oby]
#    
    objmask=np.roll(np.roll(objmask,int(u[0,obx,oby]),axis=0),int(u[1,obx,oby]),axis=1)
    mask3=np.roll(np.roll(mask3,int(u[0,obx,oby]),axis=0),int(u[1,obx,oby]),axis=1)
    uobj=np.roll(np.roll(uobj,int(u[0,obx,oby]),axis=0),int(u[1,obx,oby]),axis=1)     
    densold=densnew
    
#%%
####Plotting velocity profiles#####        
#U=u[0,obx,:]
#Re =np.sum(u)*ny/(nx*ny*visc)
#y=np.arange(ny)
#a,cov=opt.curve_fit(Velx,y,U,0.0,None) 
#plt.figure()
#plt.plot(U,y, 'r+', Velx(y,a),y) 
#plt.xlabel('Horizontal Velocity') ; plt.ylabel('y') ;plt.title( 'Poiseuille Flow, '"Re=%g"%Re)
#curve=curvature(U)
#plt.show()
#
#       
#v = np.transpose(u)
#fig, ax = plt.subplots(1,1)
#plt.imshow(v[1:-1,:,0]**2+v[1:-1,:,1]**2)
#plt.colorbar()
#Q = ax.quiver(v[1:-1,:,0], v[1:-1,:,1])
#ax.set_xlim(0, nx-1)
#ax.set_ylim(0, ny-2)
#show()

##########Animation##########
#fig, ax = plt.subplots(1,1)
#Vinx=np.transpose(Ux[0,:,:]) 
#Viny=np.transpose(Uy[0,:,:])
#im=ax.imshow(Vinx[1:-1,:]**2+Viny[1:-1,:]**2)
#
#def animate(i,im,Ux,Uy):
#    
#    Vinx=np.transpose(Ux[i,:,:])
#    Viny=np.transpose(Uy[i,:,:])
#    im.set_array(Vinx[1:-1,:]**2+Viny[1:-1,:]**2)
#    im.autoscale()    
#    return im
#    
#anim = animation.FuncAnimation(fig, animate, np.arange(maxiter),fargs=(im,Ux, Uy),
#                               interval=100, blit=False, repeat=True)    
#title('velocity profile rectangular obstacle')
#show() 
# 
