import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from pylab import *
from matplotlib import animation
#### MODEL PARAMETERS##########################################################
nx = 50 ; ny = 20 ; q = 9 ; tau = 1.85 ; maxiter = 100; du = 0.005
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
uobj=np.zeros((2,nx,ny))
Ux=np.zeros((maxiter+1,nx,ny))
Uy=np.zeros((maxiter+1,nx,ny))
S=np.zeros((2,4*lx*ly))
Fy=np.zeros((q,nx,ny))
Fx=np.zeros((q,nx,ny))
Ftot=np.zeros((2))
Ftot0=np.zeros((2))
i=0

mask = np.ones((nx,ny),dtype=bool)
mask[:,[0,-1]]=False
mask3= np.zeros((nx,ny),dtype=bool)
mask3[obx-lx:obx+lx,oby-ly:oby+ly] = True
mask[obx-lx:obx+lx,oby-ly:oby+ly] = False
objmask = np.zeros((nx,ny),dtype=bool)
objmask[obx-lx:obx+lx,[oby-ly,oby+ly-1]]=True   ; objmask[[obx-lx,obx+lx-1],oby-ly:oby+ly]=True
notbulk= ~mask
wall=np.zeros((q,nx,ny),dtype=bool)
edge = np.zeros((q,nx,ny),dtype=bool)
edgefluid=np.zeros((q,nx,ny),dtype=bool) # Will contain the crossed boundary points.
qflip = np.mod((np.arange(q) +3),8)+1 ; qflip[0]=0
Rcom=np.array([obx,oby])
utemp = np.zeros((2,nx,3*ny))

for time in range(maxiter):
    for j in range(q):
        wall[j,:,:] = np.logical_and(notbulk,
                                     np.roll(np.roll(mask,e[j,0],axis=0),e[j,1],axis=1))
        
        edge [j,:,:]= np.logical_and(mask3,
                                     np.roll(np.roll(~mask3,e[j,0],axis=0),e[j,1],axis=1))
        edgefluid[j,:,:] = np.roll(np.roll(edge[j,:,:],e[qflip[j],0],axis=0),e[qflip[j],1],axis=1)
    
    index = []
    edgeindex = [] # Contains indices of boundary points adjacent to bulk points.
    edgefluidindex = []
    for j in range(q):
        [x,y]=wall[j,:,:].nonzero()
        index.append([x,y])
        [m,n]=edge[j,:,:].nonzero() 
        edgeindex.append([m,n])    
        [k,l]=edgefluid[j,:,:].nonzero() 
        edgefluidindex.append([k,l])
    


    ub[:,mask3] = u[:,mask3]    
    eub = np.dot(e,ub.transpose(1,0,2))
#    print(densnew[5,2,2]-densold[5,2,2])
    for j in range(q): 
        densnew[j,:,:]=np.roll(np.roll(densold[j,:,:],e[j,0],axis=0),e[j,1],axis=1)
    for j in range(q):

        densnew[j,index[j][0],index[j][1]] = densnew[qflip[j],index[j][0],
                                                   index[j][1]] - 6 * weight[j]*rho[index[j][0],index[j][1]]*eub[j,index[j][0],index[j][1]]

        densnew[j,edgefluidindex[j][0],edgefluidindex[j][1]] = densnew[qflip[j],edgefluidindex[j][0],
                                                   edgefluidindex[j][1]] - 6 * weight[j]*rho[edgefluidindex[j][0],edgefluidindex[j][1]]*eub[j,edgefluidindex[j][0],edgefluidindex[j][1]]
        
        Fx[j,edgeindex[j][0],edgeindex[j][1]] = 2*(densold[j,edgefluidindex[j][0],edgefluidindex[j][1]] - densnew[[j,edgeindex[j][0],edgeindex[j][1]]]-2*(weight[j]*3*rho[edgeindex[j][0],edgeindex[j][1]]*eub[j,edgeindex[j][0],edgeindex[j][1]]))*e[j,0]
        Fy[j,edgeindex[j][0],edgeindex[j][1]] = 2*(densold[j,edgefluidindex[j][0],edgefluidindex[j][1]] - densnew[[j,edgeindex[j][0],edgeindex[j][1]]]-2*(weight[j]*3*rho[edgeindex[j][0],edgeindex[j][1]]*eub[j,edgeindex[j][0],edgeindex[j][1]]))*e[j,1]
    
    Ftot[0] = sum (Fx) ; Ftot[1] = sum(Fy)
    F = 0.5 * (Ftot+Ftot0)    
    Ftot0 = Ftot
    
#    F = np.array([0.01,0]) 
    rho = np.sum(densnew,axis=0)
    u = np.dot(e.T,densnew.transpose(1,0,2))/rho  
    u[0,mask]+=du
   
##    print(rho[1,1]) # Checks density conservation
    denseq = equilibrium(rho,u)
    for j in range (q): #Relaxation
        densnew[j,mask] = (1-1/tau)*densnew[j,mask] + denseq[j,mask]/tau
#    uold[:,mask3] = u[:,mask3]
#    unew[0,mask3]=uold[0,mask3]+2*F[0]
#    unew[1,mask3]=uold[1,mask3]+2*F[1]
#    uobj[:,mask3]=unew[:,mask3]
#    u[:,mask3]=uobj[:,mask3]
    utemp[:,:,20:40] = u
    utemp[0,:,40:60] = utemp[0,:,0:20] + 2*F[0]*mask3
    utemp[1,:,40:60] = utemp[1,:,0:20] + 2*F[1]*mask3
    utemp[:,:,0:20] = 0 ; utemp = np.roll(utemp,-20,axis=2)
    temp = utemp[:,:,20:40] ; u [:,mask3] = temp[:,mask3]    
#    u = utemp[:,:,20:40]
# Here we try to move the masks and the center of mass by 1 point to the right if the velocity is large enough.   
    S+=u[:,mask3]
    if S[0,0]>=i+1:
        i=i+1
        objmask=np.roll(objmask,1,axis=0)
        mask=np.roll(mask,1,axis=0)
        mask3=np.roll(mask3,1,axis=0)
        uobj=np.roll(uobj,1,axis=0) 
        uobj[:,~mask3]=0
        Rcom+=e[1,:]
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
v = np.transpose(u)
fig, ax = plt.subplots(1,1)
plt.imshow(v[1:-1,:,0]**2+v[1:-1,:,1]**2)
plt.colorbar()
Q = ax.quiver(v[1:-1,:,0], v[1:-1,:,1])
ax.set_xlim(0, nx-1)
ax.set_ylim(0, ny-2)
show()

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
