# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:06:04 2019

@author: Boda Li
"""

import numpy as np
import time
from numba import jit



xm=1.05e7 #in m
ym=4.5e6
eps=0.015#boundary current width and current vlocity
nx=int(np.ceil(10/eps)); ny=int(np.ceil(nx*ym/xm))
dx=xm/nx; dy=ym/ny
x=np.linspace(0,xm,nx+1,endpoint=True)
y=np.linspace(-ym,ym,2*ny+1,endpoint=True)
X,Y=np.meshgrid(x,y)
lam1=(-1+(1+4*np.pi**2*eps**2)**0.5)/2/eps;
lam2=(-1-(1+4*np.pi**2*eps**2)**0.5)/2/eps;
c1=(1-np.exp(lam2))/(np.exp(lam2)-np.exp(lam1));
c2=-(1+c1);
A= 5.025e4;#m^2/s
Psi=A*(c1*np.exp(lam1*X/xm)+c2*np.exp(lam2*X/xm)+1)*np.sin(np.pi*Y/ym);
Actual=Psi.max()
A=A*A/Actual
Psi=A*(c1*np.exp(lam1*X/xm)+c2*np.exp(lam2*X/xm)+1)*np.sin(np.pi*Y/ym);

x=np.linspace(0.5*dx,xm-0.5*dx,nx,endpoint=True)
y=np.linspace(-ym+0.5*dy,ym-0.5*dy,2*ny,endpoint=True)
X,Y=np.meshgrid(x,y)
u=-(np.pi*A/ym)*(c1*np.exp(lam1*X/xm)+c2*np.exp(lam2*X/xm)+1)*np.cos(np.pi*Y/ym)
v=(A/xm)*(c1*lam1*np.exp(lam1*X/xm)+c2*lam2*np.exp(lam2*X/xm))*np.sin(np.pi*Y/ym)
    
#Velocity for the boundary for advection (the reason for second upwind differencing method
xr=np.linspace(dx,xm,nx,endpoint=True)
xl=np.linspace(0,xm-dx,nx,endpoint=True)
ya=np.linspace(-ym+dy,ym,2*ny,endpoint=True)
yb=np.linspace(-ym,ym-dy,2*ny,endpoint=True)

Xo,Yo=np.meshgrid(xr,y)
ur=-(np.pi*A/ym)*(c1*np.exp(lam1*Xo/xm)+c2*np.exp(lam2*Xo/xm)+1)*np.cos(np.pi*Yo/ym)
Xo,Yo=np.meshgrid(xl,y)
ul=-(np.pi*A/ym)*(c1*np.exp(lam1*Xo/xm)+c2*np.exp(lam2*Xo/xm)+1)*np.cos(np.pi*Yo/ym)
Xo,Yo=np.meshgrid(x,ya)
va=(A/xm)*(c1*lam1*np.exp(lam1*Xo/xm)+c2*lam2*np.exp(lam2*Xo/xm))*np.sin(np.pi*Yo/ym)
Xo,Yo=np.meshgrid(x,yb)
vb=(A/xm)*(c1*lam1*np.exp(lam1*Xo/xm)+c2*lam2*np.exp(lam2*Xo/xm))*np.sin(np.pi*Yo/ym)
x=np.linspace(0.5*dx,xm-0.5*dx,nx,endpoint=True)
y=np.linspace(-ym+0.5*dy,ym-0.5*dy,2*ny,endpoint=True)
X,Y=np.meshgrid(x,y)
u=-(np.pi*A/ym)*(c1*np.exp(lam1*X/xm)+c2*np.exp(lam2*X/xm)+1)*np.cos(np.pi*Y/ym)
v=(A/xm)*(c1*lam1*np.exp(lam1*X/xm)+c2*lam2*np.exp(lam2*X/xm))*np.sin(np.pi*Y/ym)
    
#Velocity for the boundary for advection (the reason for second upwind differencing method
xr=np.linspace(dx,xm,nx,endpoint=True)
xl=np.linspace(0,xm-dx,nx,endpoint=True)
ya=np.linspace(-ym+dy,ym,2*ny,endpoint=True)
yb=np.linspace(-ym,ym-dy,2*ny,endpoint=True)
    
Xo,Yo=np.meshgrid(xr,y)
ur=-(np.pi*A/ym)*(c1*np.exp(lam1*Xo/xm)+c2*np.exp(lam2*Xo/xm)+1)*np.cos(np.pi*Yo/ym)
Xo,Yo=np.meshgrid(xl,y)
ul=-(np.pi*A/ym)*(c1*np.exp(lam1*Xo/xm)+c2*np.exp(lam2*Xo/xm)+1)*np.cos(np.pi*Yo/ym)
Xo,Yo=np.meshgrid(x,ya)
va=(A/xm)*(c1*lam1*np.exp(lam1*Xo/xm)+c2*lam2*np.exp(lam2*Xo/xm))*np.sin(np.pi*Yo/ym)
Xo,Yo=np.meshgrid(x,yb)
vb=(A/xm)*(c1*lam1*np.exp(lam1*Xo/xm)+c2*lam2*np.exp(lam2*Xo/xm))*np.sin(np.pi*Yo/ym)
Kx=1000*np.ones(np.shape(X))
tym=int(ny-30); typ=int(ny+31)
tx1=int(np.ceil(nx*45/130)); tx2=int(nx-np.ceil(nx*50/130));
Kx[tym-1:typ,0:tx1]=6300
Kx[tym-1:typ,tx1-1:tx2]=4000
Kx[tym-1:typ,tx2-1:nx+1]=5000
Ky=Kx
    
dt=1/(np.abs(u/dx)+np.abs(v/dy)+2*Kx/(dx**2)+2*Ky/(dy**2));
dt=0.98*dt.min()
nyrs=1
tyear=365.25*24*60*60
nt=int(np.ceil(nyrs*tyear/dt))
ntyr=(tyear/dt)
        
Dx=Kx*dt/(dx**2)
Dy=Ky*dt/(dy**2)
    
wxm=Dx+0 #np.copy +0 *1  equal in matrix is like ponter in C/C++
wxm[:,0]=0
wxp=Dx
wxp[:,nx-1]=0
wym=Dy+0
wym[0,:]=0
wyp=Dy
wyp[2*ny-1,:]=0

wxp[ur<0]=wxp[ur<0]-dt*ur[ur<0]/dx
wxm[ul>0]=wxm[ul>0]=dt*ul[ul>0]/dx
wyp[va<0]=wyp[va<0]-dt*va[va<0]/dy
wym[vb>0]=wym[vb>0]+dt*vb[vb>0]/dy
w0=1-wxp-wxm-wyp-wym

C=np.zeros((2*ny,nx)); 
C33=np.zeros((2*ny,nx)); 
C34=np.zeros((2*ny,nx)); 
C35=np.zeros((2*ny,nx)); 
C36=np.zeros((2*ny,nx)); 
Cnew=C+0; 
Cnew33=C33+0;
Cnew34=C34+0;
Cnew35=C35+0;
Cnew36=C36+0;

alpha18=0.982
Csto=np.zeros((2*ny,nx,20))
C33sto=np.zeros((2*ny,nx,20))
C34sto=np.zeros((2*ny,nx,20))
C35sto=np.zeros((2*ny,nx,20))
C36sto=np.zeros((2*ny,nx,20))

Cbound=241 #umol/kg
nw=100
J=1.0/tyear*dt
C=Cbound*np.ones(np.shape(C))
C33=Cbound*1.000372024695310*np.ones(np.shape(C33))
C34=Cbound*1.000699*np.ones(np.shape(C34))
C35=Cbound*1.002052334599618*np.ones(np.shape(C35))
C36=Cbound*1.003345207262840*np.ones(np.shape(C36))

theta33=0.5204;
theta35=1.538;
theta36=2.049; 

J34=J*alpha18;
J33=J*alpha18**theta33;
J36=J*alpha18**theta36;
J35=J*alpha18**theta35;

@jit(parallel=True)
def Advdiff(nx,ny,nt,nw,Cbound,Cnew,Cnew33,Cnew34,Cnew35,Cnew36,C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto,ntyr,w0,wxm,wxp,wym,wyp,J,J33,J34,J35,J36):
#tic=time.time()

    ry=2*ny-1
    rx=nx-1
    nyr=1
    

    for i in range(nt):
    
        C[0:nw,0:nx]=Cbound
        C33[0:nw,0:nx]=Cbound*1.000372024695310
        C34[0:nw,0:nx]=Cbound*1.000699
        C35[0:nw,0:nx]=Cbound*1.002052334599618
        C36[0:nw,0:nx]=Cbound*1.003345207262840
    
        #regular iteration for 32
        Cnew[1:ry,1:rx]=w0[1:ry,1:rx]*C[1:ry,1:rx]+wxm[1:ry,1:rx]*C[1:ry,0:rx-1]+wxp[1:ry,1:rx]*C[1:ry,2:rx+1]+wym[1:ry,1:rx]*C[0:ry-1,1:rx]+wyp[1:ry,1:rx]*C[2:ry+1,1:rx]
        #for each wall
        Cnew[1:ry,0]=w0[1:ry,0]*C[1:ry,0]+wxp[1:ry,0]*C[1:ry,1]+wym[1:ry,0]*C[0:ry-1,0]+wyp[1:ry,0]*C[2:ry+1,0]
        Cnew[1:ry,rx]=w0[1:ry,rx]*C[1:ry,rx]+wxm[1:ry,rx]*C[1:ry,rx-1]+wym[1:ry,rx]*C[0:ry-1,rx]+wyp[1:ry,rx]*C[2:ry+1,rx]
        Cnew[0,1:rx]=w0[0,1:rx]*C[0,1:rx]+wxm[0,1:rx]*C[0,0:rx-1]+wxp[0,1:rx]*C[0,2:rx+1]+wyp[0,1:rx]*C[1,1:rx]
        Cnew[ry,1:rx]=w0[ry,1:rx]*C[ry,1:rx]+wxm[ry,1:rx]*C[ry,0:rx-1]+wxp[ry,1:rx]*C[ry,2:rx+1]+wym[ry,1:rx]*C[ry-1,1:rx]
        #for corner
        Cnew[0,0]=w0[0,0]*C[0,0]+wxp[0,0]*C[0,1]+wyp[0,0]*C[1,0]
        Cnew[0,rx]=w0[0,rx]*C[0,rx]+wxm[0,rx]*C[0,rx-1]+wyp[0,rx]*C[1,rx]
        Cnew[ry,0]=w0[ry,0]*C[ry,0]+wxp[ry,0]*C[ry,1]+wym[ry,0]*C[ry-1,0]
        Cnew[ry,rx]=w0[ry,rx]*C[ry,rx]+wxm[ry,rx]*C[ry,rx-1]+wym[ry,rx]*C[ry-1,rx]
    
        # for 33
        Cnew33[1:ry,1:rx]=w0[1:ry,1:rx]*C33[1:ry,1:rx]+wxm[1:ry,1:rx]*C33[1:ry,0:rx-1]+wxp[1:ry,1:rx]*C33[1:ry,2:rx+1]+wym[1:ry,1:rx]*C33[0:ry-1,1:rx]+wyp[1:ry,1:rx]*C33[2:ry+1,1:rx]
        #for each wall
        Cnew33[1:ry,0]=w0[1:ry,0]*C33[1:ry,0]+wxp[1:ry,0]*C33[1:ry,1]+wym[1:ry,0]*C33[0:ry-1,0]+wyp[1:ry,0]*C33[2:ry+1,0]
        Cnew33[1:ry,rx]=w0[1:ry,rx]*C33[1:ry,rx]+wxm[1:ry,rx]*C33[1:ry,rx-1]+wym[1:ry,rx]*C33[0:ry-1,rx]+wyp[1:ry,rx]*C33[2:ry+1,rx]
        Cnew33[0,1:rx]=w0[0,1:rx]*C33[0,1:rx]+wxm[0,1:rx]*C33[0,0:rx-1]+wxp[0,1:rx]*C33[0,2:rx+1]+wyp[0,1:rx]*C33[1,1:rx]
        Cnew33[ry,1:rx]=w0[ry,1:rx]*C33[ry,1:rx]+wxm[ry,1:rx]*C33[ry,0:rx-1]+wxp[ry,1:rx]*C33[ry,2:rx+1]+wym[ry,1:rx]*C33[ry-1,1:rx]
        #for corner
        Cnew33[0,0]=w0[0,0]*C33[0,0]+wxp[0,0]*C33[0,1]+wyp[0,0]*C33[1,0]
        Cnew33[0,rx]=w0[0,rx]*C33[0,rx]+wxm[0,rx]*C33[0,rx-1]+wyp[0,rx]*C33[1,rx]
        Cnew33[ry,0]=w0[ry,0]*C33[ry,0]+wxp[ry,0]*C33[ry,1]+wym[ry,0]*C33[ry-1,0]
        Cnew33[ry,rx]=w0[ry,rx]*C33[ry,rx]+wxm[ry,rx]*C33[ry,rx-1]+wym[ry,rx]*C33[ry-1,rx]
        
        #regular iteration for 34
        Cnew34[1:ry,1:rx]=w0[1:ry,1:rx]*C34[1:ry,1:rx]+wxm[1:ry,1:rx]*C34[1:ry,0:rx-1]+wxp[1:ry,1:rx]*C34[1:ry,2:rx+1]+wym[1:ry,1:rx]*C34[0:ry-1,1:rx]+wyp[1:ry,1:rx]*C34[2:ry+1,1:rx]
        #for each wall
        Cnew34[1:ry,0]=w0[1:ry,00]*C34[1:ry,0]+wxp[1:ry,0]*C34[1:ry,1]+wym[1:ry,0]*C34[0:ry-1,0]+wyp[1:ry,0]*C34[2:ry+1,0]
        Cnew34[1:ry,rx]=w0[1:ry,rx]*C34[1:ry,rx]+wxm[1:ry,rx]*C34[1:ry,rx-1]+wym[1:ry,rx]*C34[0:ry-1,rx]+wyp[1:ry,rx]*C34[2:ry+1,rx]
        Cnew34[0,1:rx]=w0[0,1:rx]*C34[0,1:rx]+wxm[0,1:rx]*C34[0,0:rx-1]+wxp[0,1:rx]*C34[0,2:rx+1]+wyp[0,1:rx]*C34[1,1:rx]
        Cnew34[ry,1:rx]=w0[ry,1:rx]*C34[ry,1:rx]+wxm[ry,1:rx]*C34[ry,0:rx-1]+wxp[ry,1:rx]*C34[ry,2:rx+1]+wym[ry,1:rx]*C34[ry-1,1:rx]
        #for corner
        Cnew34[0,0]=w0[0,0]*C34[0,0]+wxp[0,0]*C34[0,1]+wyp[0,0]*C34[1,0]
        Cnew34[0,rx]=w0[0,rx]*C34[0,rx]+wxm[0,rx]*C34[0,rx-1]+wyp[0,rx]*C34[1,rx]
        Cnew34[ry,0]=w0[ry,0]*C34[ry,0]+wxp[ry,0]*C34[ry,1]+wym[ry,0]*C34[ry-1,0]
        Cnew34[ry,rx]=w0[ry,rx]*C34[ry,rx]+wxm[ry,rx]*C34[ry,rx-1]+wym[ry,rx]*C34[ry-1,rx]
    
        #regular iteration for 35
        Cnew35[1:ry,1:rx]=w0[1:ry,1:rx]*C35[1:ry,1:rx]+wxm[1:ry,1:rx]*C35[1:ry,0:rx-1]+wxp[1:ry,1:rx]*C35[1:ry,2:rx+1]+wym[1:ry,1:rx]*C35[0:ry-1,1:rx]+wyp[1:ry,1:rx]*C35[2:ry+1,1:rx]
        #for each wall
        Cnew35[1:ry,0]=w0[1:ry,0]*C35[1:ry,0]+wxp[1:ry,0]*C35[1:ry,1]+wym[1:ry,0]*C35[0:ry-1,0]+wyp[1:ry,0]*C35[2:ry+1,0]
        Cnew35[1:ry,rx]=w0[1:ry,rx]*C35[1:ry,rx]+wxm[1:ry,rx]*C35[1:ry,rx-1]+wym[1:ry,rx]*C35[0:ry-1,rx]+wyp[1:ry,rx]*C35[2:ry+1,rx]
        Cnew35[0,1:rx]=w0[0,1:rx]*C35[0,1:rx]+wxm[0,1:rx]*C35[0,0:rx-1]+wxp[0,1:rx]*C35[0,2:rx+1]+wyp[0,1:rx]*C35[1,1:rx]
        Cnew35[ry,1:rx]=w0[ry,1:rx]*C35[ry,1:rx]+wxm[ry,1:rx]*C35[ry,0:rx-1]+wxp[ry,1:rx]*C35[ry,2:rx+1]+wym[ry,1:rx]*C35[ry-1,1:rx]
        #for corner
        Cnew35[0,0]=w0[0,0]*C35[0,0]+wxp[0,0]*C35[0,1]+wyp[0,0]*C35[1,0]
        Cnew35[0,rx]=w0[0,rx]*C35[0,rx]+wxm[0,rx]*C35[0,rx-1]+wyp[0,rx]*C35[1,rx]
        Cnew35[ry,0]=w0[ry,0]*C35[ry,0]+wxp[ry,0]*C35[ry,1]+wym[ry,0]*C35[ry-1,0]
        Cnew35[ry,rx]=w0[ry,rx]*C35[ry,rx]+wxm[ry,rx]*C35[ry,rx-1]+wym[ry,rx]*C35[ry-1,rx]
    
        #regular iteration for 36
        Cnew36[1:ry,1:rx]=w0[1:ry,1:rx]*C[1:ry,1:rx]+wxm[1:ry,1:rx]*C[1:ry,0:rx-1]+wxp[1:ry,1:rx]*C[1:ry,2:rx+1]+wym[1:ry,1:rx]*C[0:ry-1,1:rx]+wyp[1:ry,1:rx]*C[2:ry+1,1:rx]
        #for each wall
        Cnew36[1:ry,0]=w0[1:ry,0]*C36[1:ry,0]+wxp[1:ry,0]*C36[1:ry,1]+wym[1:ry,0]*C36[0:ry-1,0]+wyp[1:ry,0]*C36[2:ry+1,0]
        Cnew36[1:ry,rx]=w0[1:ry,rx]*C36[1:ry,rx]+wxm[1:ry,rx]*C36[1:ry,rx-1]+wym[1:ry,rx]*C36[0:ry-1,rx]+wyp[1:ry,rx]*C36[2:ry+1,rx]
        Cnew36[0,1:rx]=w0[0,1:rx]*C36[0,1:rx]+wxm[0,1:rx]*C36[0,0:rx-1]+wxp[0,1:rx]*C36[0,2:rx+1]+wyp[0,1:rx]*C36[1,1:rx]
        Cnew36[ry,1:rx]=w0[ry,1:rx]*C36[ry,1:rx]+wxm[ry,1:rx]*C36[ry,0:rx-1]+wxp[ry,1:rx]*C36[ry,2:rx+1]+wym[ry,1:rx]*C36[ry-1,1:rx]
        #for corner
        Cnew36[0,0]=w0[0,0]*C36[0,0]+wxp[0,0]*C36[0,1]+wyp[0,0]*C36[1,0]
        Cnew36[0,rx]=w0[0,rx]*C36[0,rx]+wxm[0,rx]*C36[0,rx-1]+wyp[0,rx]*C36[1,rx]
        Cnew36[ry,0]=w0[ry,0]*C36[ry,0]+wxp[ry,0]*C36[ry,1]+wym[ry,0]*C36[ry-1,0]
        Cnew36[ry,rx]=w0[ry,rx]*C36[ry,rx]+wxm[ry,rx]*C36[ry,rx-1]+wym[ry,rx]*C36[ry-1,rx]
    
        Cnew33[Cnew>0]=Cnew33[Cnew>0]-J33*Cnew33[Cnew>0]/Cnew[Cnew>0]
        Cnew34[Cnew>0]=Cnew34[Cnew>0]-J34*Cnew34[Cnew>0]/Cnew[Cnew>0]
        Cnew35[Cnew>0]=Cnew35[Cnew>0]-J35*Cnew35[Cnew>0]/Cnew[Cnew>0]
        Cnew36[Cnew>0]=Cnew36[Cnew>0]-J36*Cnew36[Cnew>0]/Cnew[Cnew>0]
        
        Cnew=Cnew-J
        Cnew[Cnew<0]=0
        Cnew34[Cnew==0]=0
        Cnew33[Cnew==0]=0
        Cnew36[Cnew==0]=0
        Cnew35[Cnew==0]=0
        Cnew34[Cnew34<0]=0
        Cnew33[Cnew33<0]=0
        Cnew36[Cnew36<0]=0
        Cnew35[Cnew35<0]=0
    
        C=Cnew;
        C34=Cnew34;
        C33=Cnew33;
        C36=Cnew36;
        C35=Cnew35;
        #Sums(i,1)=sum(sum(C));
   
        if np.floor(i/ntyr/50)>=nyr:
            nyr=np.floor(i/ntyr/50)
            Csto[:,:,nyr]=C; 
            C34sto[:,:,nyr]=C34
            C33sto[:,:,nyr]=C33
            C35sto[:,:,nyr]=C35 
            C36sto[:,:,nyr]=C36


    #toc=time.time()
    #print(toc-tic)      
    return C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto
tic=time.time()
[C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto]=Advdiff(nx,ny,1,nw,Cbound,Cnew,Cnew33,Cnew34,Cnew35,Cnew36,C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto,ntyr,w0,wxm,wxp,wym,wyp,J,J33,J34,J35,J36)
toc=time.time()
print(toc-tic)
tic=time.time()
[C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto]=Advdiff(nx,ny,1,nw,Cbound,Cnew,Cnew33,Cnew34,Cnew35,Cnew36,C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto,ntyr,w0,wxm,wxp,wym,wyp,J,J33,J34,J35,J36)
toc=time.time()
print(toc-tic)
