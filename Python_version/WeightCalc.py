# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 18:56:24 2019

@author: Boda Li
"""

import numpy as np
#from numba import jit

#@jit
def Weight(X,nx,ny,u,v,dx,dy,ur,ul,va,vb):
#tic=time.time()
    Kx=1000*np.ones(np.shape(X))
    tym=int(ny-30); typ=int(ny+31)
    tx1=int(np.ceil(nx*45/130)); tx2=int(nx-np.ceil(nx*50/130));
    Kx[tym-1:typ,0:tx1]=6300
    Kx[tym-1:typ,tx1-1:tx2]=4000
    Kx[tym-1:typ,tx2-1:nx+1]=5000
    Ky=Kx
    
    dt=1/(np.abs(u/dx)+np.abs(v/dy)+2*Kx/(dx**2)+2*Ky/(dy**2));
    dt=0.98*dt.min()
    nyrs=3
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
#toc=time.time()
#print(toc-tic)     
    return w0,wxp,wxm,wyp,wym,nt,tyear,dt,ntyr

