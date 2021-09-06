# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 17:31:53 2019

@author: Boda Li
"""

import numpy as np
#from numba import jit

#@jit
#import matplotlib.pyplot as plt #draw contour plot
def Psi():
#tic=time.time()
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
#toc=time.time()
#print(toc-tic)   
    return A,dx,dy,xm,ym,nx,ny,c1,c2,lam1,lam2
#contours=plt.contour(X,Y,Psi,10)
#plt.clabel(contours,inline=True)
#plt.savefig('Fig1.pdf')