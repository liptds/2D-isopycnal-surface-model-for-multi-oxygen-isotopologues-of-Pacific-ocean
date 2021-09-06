# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:33:19 2019

@author: Boda Li
"""
import numpy as np
#import matplotlib.pyplot as plt
#from numba import jit

#@jit

def Vel(A,dx,dy,xm,ym,nx,ny,c1,c2,lam1,lam2):
#tic=time.time()
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
#toc=time.time()
#print(toc-tic) 
    return X,Y,u,v,ur,ul,va,vb

    #contours=plt.contour(X,Y,u,25)
    #plt.clabel(contours,inline=True)
    #plt.savefig('Fig2.pdf')
    
    #contours=plt.contour(X,Y,v,25)
    #plt.clabel(contours,inline=True)
    #plt.savefig('Fig3.pdf')