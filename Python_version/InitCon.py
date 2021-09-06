# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:36:26 2019

@author: Boda Li
"""

import numpy as np
#from numba import jit

#@jit
#tic=time.time()
def Initcon(nx,ny,tyear,dt):
#initial concentration and isotopic fractionation
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
#R34=0.004110354;
#R33=0.000769476;
#R35=1.58298E-06;
#R36=4.23077E-06;
#L32 = 1/(1+R34+R33+R35+R36);
#H33 = R33/(1+R34+R33+R35+R36);
#H34 = R34/(1+R34+R33+R35+R36);
#H35 = R35/(1+R34+R33+R35+R36);
#H36 = R36/(1+R34+R33+R35+R36);
#R17 = (R33 + 2*(R33/2)^2 + R35)/(2 + R33 + (R34 -(R33/2)^2)); 
#R18 = (R34 - (R33/2)^2 + R35 + 2*R36)/(2 + R33 + (R34 -(R33/2)^2));

    J34=J*alpha18;
    J33=J*alpha18**theta33;
    J36=J*alpha18**theta36;
    J35=J*alpha18**theta35;
#toc=time.time()
#print(toc-tic)    
    return Cbound,C,C33,C34,C35,C36,Cnew,Cnew33,Cnew34,Cnew35,Cnew36,Csto,C33sto,C34sto,C35sto,C36sto,J,J33,J34,J35,J36,nw