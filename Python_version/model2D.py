# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 21:23:05 2019

@author: Boda Li
"""

from PsiDef import Psi
from VelCalc import Vel
from WeightCalc import Weight
from InitCon import Initcon
import time
from AdvDiff import Advdiff
tic=time.time()
[A,dx,dy,xm,ym,nx,ny,c1,c2,lam1,lam2]=Psi()
[X,Y,u,v,ur,ul,va,vb]=Vel(A,dx,dy,xm,ym,nx,ny,c1,c2,lam1,lam2)
[w0,wxp,wxm,wyp,wym,nt,tyear,dt,ntyr]=Weight(X,nx,ny,u,v,dx,dy,ur,ul,va,vb)
[Cbound,C,C33,C34,C35,C36,Cnew,Cnew33,Cnew34,Cnew35,Cnew36,Csto,C33sto,C34sto,C35sto,C36sto,J,J33,J34,J35,J36,nw]=Initcon(nx,ny,tyear,dt)
[C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto]=Advdiff(nx,ny,nt,nw,Cbound,Cnew,Cnew33,Cnew34,Cnew35,Cnew36,C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto,ntyr,w0,wxm,wxp,wym,wyp,J,J33,J34,J35,J36)

toc=time.time()
print(toc-tic)


