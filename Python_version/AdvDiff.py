# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 16:28:35 2019

@author: Boda Li
"""



import numpy as np

#from numba import jit

#@jit
def Advdiff(nx,ny,nt,nw,Cbound,Cnew,Cnew33,Cnew34,Cnew35,Cnew36,C,C33,C34,C35,C36,Csto,C33sto,C34sto,C35sto,C36sto,ntyr,w0,wxm,wxp,wym,wyp,J,J33,J34,J35,J36):
#tic=time.time()

    ry=2*ny-1
    rx=nx-1
    nyr=1
    

    for i in range(0,nt):
    
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