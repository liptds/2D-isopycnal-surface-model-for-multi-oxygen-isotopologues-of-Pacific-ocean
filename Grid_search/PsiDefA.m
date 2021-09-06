xm=1.05e7; ym=4.5e6; %domain dimensions
eps=0.015; %boundary current 
nx=ceil(10/eps); ny=ceil(nx*ym/xm); %rectangular grid
dx=xm/nx;dy=ym/ny; 
x=0:dx:xm;y=-ym:dy:ym;
[X,Y]=meshgrid(x,y);
lam1=(-1+(1+4*pi^2*eps^2)^0.5)/2/eps;
lam2=(-1-(1+4*pi^2*eps^2)^0.5)/2/eps;
c1=(1-exp(lam2))/(exp(lam2)-exp(lam1));
c2=-(1+c1);