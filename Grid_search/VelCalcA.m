x=[0.5*dx:dx:xm-0.5*dx];
y=[-ym+0.5*dy:dy:ym-0.5*dy];
[X,Y]=meshgrid(x,y);
u=-(pi*A/ym)*(c1*exp(lam1*X/xm)+c2*exp(lam2*X/xm)+1).*cos(pi*Y/ym);
v=(A/xm)*(c1*lam1*exp(lam1*X/xm)+c2*lam2*exp(lam2*X/xm)).*sin(pi*Y/ym);

%Velocity for the boundary for advection (the reason for second upwind
%differencing method
xr=[dx:dx:xm];xl=[0:dx:xm-dx];
ya=[-ym+dy:dy:ym];yb=[-ym:dy:ym-dy];
[Xo,Yo]=meshgrid(xr,y);
ur=-(pi*A/ym)*(c1*exp(lam1*Xo/xm)+c2*exp(lam2*Xo/xm)+1).*cos(pi*Yo/ym);
[Xo,Yo]=meshgrid(xl,y);
ul=-(pi*A/ym)*(c1*exp(lam1*Xo/xm)+c2*exp(lam2*Xo/xm)+1).*cos(pi*Yo/ym);
[Xo,Yo]=meshgrid(x,ya);
va=(A/xm)*(c1*lam1*exp(lam1*Xo/xm)+c2*lam2*exp(lam2*Xo/xm)).*sin(pi*Yo/ym);
[Xo,Yo]=meshgrid(x,yb);
vb=(A/xm)*(c1*lam1*exp(lam1*Xo/xm)+c2*lam2*exp(lam2*Xo/xm)).*sin(pi*Yo/ym);