%finite time von Neuman stability analysis
Kx=D*1000*ones(size(X));
ty=[ny-30:ny+31];
tx1=[1:ceil(nx*45/130)];
tx2=[ceil(nx*45/130):nx-ceil(nx*50/130)];
tx3=[nx-ceil(nx*50/130):nx];
Kx(ty,tx1)=6300*0.8;
Kx(ty,tx2)=4000*0.8;
Kx(ty,tx3)=5000*0.8;
Ky=Kx; %isotrpic for diffusion coefficient 1000m^2/s
dt=min(min(1./(abs(u/dx)+abs(v/dy)+2*Kx/dx.^2+2*Ky/dy.^2)));
dt=0.98*dt; %safety reuduction
nyrs=500;
tyear=365.25*24*60*60;
nt=ceil(nyrs*tyear/dt); 
ntyr=(tyear/dt);

%Advection diffusion
Dx=Kx*dt/(dx^2);
Dy=Ky*dt/(dy^2);

wxm=Dx;
wxm(:,1)=0; %i-1>0 deal with the bottom boundary
wxp=Dx;
wxp(:,nx)=0;

wym=Dy;
wym(1,:)=0;
wyp=Dy;
wyp(2*ny,:)=0;

wxp(ur<0)=wxp(ur<0)-dt*ur(ur<0)/dx;
wxm(ul>0)=wxm(ul>0)+dt*ul(ul>0)/dx;
wyp(va<0)=wyp(va<0)-dt*va(va<0)/dy;
wym(vb>0)=wym(vb>0)+dt*vb(vb>0)/dy;
w0=1-wxp-wxm-wyp-wym;
