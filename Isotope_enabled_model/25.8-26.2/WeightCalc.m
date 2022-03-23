%finite time von Neuman stability analysis
Kx=1200*ones(size(X));
ty=[ny-30:ny+31];
ty2 = [ny-62:ny+63]; 
tx1=[1:ceil(nx*45/130)];
tx2=[ceil(nx*45/130):nx-ceil(nx*50/130)];
tx3=[nx-ceil(nx*50/130):nx];
Kx(ty,tx1)=6300*0.8; % Surface eddy diffusivity x scaling factor from Cole et al. (2015)
Kx(ty,tx2)=4000*0.8; % Surface eddy diffusivity x scaling factor from Cole et al. (2015)
Kx(ty,tx3)=5000*0.8; % Surface eddy diffusivity x scaling factor from Cole et al. (2015)
Ky=Kx; %isotropic eddy diffusivity
dt=min(min(1./(abs(u/dx)+abs(v/dy)+2*Kx/dx.^2+2*Ky/dy.^2)));
dt=0.98*dt; %safety reuduction
nyrs=1000;
tyear=365.25*24*60*60;
nt=ceil(nyrs*tyear/dt); 
ntyr=(tyear/dt);

%Advection diffusion calculation
Dx=Kx*dt/(dx^2);
Dy=Ky*dt/(dy^2);

wxm=Dx;
wxm(:,1)=0; % i-1>0 deal with the bottom boundary
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


