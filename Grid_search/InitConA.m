
%initial concentration and isotopic fractionation
C=zeros(2*ny,nx); 

Cnew=C; 


Cbound=273; %umol/kg
nw=50;
nw2=454;
% np=6;
J=J_res/tyear*dt;%umol/kg/dt