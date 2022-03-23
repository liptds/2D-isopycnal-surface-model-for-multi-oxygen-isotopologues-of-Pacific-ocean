
%initial concentration and isotopic fractionation factors
C=zeros(2*ny,nx); 
C33=zeros(2*ny,nx); 
C34=zeros(2*ny,nx); 
C35=zeros(2*ny,nx); 
C36=zeros(2*ny,nx); 
Cnew=C; 
Cnew33=C33;
Cnew34=C34;
Cnew35=C35;
Cnew36=C36;

alpha18=0.982;

Cbound=317; % O2 concentration when ventilated, umol/kg (273 for 25.8-26.2 
            % surface and 317 for 26.5-26.9 surface)
nw=40; % number of boxes in Y-dimension of Southern ventilation region
J1=1.6/tyear*dt; % O2 consumption rate in gyre, umol/kg/dt
J2 = 20/tyear*dt; % Additional O2 consumption rate in equatorial region, umol/kg/dt

% Ventilated (solubility equilibrium) endmember
R18eq=1+(-0.730+427/(273.15+14))/1000;
R17eq=exp(8*10^-6+0.518*log(R18eq));
R35eq=R18eq*R17eq*(0.98/1000+1);
R36eq=R18eq^2*(1.944/1000+1);

% Implicit photosynthesis endmember
% R18eq=1+(-0.730+427/(273.15+14)-0.4)/1000;
% R17eq=exp(30*10^-6+0.518*log(R18eq));
% R35eq=R18eq*R17eq*(0.98/1000+1);
% R36eq=R18eq^2*(1.90/1000+1);

% Mass-dependent fractionation exponents
theta33=0.520;
theta35=1.534;
theta36=2.048;

% Isotopologue-specific respiration rates
J34=J1*alpha18;
J33=J1*alpha18^theta33;
J36=J1*alpha18^theta36;
J35=J1*alpha18^theta35;

J342=J2*alpha18;
J332=J2*alpha18^theta33;
J362=J2*alpha18^theta36;
J352=J2*alpha18^theta35;