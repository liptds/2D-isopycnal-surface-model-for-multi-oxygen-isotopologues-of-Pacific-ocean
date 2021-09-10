
%initial concentration and isotopic fractionation
C=zeros(2*ny,nx); 
C33=zeros(2*ny,nx); 
C34=zeros(2*ny,nx); 
C35=zeros(2*ny,nx); 
C36=zeros(2*ny,nx); 
Cnew=C; 
% Cnew33=C33;
% Cnew34=C34;
% Cnew35=C35;
% Cnew36=C36;



alpha18=0.982;%0.982,0.990
%storage every 50 years
% Csto=zeros(2*ny,nx,20); 
% C33sto=zeros(2*ny,nx,20); 
% C34sto=zeros(2*ny,nx,20); 
% C35sto=zeros(2*ny,nx,20); 
% C36sto=zeros(2*ny,nx,20); 


%test for sharp boundary
%C=zeros(2*ny,nx);
%C(:,334:667)=1;
%Sums=zeros(51,1);

Cbound=273; %umol/kg
nw=50;
nw2=454;
% np=6;
J=2.4/tyear*dt;%umol/kg/dt
% C=Cbound*ones(2*ny,nx);
% C33=Cbound*1.000372024695310*ones(2*ny,nx);
% C34=Cbound*1.000699*ones(2*ny,nx);
% C35=Cbound*1.002052334599618*ones(2*ny,nx);
% C36=Cbound*1.003345207262840*ones(2*ny,nx);

% fraction=0.5;
% R18eq=1.000699;
% R17eq=1.000372024695310;
% R35eq=1.002052334599618;
% R36eq=1.003345207262840;
R18eq=1+(-0.730+427/(273.15+14))/1000;
R17eq=exp(8*10^-6+0.518*log(R18eq));
R35eq=R18eq*R17eq*(0.98/1000+1);
R36eq=R18eq^2*(1.944/1000+1);
R18p = 1 - 20.172/1000;
R17p = 1 - 10.275/1000;
% R18p=1-20.014/1000;
% R17p=1-10.222/1000;
R35p=R17p*R18p*(1-0.2/1000);
R36p=R18p^2*(1-0.4/1000);
% R18mix=R18eq*fraction+R18p*(1-fraction);
% R17mix=R17eq*fraction+R17p*(1-fraction);
% R36mix=R36eq*fraction+R36p*(1-fraction);
% R35mix=R35eq*fraction+R35p*(1-fraction);
% D17Omix=1e6*(log(R17mix)-0.518*log(R18mix));

% R18mix=1-1.37/1000;
% R17mix=exp(0.518*log(R18mix)+155.407*10^(-6));
% R36mix=R18mix^2*(1+1.267/1000);
% R35mix=R17mix*R18mix*(1+0.709/1000);

theta33=0.516;%0.5204,0.5234
theta35=1.534;%1.538,1.555
theta36=2.048;%2.049,2.061
% R34=0.004110354;
% R33=0.000769476;
% R35=1.58298E-06;
% R36=4.23077E-06;
% L32 = 1/(1+R34+R33+R35+R36);
% H33 = R33/(1+R34+R33+R35+R36);
% H34 = R34/(1+R34+R33+R35+R36);
% H35 = R35/(1+R34+R33+R35+R36);
% H36 = R36/(1+R34+R33+R35+R36);
% R17 = (R33 + 2*(R33/2)^2 + R35)/(2 + R33 + (R34 -(R33/2)^2)); 
% R18 = (R34 - (R33/2)^2 + R35 + 2*R36)/(2 + R33 + (R34 -(R33/2)^2));

J34=J*alpha18;
J33=J*alpha18^theta33;
J36=J*alpha18^theta36;
J35=J*alpha18^theta35;
