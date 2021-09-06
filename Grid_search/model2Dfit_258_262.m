clear 
tic
%for A=4*10^5:1*10^5:8*10^5
% A=2*10^5;
for A = 2*10^5:2*10^4:3*10^5
     for D=0.8:0.2:1.2
         for J_res=2.0:0.2:3.0
             clearvars -except A D J_res
             load('Initial258_262.mat')
%A=5*10^5;
%D=1;
%J_res=0.5;
        PsiDefA
        VelCalcA
        WeightCalcA
        InitConA
        AdvDiffGPUA
        save(sprintf('ResultA%.5gD%.1fJ%1.1f.mat_sigmatheta_25.8_26.2.mat',A,D,J_res),'A','D','J_res','C');

%figure
%[z1,z2]=contour(x,y,C,320:-20:0);
%clabel(z1,z2)        
         end
     end
end
toc 