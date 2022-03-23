% Two-dimensional isopycnal advection-diffusion-reaction model for the 
% simulation of oxygen isotopologues in the Pacific ocean. Runs five
% programs in sequence: PsiDef, VelCalc, WeightCalc, InitCon, and
% AdvDiffGPU. Program utilizes hardware Graphics Processing Unit 
% acceleration via gpuArray in AdvDiffGPU.
%
% This set of programs contains the initial settings and parameters for the
% 25.8-26.2 potential density surface, including code for the "explicit 
% addition" method of simulating the effects of photosynthesis. To run in
% respiration-only mode, comment out Lines 43-47 in AdvDiffGPU.
% 
% Cite as: Li, B., H. Hu, W. M. Berelson, J. F. Adkins, and L. Y. Yeung
%   (2022). On the use of dissolved oxygen isotopologues as biogeochemical
%   tracers in the Pacific Ocean. Submitted to Journal of Geophysical
%   Research, Oceans. ESSOAr Preprint doi: 10.1002/essoar.10510808.2

clear 
PsiDef
VelCalc
WeightCalc
InitCon
AdvDiffGPU

% Convert simulation results to delta notation
d18O=NaN(size(C));
D17O=NaN(size(C));
D36=NaN(size(C));
d18O(C>0)=(C34(C>0)./C(C>0)-1)*1000;
D17O(C>0)=1e6*(log(C33(C>0)./C(C>0))-0.518*log(C34(C>0)./C(C>0)));
D36(C>0)=(C36(C>0)./C(C>0)./(C34(C>0)./C(C>0)).^2-1)*1000;
D35=(C35./C./(C34./C)./(C33./C)-1)*1000;
O2sat=C/Cbound;

% Contour plots of results
figure
[c,h]=contour(x,y,O2sat,10);
clabel(c,h)

figure
[c,h]=contour(x,y,d18O,0:1:10);
clabel(c,h)

figure
[c,h]=contour(x,y,D17O,0:10:100);
clabel(c,h)

figure
[c,h]=contour(x,y,D36, 2.5:-0.15:1.3);
clabel(c,h)