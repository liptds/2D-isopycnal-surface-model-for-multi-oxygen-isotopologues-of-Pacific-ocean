ry=[2:2*ny-1]; rx=[2:nx-1];

tic
nyr=1;
nyrg=gpuArray(1);
ntyrg=gpuArray(ntyr);

Cg=gpuArray(C);
C33g=gpuArray(C33);
C34g=gpuArray(C34);
C35g=gpuArray(C35);
C36g=gpuArray(C36);
Cnewg=gpuArray(Cnew);
Cnew33g=gpuArray(C33);
Cnew34g=gpuArray(C34);
Cnew35g=gpuArray(C35);
Cnew36g=gpuArray(C36);
w0g=gpuArray(w0);
wxmg=gpuArray(wxm);
wxpg=gpuArray(wxp);
wymg=gpuArray(wym);
wypg=gpuArray(wyp);

time1=toc
tic
for i=1:nt
%   boundary conditions
    % Southern ventilation region 
    j=1:nw;
    Cg(j,:)=Cbound;
    C33g(j,:)=Cbound*R17eq;
    C34g(j,:)=Cbound*R18eq;
    C35g(j,:)=Cbound*R35eq;
    C36g(j,:)=Cbound*R36eq;

    % Northern ventilation region
    k=8:10;
    Cg(571:576,k)=Cbound;
    C33g(571:576,k)=Cbound*R17eq;
    C34g(571:576,k)=Cbound*R18eq;
    C35g(571:576,k)=Cbound*R35eq;
    C36g(571:576,k)=Cbound*R36eq;

% regular stepping iteration
    Cnewg(ry,rx)=w0g(ry,rx).*Cg(ry,rx)+wxmg(ry,rx).*Cg(ry,rx-1)+wxpg(ry,rx).*Cg(ry,rx+1)+wymg(ry,rx).*Cg(ry-1,rx)+wypg(ry,rx).*Cg(ry+1,rx);
    %for each wall
    Cnewg(ry,1)=w0g(ry,1).*Cg(ry,1)+wxpg(ry,1).*Cg(ry,2)+wymg(ry,1).*Cg(ry-1,1)+wypg(ry,1).*Cg(ry+1,1);
    Cnewg(ry,nx)=w0g(ry,nx).*Cg(ry,nx)+wxmg(ry,nx).*Cg(ry,nx-1)+wymg(ry,nx).*Cg(ry-1,nx)+wypg(ry,nx).*Cg(ry+1,nx);
    Cnewg(1,rx)=w0g(1,rx).*Cg(1,rx)+wxmg(1,rx).*Cg(1,rx-1)+wxpg(1,rx).*Cg(1,rx+1)+wypg(1,rx).*Cg(2,rx);
    Cnewg(2*ny,rx)=w0g(2*ny,rx).*Cg(2*ny,rx)+wxmg(2*ny,rx).*Cg(2*ny,rx-1)+wxpg(2*ny,rx).*Cg(2*ny,rx+1)+wymg(2*ny,rx).*Cg(2*ny-1,rx);
    %for corner
    Cnewg(1,1)=w0g(1,1).*Cg(1,1)+wxpg(1,1).*Cg(1,2)+wypg(1,1).*Cg(2,1);
    Cnewg(1,nx)=w0g(1,nx).*Cg(1,nx)+wxmg(1,nx).*Cg(1,nx-1)+wypg(1,nx).*Cg(2,nx);
    Cnewg(2*ny,1)=w0g(2*ny,1).*Cg(2*ny,1)+wxpg(2*ny,1).*Cg(2*ny,2)+wymg(2*ny,1).*Cg(2*ny-1,1);
    Cnewg(2*ny,nx)=w0g(2*ny,nx).*Cg(2*ny,nx)+wxmg(2*ny,nx).*Cg(2*ny,nx-1)+wymg(2*ny,nx).*Cg(2*ny-1,nx);
    
    %for other isotoplogues
    Cnew34g(ry,rx)=w0g(ry,rx).*C34g(ry,rx)+wxmg(ry,rx).*C34g(ry,rx-1)+wxpg(ry,rx).*C34g(ry,rx+1)+wymg(ry,rx).*C34g(ry-1,rx)+wypg(ry,rx).*C34g(ry+1,rx);
    %for each wall
    Cnew34g(ry,1)=w0g(ry,1).*C34g(ry,1)+wxpg(ry,1).*C34g(ry,2)+wymg(ry,1).*C34g(ry-1,1)+wypg(ry,1).*C34g(ry+1,1);
    Cnew34g(ry,nx)=w0g(ry,nx).*C34g(ry,nx)+wxmg(ry,nx).*C34g(ry,nx-1)+wymg(ry,nx).*C34g(ry-1,nx)+wypg(ry,nx).*C34g(ry+1,nx);
    Cnew34g(1,rx)=w0g(1,rx).*C34g(1,rx)+wxmg(1,rx).*C34g(1,rx-1)+wxpg(1,rx).*C34g(1,rx+1)+wypg(1,rx).*C34g(2,rx);
    Cnew34g(2*ny,rx)=w0g(2*ny,rx).*C34g(2*ny,rx)+wxmg(2*ny,rx).*C34g(2*ny,rx-1)+wxpg(2*ny,rx).*C34g(2*ny,rx+1)+wymg(2*ny,rx).*C34g(2*ny-1,rx);
    %for corner
    Cnew34g(1,1)=w0g(1,1).*C34g(1,1)+wxpg(1,1).*C34g(1,2)+wypg(1,1).*C34g(2,1);
    Cnew34g(1,nx)=w0g(1,nx).*C34g(1,nx)+wxmg(1,nx).*C34g(1,nx-1)+wypg(1,nx).*C34g(2,nx);
    Cnew34g(2*ny,1)=w0g(2*ny,1).*C34g(2*ny,1)+wxpg(2*ny,1).*C34g(2*ny,2)+wymg(2*ny,1).*C34g(2*ny-1,1);
    Cnew34g(2*ny,nx)=w0g(2*ny,nx).*C34g(2*ny,nx)+wxmg(2*ny,nx).*C34g(2*ny,nx-1)+wymg(2*ny,nx).*C34g(2*ny-1,nx);
   
    Cnew33g(ry,rx)=w0g(ry,rx).*C33g(ry,rx)+wxmg(ry,rx).*C33g(ry,rx-1)+wxpg(ry,rx).*C33g(ry,rx+1)+wymg(ry,rx).*C33g(ry-1,rx)+wypg(ry,rx).*C33g(ry+1,rx);
    %for each wall
    Cnew33g(ry,1)=w0g(ry,1).*C33g(ry,1)+wxpg(ry,1).*C33g(ry,2)+wymg(ry,1).*C33g(ry-1,1)+wypg(ry,1).*C33g(ry+1,1);
    Cnew33g(ry,nx)=w0g(ry,nx).*C33g(ry,nx)+wxmg(ry,nx).*C33g(ry,nx-1)+wymg(ry,nx).*C33g(ry-1,nx)+wypg(ry,nx).*C33g(ry+1,nx);
    Cnew33g(1,rx)=w0g(1,rx).*C33g(1,rx)+wxmg(1,rx).*C33g(1,rx-1)+wxpg(1,rx).*C33g(1,rx+1)+wypg(1,rx).*C33g(2,rx);
    Cnew33g(2*ny,rx)=w0g(2*ny,rx).*C33g(2*ny,rx)+wxmg(2*ny,rx).*C33g(2*ny,rx-1)+wxpg(2*ny,rx).*C33g(2*ny,rx+1)+wymg(2*ny,rx).*C33g(2*ny-1,rx);
    %for corner
    Cnew33g(1,1)=w0g(1,1).*C33g(1,1)+wxpg(1,1).*C33g(1,2)+wypg(1,1).*C33g(2,1);
    Cnew33g(1,nx)=w0g(1,nx).*C33g(1,nx)+wxmg(1,nx).*C33g(1,nx-1)+wypg(1,nx).*C33g(2,nx);
    Cnew33g(2*ny,1)=w0g(2*ny,1).*C33g(2*ny,1)+wxpg(2*ny,1).*C33g(2*ny,2)+wymg(2*ny,1).*C33g(2*ny-1,1);
    Cnew33g(2*ny,nx)=w0g(2*ny,nx).*C33g(2*ny,nx)+wxmg(2*ny,nx).*C33g(2*ny,nx-1)+wymg(2*ny,nx).*C33g(2*ny-1,nx);
    
    Cnew36g(ry,rx)=w0g(ry,rx).*C36g(ry,rx)+wxmg(ry,rx).*C36g(ry,rx-1)+wxpg(ry,rx).*C36g(ry,rx+1)+wymg(ry,rx).*C36g(ry-1,rx)+wypg(ry,rx).*C36g(ry+1,rx);
    %for each wall
    Cnew36g(ry,1)=w0g(ry,1).*C36g(ry,1)+wxpg(ry,1).*C36g(ry,2)+wymg(ry,1).*C36g(ry-1,1)+wypg(ry,1).*C36g(ry+1,1);
    Cnew36g(ry,nx)=w0g(ry,nx).*C36g(ry,nx)+wxmg(ry,nx).*C36g(ry,nx-1)+wymg(ry,nx).*C36g(ry-1,nx)+wypg(ry,nx).*C36g(ry+1,nx);
    Cnew36g(1,rx)=w0g(1,rx).*C36g(1,rx)+wxmg(1,rx).*C36g(1,rx-1)+wxpg(1,rx).*C36g(1,rx+1)+wypg(1,rx).*C36g(2,rx);
    Cnew36g(2*ny,rx)=w0g(2*ny,rx).*C36g(2*ny,rx)+wxmg(2*ny,rx).*C36g(2*ny,rx-1)+wxpg(2*ny,rx).*C36g(2*ny,rx+1)+wymg(2*ny,rx).*C36g(2*ny-1,rx);
    %for corner
    Cnew36g(1,1)=w0g(1,1).*C36g(1,1)+wxpg(1,1).*C36g(1,2)+wypg(1,1).*C36g(2,1);
    Cnew36g(1,nx)=w0g(1,nx).*C36g(1,nx)+wxmg(1,nx).*C36g(1,nx-1)+wypg(1,nx).*C36g(2,nx);
    Cnew36g(2*ny,1)=w0g(2*ny,1).*C36g(2*ny,1)+wxpg(2*ny,1).*C36g(2*ny,2)+wymg(2*ny,1).*C36g(2*ny-1,1);
    Cnew36g(2*ny,nx)=w0g(2*ny,nx).*C36g(2*ny,nx)+wxmg(2*ny,nx).*C36g(2*ny,nx-1)+wymg(2*ny,nx).*C36g(2*ny-1,nx);
    
    Cnew35g(ry,rx)=w0g(ry,rx).*C35g(ry,rx)+wxmg(ry,rx).*C35g(ry,rx-1)+wxpg(ry,rx).*C35g(ry,rx+1)+wymg(ry,rx).*C35g(ry-1,rx)+wypg(ry,rx).*C35g(ry+1,rx);
    %for each wall
    Cnew35g(ry,1)=w0g(ry,1).*C35g(ry,1)+wxpg(ry,1).*C35g(ry,2)+wymg(ry,1).*C35g(ry-1,1)+wypg(ry,1).*C35g(ry+1,1);
    Cnew35g(ry,nx)=w0g(ry,nx).*C35g(ry,nx)+wxmg(ry,nx).*C35g(ry,nx-1)+wymg(ry,nx).*C35g(ry-1,nx)+wypg(ry,nx).*C35g(ry+1,nx);
    Cnew35g(1,rx)=w0g(1,rx).*C35g(1,rx)+wxmg(1,rx).*C35g(1,rx-1)+wxpg(1,rx).*C35g(1,rx+1)+wypg(1,rx).*C35g(2,rx);
    Cnew35g(2*ny,rx)=w0g(2*ny,rx).*C35g(2*ny,rx)+wxmg(2*ny,rx).*C35g(2*ny,rx-1)+wxpg(2*ny,rx).*C35g(2*ny,rx+1)+wymg(2*ny,rx).*C35g(2*ny-1,rx);
    %for corner
    Cnew35g(1,1)=w0g(1,1).*C35g(1,1)+wxpg(1,1).*C35g(1,2)+wypg(1,1).*C35g(2,1);
    Cnew35g(1,nx)=w0g(1,nx).*C35g(1,nx)+wxmg(1,nx).*C35g(1,nx-1)+wypg(1,nx).*C35g(2,nx);
    Cnew35g(2*ny,1)=w0g(2*ny,1).*C35g(2*ny,1)+wxpg(2*ny,1).*C35g(2*ny,2)+wymg(2*ny,1).*C35g(2*ny-1,1);
    Cnew35g(2*ny,nx)=w0g(2*ny,nx).*C35g(2*ny,nx)+wxmg(2*ny,nx).*C35g(2*ny,nx-1)+wymg(2*ny,nx).*C35g(2*ny-1,nx);
    
% gyre respiration calculation   
    Cnew34g(Cnewg>0)=Cnew34g(Cnewg>0)-J34.*Cnew34g(Cnewg>0)./Cnewg(Cnewg>0);
    Cnew33g(Cnewg>0)=Cnew33g(Cnewg>0)-J33.*Cnew33g(Cnewg>0)./Cnewg(Cnewg>0);
    Cnew36g(Cnewg>0)=Cnew36g(Cnewg>0)-J36.*Cnew36g(Cnewg>0)./Cnewg(Cnewg>0);
    Cnew35g(Cnewg>0)=Cnew35g(Cnewg>0)-J35.*Cnew35g(Cnewg>0)./Cnewg(Cnewg>0);
    Cnewg=Cnewg-J1;
    
% equatorial respiration calculation
    Cnew34g(ty2,tx2)=Cnew34g(ty2,tx2)-J342.*Cnew34g(ty2,tx2)./Cnewg(ty2,tx2);
    Cnew33g(ty2,tx2)=Cnew33g(ty2,tx2)-J332.*Cnew33g(ty2,tx2)./Cnewg(ty2,tx2);
    Cnew36g(ty2,tx2)=Cnew36g(ty2,tx2)-J362.*Cnew36g(ty2,tx2)./Cnewg(ty2,tx2);
    Cnew35g(ty2,tx2)=Cnew35g(ty2,tx2)-J352.*Cnew35g(ty2,tx2)./Cnewg(ty2,tx2);
    Cnewg(ty2,tx2)=Cnewg(ty2,tx2)-J2;
    Cnew34g(ty2,tx3)=Cnew34g(ty2,tx3)-J342.*Cnew34g(ty2,tx3)./Cnewg(ty2,tx3);
    Cnew33g(ty2,tx3)=Cnew33g(ty2,tx3)-J332.*Cnew33g(ty2,tx3)./Cnewg(ty2,tx3);
    Cnew36g(ty2,tx3)=Cnew36g(ty2,tx3)-J362.*Cnew36g(ty2,tx3)./Cnewg(ty2,tx3);
    Cnew35g(ty2,tx3)=Cnew35g(ty2,tx3)-J352.*Cnew35g(ty2,tx3)./Cnewg(ty2,tx3);
    Cnewg(ty2,tx3)=Cnewg(ty2,tx3)-J2;
    Cnewg(Cnewg<0)=0;
    
% Set O2 concentration to zero if it is less than zero
    Cnew34g(Cnewg==0)=0;
    Cnew33g(Cnewg==0)=0;
    Cnew36g(Cnewg==0)=0;
    Cnew35g(Cnewg==0)=0;
    Cnew34g(Cnew34g<0)=0;
    Cnew33g(Cnew33g<0)=0;
    Cnew36g(Cnew36g<0)=0;
    Cnew35g(Cnew35g<0)=0;
    
% Set starting condtion for next iteration
    Cg=Cnewg;
    C34g=Cnew34g;
    C33g=Cnew33g;
    C36g=Cnew36g;
    C35g=Cnew35g;
end
time2=toc
tic
[C,C33,C34,C35,C36]=gather(Cg,C33g,C34g,C35g,C36g);
time3=toc