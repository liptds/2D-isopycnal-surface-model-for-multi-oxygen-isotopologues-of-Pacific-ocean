ry=[2:2*ny-1]; rx=[2:nx-1];
%Z=zeros(2*ny,nx,2); 
%Csto=cat(3,Csto,Z);
%C33sto=cat(3,C33sto,Z);
%C34sto=cat(3,C34sto,Z);
%C35sto=cat(3,C35sto,Z);
%C36sto=cat(3,C36sto,Z);

tic
%nx=gpuArray(nx);ny=gpuArray(ny);
%ry=[2:2*ny-1]; rx=[2:nx-1];
% nyr=1;
% nyrg=gpuArray(1);
% ntyrg=gpuArray(ntyr);
%nwg=gpuArray(nw);
%ntg=gpuArray(nt);
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
% Jg=gpuArray(J);
% J33g=gpuArray(J33);
% J34g=gpuArray(J34);
% J35g=gpuArray(J35);
% J36g=gpuArray(J36);
% Cstog=gpuArray(Csto);
% C33stog=gpuArray(C33sto);
% C34stog=gpuArray(C34sto);
% C35stog=gpuArray(C35sto);
% C36stog=gpuArray(C36sto);
% Cboundg=gpuArray(Cbound);
%ntyrg=gpuArray(ntyr);
%zero=zeros(2*ny,nx,10);
%Csto=cat(3,Csto,zero);
%C33sto=cat(3,C33sto,zero);
%C34sto=cat(3,C34sto,zero);
%C35sto=cat(3,C35sto,zero);
%C36sto=cat(3,C36sto,zero);
time1=toc
tic
for i=1:nt
    %conservative test
    
    %boundary condition
    j=1:nw;
    Cg(j,:)=Cbound;
    C33g(j,:)=Cbound*R17eq;
    C34g(j,:)=Cbound*R18eq;
    C35g(j,:)=Cbound*R35eq;
    C36g(j,:)=Cbound*R36eq;
    %C(j,:)=Cbound;
    k=1:nw2;
    Cg(571:572,k)=Cbound;
    C33g(571:572,k)=Cbound*R17eq;
    C34g(571:572,k)=Cbound*R18eq;
    C35g(571:572,k)=Cbound*R35eq;
    C36g(571:572,k)=Cbound*R36eq;
    
    % photosynthesis from SPOT 500km-800km ont the east coast of Pacific
%     Cg(559:568,(1:20))=Cbound*1.6;%158:175,(nx-3:nx)%138:195,(nx-30):nx
%     C33g(559:568,(1:20))=Cbound*1*R17eq+Cbound*0.6*R17p;
%     C34g(559:568,(1:20))=Cbound*1*R18eq+Cbound*0.6*R18p;
%     C35g(559:568,(1:20))=Cbound*1*R35eq+Cbound*0.6*R35p;
%     C36g(559:568,(1:20))=Cbound*1*R36eq+Cbound*0.6*R36p;
%     Cg(158:175,(nx-3:nx))=Cbound*1.04;%158:175,(nx-3:nx)%138:195,(nx-30):nx
%     C33g(158:175,(nx-3:nx))=Cbound*1.04*R17p;
%     C34g(158:175,(nx-3:nx))=Cbound*1.04*R18p;
%     C35g(158:175,(nx-3:nx))=Cbound*1.04*R35p;
%     C36g(158:175,(nx-3:nx))=Cbound*1.04*R36p;
%     k=170:(170+np);
%     Cg(k,nx)=138.13;%from SPOT
%     C33g(k,nx)=138.13*1.004479594;%from SPOT
%     C34g(k,nx)=138.13*1.008491342;
%     C35g(k,nx)=138.13*1.014135913;
%     C36g(k,nx)=138.13*1.019234189;
    
    %peru high photo
%     l=51:67;
%     Cg(l,(nx-5):nx)=241*1.1030;
%     C33g(l,(nx-5):nx)=241*1.1030*R17mix;
%     C34g(l,(nx-5):nx)=241*1.1030*R18mix;
%     C35g(l,(nx-5):nx)=241*1.1030*R35mix;
%     C36g(l,(nx-5):nx)=241*1.1030*R36mix;

%     Cg(l,nx)=241*1.1030;
%     C33g(l,nx)=241*1.1030*R17mix;
%     C34g(l,nx)=241*1.1030*R18mix;
%     C35g(l,nx)=241*1.1030*R35mix;
%     C36g(l,nx)=241*1.1030*R36mix;
    
    %regular iteration
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
    
%     %respiration
    
    %respiration is related to (old code slow)
    %Cnew34g(Cnew>0&Cnew34g>J34.*Cnew34./Cnew)=Cnew34(Cnew>0&Cnew34>J34.*Cnew34./Cnew)-J34.*Cnew34(Cnew>0&Cnew34>J34.*Cnew34./Cnew)./Cnew(Cnew>0&Cnew34>J34.*Cnew34./Cnew);
    %Cnew34(Cnew>0&Cnew34<=J34.*Cnew34./Cnew)=0;
    %Cnew33(Cnew>0&Cnew33>J33.*Cnew33./Cnew)=Cnew33(Cnew>0&Cnew33>J33.*Cnew33./Cnew)-J33.*Cnew33(Cnew>0&Cnew33>J33.*Cnew33./Cnew)./Cnew(Cnew>0&Cnew33>J33.*Cnew33./Cnew);
    %Cnew33(Cnew>0&Cnew33<=J33.*Cnew33./Cnew)=0;
    %Cnew36(Cnew>0&Cnew36>J36.*Cnew36./Cnew)=Cnew36(Cnew>0&Cnew36>J36.*Cnew36./Cnew)-J36.*Cnew36(Cnew>0&Cnew36>J36.*Cnew36./Cnew)./Cnew(Cnew>0&Cnew36>J36.*Cnew36./Cnew);
    %Cnew36(Cnew>0&Cnew36<=J36.*Cnew36./Cnew)=0;
    %Cnew35(Cnew>0&Cnew35>J35.*Cnew33./Cnew)=Cnew35(Cnew>0&Cnew35>J35.*Cnew33./Cnew)-J35.*Cnew33(Cnew>0&Cnew35>J35.*Cnew33./Cnew)./Cnew(Cnew>0&Cnew35>J35.*Cnew33./Cnew);
    %Cnew35(Cnew>0&Cnew35<=J35.*Cnew35./Cnew)=0;
    
    %respiration new code for faster
    Cnew34g(Cnewg>0)=Cnew34g(Cnewg>0)-J34.*Cnew34g(Cnewg>0)./Cnewg(Cnewg>0);
    Cnew33g(Cnewg>0)=Cnew33g(Cnewg>0)-J33.*Cnew33g(Cnewg>0)./Cnewg(Cnewg>0);
    Cnew36g(Cnewg>0)=Cnew36g(Cnewg>0)-J36.*Cnew36g(Cnewg>0)./Cnewg(Cnewg>0);
    Cnew35g(Cnewg>0)=Cnew35g(Cnewg>0)-J35.*Cnew35g(Cnewg>0)./Cnewg(Cnewg>0);
    Cnewg=Cnewg-J;
    Cnewg(Cnewg<0)=0;
%     
    Cnew34g(Cnewg==0)=0;
    Cnew33g(Cnewg==0)=0;
    Cnew36g(Cnewg==0)=0;
    Cnew35g(Cnewg==0)=0;
    Cnew34g(Cnew34g<0)=0;
    Cnew33g(Cnew33g<0)=0;
    Cnew36g(Cnew36g<0)=0;
    Cnew35g(Cnew35g<0)=0;
    
    %Cnew34(Cnew==0|Cnew34<0)=0;
    %Cnew33(Cnew==0|Cnew33<0)=0;
    %Cnew36(Cnew==0|Cnew36<0)=0;
    %Cnew35(Cnew==0|Cnew35<0)=0;
    %Cnew(Cnew>J)=Cnew(Cnew>J)-J;
    %Cnew(Cnew<=J)=0;
    
    %Cnew34(Cnew==0)=0;
    %Cnew33(Cnew==0)=0;
    %Cnew36(Cnew==0)=0;
    %Cnew35(Cnew==0)=0;
  
    
    
    
    Cg=Cnewg;
    C34g=Cnew34g;
    C33g=Cnew33g;
    C36g=Cnew36g;
    C35g=Cnew35g;
    %Sums(i,1)=sum(sum(C));
   
%     if floor(i/ntyr/50)>=nyr
%         
%         Cstog(:,:,nyr)=Cg(:,:); 
%         C34stog(:,:,nyr)=C34g(:,:); 
%         C33stog(:,:,nyr)=C33g(:,:); 
%         C35stog(:,:,nyr)=C35g(:,:); 
%         C36stog(:,:,nyr)=C36g(:,:); 
%         nyr=nyr+1; 
%         %Sums(nyr+1,1)=sum(sum(C));
%     end

    
end
time2=toc
tic
% [C]=gather(Cg);
[C,C33,C34,C35,C36]=gather(Cg,C33g,C34g,C35g,C36g);
% [Csto,C33sto,C34sto,C35sto,C36sto]=gather(Cstog,C33stog,C34stog,C35stog,C36stog);
time3=toc