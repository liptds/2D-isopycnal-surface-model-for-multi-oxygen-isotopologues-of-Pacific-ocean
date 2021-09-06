ry=[2:2*ny-1]; rx=[2:nx-1];


tic

Cg=gpuArray(C);
Cnewg=gpuArray(Cnew);

w0g=gpuArray(w0);
wxmg=gpuArray(wxm);
wxpg=gpuArray(wxp);
wymg=gpuArray(wym);
wypg=gpuArray(wyp);

time1=toc
tic
for i=1:nt
    %conservative test
    
    %boundary condition
    j=1:nw;
    Cg(j,:)=Cbound;

    %C(j,:)=Cbound;
    k=1:nw2;
    Cg(571:572,k)=Cbound;

    
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
    
 
    Cnewg=Cnewg-J;
    Cnewg(Cnewg<0)=0;

    
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
[C]=gather(Cg);

% [Csto,C33sto,C34sto,C35sto,C36sto]=gather(Cstog,C33stog,C34stog,C35stog,C36stog);
time3=toc