clear 
%load('Full resolutionJ2.mat')
PsiDef
VelCalc
WeightCalc
InitCon
AdvDiffGPU



d18O=NaN(size(C));
D17O=NaN(size(C));
D36=NaN(size(C));
d18O(C>0)=(C34(C>0)./C(C>0)-1)*1000;
D17O(C>0)=1e6*(log(C33(C>0)./C(C>0))-0.518*log(C34(C>0)./C(C>0)));
D36(C>0)=(C36(C>0)./C(C>0)./(C34(C>0)./C(C>0)).^2-1)*1000;
D35=(C35./C./(C34./C)./(C33./C)-1)*1000;
O2sat=C/Cbound;

%load('Sampleoceandata.mat')


figure
[c,h]=contour(x,y,O2sat);
clabel(c,h)

%figure
%[c,h]=contour(x,y,D17O);
%clabel(c,h)

%figure
%[c,h]=contour(x,y,D36);
%clabel(c,h)

%figure
%[c,h]=contour(x,y,D35);
%clabel(c,h)

%figure
%plot(O2sat,d18O,'Marker','.','Color','red')
%hold on
%axis([0 1.1 -5 20])
%plot(SPOT91416(:,18),SPOT91416(:,19),'o','MarkerSize',8)
%plot(SPOT41217(:,18),SPOT41217(:,19),'o','MarkerSize',8)
%plot(CDISK4S1(:,18),CDISK4S1(:,19),'o','MarkerSize',8)
%plot(CDISK4S2(:,18),CDISK4S2(:,19),'o','MarkerSize',8)
%plot(CDISK4S3(:,18),CDISK4S3(:,19),'o','MarkerSize',8)
%plot(CDISK4S4(:,18),CDISK4S4(:,19),'o','MarkerSize',8)
%plot(CDISK4S5(:,18),CDISK4S5(:,19),'o','MarkerSize',8)
%plot(CDISK4S7(:,18),CDISK4S7(:,19),'o','MarkerSize',8)
%xlabel('O2sat')
%ylabel('d18O')
%plot(O2sat,D17O,'.')
%hold
%axis([0.05 1 -50 200])
%figure
%plot(O2sat,D17O,'Marker','.','Color','red')
%hold on
%axis([0 1.1 -5 20])
%axis([0 1.1 0 150])
%plot(SPOT91416(:,18),SPOT91416(:,20),'o','MarkerSize',8)
%plot(SPOT41217(:,18),SPOT41217(:,20),'o','MarkerSize',8)
%plot(CDISK4S1(:,18),CDISK4S1(:,20),'o','MarkerSize',8)
%plot(CDISK4S2(:,18),CDISK4S2(:,20),'o','MarkerSize',8)
%plot(CDISK4S3(:,18),CDISK4S3(:,20),'o','MarkerSize',8)
%plot(CDISK4S4(:,18),CDISK4S4(:,20),'o','MarkerSize',8)
%plot(CDISK4S5(:,18),CDISK4S5(:,20),'o','MarkerSize',8)
%plot(CDISK4S7(:,18),CDISK4S7(:,20),'o','MarkerSize',8)
%xlabel('O2sat')
%ylabel('D17O')


%figure
%plot(O2sat,D36,'Marker','.','Color','red')
%hold on
%axis([0 1.1 1.4 3.3])
%plot(SPOT91416(:,18),SPOT91416(:,21),'o','MarkerSize',8)
%plot(SPOT41217(:,18),SPOT41217(:,21),'o','MarkerSize',8)
%plot(CDISK4S1(:,18),CDISK4S1(:,21),'o','MarkerSize',8)
%plot(CDISK4S2(:,18),CDISK4S2(:,21),'o','MarkerSize',8)
%plot(CDISK4S3(:,18),CDISK4S3(:,21),'o','MarkerSize',8)
%plot(CDISK4S4(:,18),CDISK4S4(:,21),'o','MarkerSize',8)
%plot(CDISK4S5(:,18),CDISK4S5(:,21),'o','MarkerSize',8)
%plot(CDISK4S7(:,18),CDISK4S7(:,21),'o','MarkerSize',8)
%xlabel('O2sat')
%ylabel('D36')
%figure
%plot(O2sat,D36)

%figure
%plot(O2sat,D35)

%figure
%plot(d18O,D17O)
%figure
%plot(D36,D17O,'Marker','.','Color','red')
%hold on
%axis([0 1.1 1.4 3.3])
%plot(SPOT91416(:,21),SPOT91416(:,20),'o','MarkerSize',8,'color','green')
%plot(SPOT41217(:,21),SPOT41217(:,20),'o','MarkerSize',8,'color','green')
%plot(CDISK4S1(:,21),CDISK4S1(:,20),'o','MarkerSize',8)
%plot(CDISK4S2(:,21),CDISK4S2(:,20),'o','MarkerSize',8)
%plot(CDISK4S3(:,21),CDISK4S3(:,20),'o','MarkerSize',8)
%plot(CDISK4S4(:,21),CDISK4S4(:,20),'o','MarkerSize',8)
%plot(CDISK4S5(:,21),CDISK4S5(:,20),'o','MarkerSize',8)
%plot(CDISK4S7(:,21),CDISK4S7(:,20),'o','MarkerSize',8)
%xlabel('D36')
%ylabel('D17O')
%figure
%plot(d18O,D36)

%figure
%plot(1:10,M)
%hold on
%axis([1 11 -2000 0])

%M=(sum(sum(Csto(:,:,1:10)))/sum(sum(Csto(:,:,10)))-1)*10^6;
%M=reshape(M,[1,10]);
