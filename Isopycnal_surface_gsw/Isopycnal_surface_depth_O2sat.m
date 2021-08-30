clear

min = 25.8;
max = 26.2;

%25.8-26.2%26.5-27.2

%read corresponding data from netcdf file
%O2 is in ml/l 

O2=ncread('woa13_all_o00_01.nc','o_an');
temp=ncread('woa13_decav_t00_01.nc','t_an');
S=ncread('woa13_decav_s00_01.nc','s_an');
depth_index=ncread('woa13_all_o00_01.nc','depth');

%do slicing for the data to only get the Pacific ocean data not including
%the subpolar gyre)
O2_pacific=NaN(360,180,102);
temp_pacific=NaN(360,180,102);
S_pacific=NaN(360,180,102);

a=24/21;
for i=0:21
    O2_pacific(ceil(331+a*i):360,135+i,1:102)=O2(ceil(331+a*i):360,135+i,1:102);
    temp_pacific(ceil(331+a*i):360,135+i,1:102)=temp(ceil(331+a*i):360,135+i,1:102);
    S_pacific(ceil(331+a*i):360,135+i,1:102)=S(ceil(331+a*i):360,135+i,1:102);
end
O2_pacific(1:80,107:156,1:102)=O2(1:80,107:156,1:102);
temp_pacific(1:80,107:156,1:102)=temp(1:80,107:156,1:102);
S_pacific(1:80,107:156,1:102)=S(1:80,107:156,1:102);



a=26/21;
for i=0:21
    O2_pacific(ceil(305+a*i):360,114+i,1:102)=O2(ceil(305+a*i):360,114+i,1:102);
    temp_pacific(ceil(305+a*i):360,114+i,1:102)=temp(ceil(305+a*i):360,114+i,1:102);
    S_pacific(ceil(305+a*i):360,114+i,1:102)=S(ceil(305+a*i):360,114+i,1:102);
end
O2_pacific(1:80,114:114+21,1:102)=O2(1:80,114:114+21,1:102);
temp_pacific(1:80,114:114+21,1:102)=temp(1:80,114:114+21,1:102);
S_pacific(1:80,114:114+21,1:102)=S(1:80,114:114+21,1:102);

O2_pacific(305:360,107:114,1:102)=O2(305:360,107:114,1:102);
temp_pacific(305:360,107:114,1:102)=temp(305:360,107:114,1:102);
S_pacific(305:360,107:114,1:102)=S(305:360,107:114,1:102);

O2_pacific(1:80,107:114,1:102)=O2(1:80,107:114,1:102);
temp_pacific(1:80,107:114,1:102)=temp(1:80,107:114,1:102);
S_pacific(1:80,107:114,1:102)=S(1:80,107:114,1:102);

O2_pacific(305:360,95:95+14,1:102)=O2(305:360,95:95+14,1:102);
temp_pacific(305:360,95:95+14,1:102)=temp(305:360,95:95+14,1:102);
S_pacific(305:360,95:95+14,1:102)=S(305:360,95:95+14,1:102);


a=22/12;
for i=0:12
    O2_pacific(1:(180-78-ceil(a*i)),95+i,1:102)=O2(1:(180-78-ceil(a*i)),95+i,1:102);
    temp_pacific(1:(180-78-ceil(a*i)),95+i,1:102)=temp(1:(180-78-ceil(a*i)),95+i,1:102);
    S_pacific(1:(180-78-ceil(a*i)),95+i,1:102)=S(1:(180-78-ceil(a*i)),95+i,1:102);
end
O2_pacific(305:360,90:95,1:102)=O2(305:360,90:95,1:102);
S_pacific(305:360,90:95,1:102)=S(305:360,90:95,1:102);
temp_pacific(305:360,90:95,1:102)=temp(305:360,90:95,1:102);
O2_pacific(1:102,90:95,1:102)=O2(1:102,90:95,1:102);
temp_pacific(1:102,90:95,1:102)=temp(1:102,90:95,1:102);
S_pacific(1:102,90:95,1:102)=S(1:102,90:95,1:102);

O2_pacific(305:360,24:90,1:102)=O2(305:360,24:90,1:102);
temp_pacific(305:360,24:90,1:102)=temp(305:360,24:90,1:102);
S_pacific(305:360,24:90,1:102)=S(305:360,24:90,1:102);

O2_pacific(1:110,24:90,1:102)=O2(1:110,24:90,1:102);
temp_pacific(1:110,24:90,1:102)=temp(1:110,24:90,1:102);
S_pacific(1:110,24:90,1:102)=S(1:110,24:90,1:102);



%calculate sigma theta for each latitude and longtitude
% O2_pacific=reshape(O2_pacific,[],1);
temp_pacific=reshape(temp_pacific,[],1);
S_pacific=reshape(S_pacific,[],1);


lon=(1:360)';
lon=repmat(lon,180*102,1);
lat=ones(360,1);
for i=2:180
    lat_t=i*ones(360,1);
    lat=cat(1,lat,lat_t);
end
lat=repmat(lat,102,1);
d=ones(360*180,1);
for i=2:102
    d_t=i*ones(360*180,1);
    d=cat(1,d,d_t);
end

depth=depth_index(d);%in m


%Calcaulate sigma_theta
rho=1.025*10^3; %kg/m^3
g=9.80665;
pressure=rho*g*depth/10^4;
lat=lat-90;
[SA_pacific, in_ocean] = gsw_SA_from_SP(S_pacific,pressure,lon,lat);
temp_con = gsw_CT_from_t(SA_pacific,temp_pacific,pressure);
sigma_theta = gsw_sigma0(SA_pacific,temp_con);
lat=reshape(lat,360,180,102);
lat_1=lat(:,:,1);
lat_1=reshape(lat_1,360,180);


p1 = 0.1;
p2 = 0.3;
p3 = 0.5;
p4 = 0.7;
p5 = 0.9;
sigma_theta=reshape(sigma_theta,360,180,102);
O2_isopycnal=NaN(360,180,102);
O2_isopycnal((sigma_theta<=max)&(sigma_theta>=min))=O2_pacific((sigma_theta<=max)&(sigma_theta>=min));
temp_isopycnal=NaN(360,180,102);
temp_isopycnal((sigma_theta<=max)&(sigma_theta>=min))=temp_pacific((sigma_theta<=max)&(sigma_theta>=min));
z=nanmean(O2_isopycnal,3);
z = z(lat_1>=-42 & lat_1<=42); %-53 50
lat_1 = lat_1(lat_1>=-42 & lat_1<=42);
% Normal one
% Area=sum(sum((~isnan(z).*cosd(lat_1))));
% total=nansum(nansum(z.*cosd(lat_1)));
% meanO2=total/Area;
% for 10-90%
Area=sum(sum((~isnan(z).*cosd(lat_1))));
total=nansum(nansum(z.*cosd(lat_1)));
meanO2=total/Area;



% Area=sum(sum(~isnan(z))); Calculate the relative area for O2 sat
O2_full=273/1.025/44.661;
% O2_full=317/1.025/44.661;
percent1=sum(sum((z>=p1*O2_full).*cosd(lat_1)))/Area*1000;
percent2=sum(sum((z>=p2*O2_full).*cosd(lat_1)))/Area*1000;
percent3=sum(sum((z>=p3*O2_full).*cosd(lat_1)))/Area*1000;
percent4=sum(sum((z>=p4*O2_full).*cosd(lat_1)))/Area*1000;
percent5=sum(sum((z>=p5*O2_full).*cosd(lat_1)))/Area*1000;
% 


z=NaN(360*180*102,1);
z((sigma_theta<=max)&(sigma_theta>=min))=O2_pacific((sigma_theta<=max)&(sigma_theta>=min));
y=NaN(360*180*102,1);
y((lat>=-45)&(lat<=-42))=z((lat>=-45)&(lat<=-42));
nanmean(y(:)); 
z=y(:);

figure
sf1 = subplot(1,2,1)
sf1.Position = [0.09 0.11 0.46 0.815]; 
z=depth((sigma_theta<=max)&(sigma_theta>=min));
size(z)
x=lon((sigma_theta<=max)&(sigma_theta>=min));
y=lat((sigma_theta<=max)&(sigma_theta>=min));
x(x>180)=x(x>180)-360;
gx(1)=geoscatter(y,x-180, 60,z,'Marker','.');


geobasemap none
t = text(55,-225,'A');
t.FontSize = 12;
t.FontWeight = 'bold';
set(gca,'FontSize',8,'LineWidth',1.5,'TickLength',[0.015, 0.025])
hcb1=colorbar;
hcb1.Limits=[0,550];
hcb1.TickLabels = {'','50','','','','','','','','','500',''};
hcb1.TickLength=0.01;
hcb1.LineWidth=1.5;
colormap(gca,cbrewer('div','RdBu',11));
title(hcb1,'Depth(m)','FontSize',7)


sf2 = subplot(1,2,2)
sf2.Position = [0.53 0.11 0.46 0.815];
z=O2_pacific((sigma_theta<=max)&(sigma_theta>=min))/O2_full*100;
x=lon((sigma_theta<=max)&(sigma_theta>=min));
y=lat((sigma_theta<=max)&(sigma_theta>=min));
x(x>180)=x(x>180)-360;
gx(2)=geoscatter(y,x-180,60,z,'Marker','.');
geobasemap none
t = text(55,-225,'B');
t.FontSize = 12;
t.FontWeight = 'bold';
geolimits(gx(2).Parent,[-60 60],[122 292]);
hold on
set(gca,'FontSize',8,'LineWidth',1.5,'TickLength',[0.015, 0.025])
hcb2=colorbar;
hcb2.Limits=[0,100];
hcb2.TickLabels = {'0','','','','','50','','','','','100'};
hcb2.TickLength=0.01;
hcb2.LineWidth=1.5;
colormap(gca,cbrewer('seq','YlOrRd',10));
caxis([0, 100]);
title(hcb2,'Dissolved O_{2}sat(%)','FontSize',7)

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4.5 3];


gx(2).Parent.LatitudeLabel.String = '';
gx(2).Parent.LatitudeAxis.TickLabels=[];
geolimits(gx(1).Parent,[-15 15],[-240 -70]);
geolimits(gx(2).Parent,[-15 15],[-240 -70]);

% print('fig3_isopycnal_plot','-djpeg','-r600')
% print('fig3_isopycnal_plot','-dpdf')
% savefig('fig3_isopycnal_plot.fig')
hold off

 
