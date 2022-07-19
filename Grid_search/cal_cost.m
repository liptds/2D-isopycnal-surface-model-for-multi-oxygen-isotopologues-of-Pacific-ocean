clear
O2_full = 273;
ans = [];
xm=1.05e7; ym=4.5e6; %domain dimensions
eps=0.015; %boundary current 
nx=ceil(10/eps); ny=ceil(nx*ym/xm); %rectangular grid
dx=xm/nx;dy=ym/ny; 
x=[0.5*dx:dx:xm-0.5*dx];
y=[-ym+0.5*dy:dy:ym-0.5*dy];
[X,Y]=meshgrid(x,y);
r = 6378.1*10^3;


%For general search
for A_run = 1.6:0.8:3.2
    for J_bulk = 0.5: 0.5: 3.0
        for J_eq = 5:10:55

        load(sprintf('O2satA%.5gJbulk%1.1f_Jeq%1.1f_larger_eq_2000kmwinter.mat',A_run,J_bulk,J_eq));
	  % cos is for latitude correction
        Area = sum(sum((C>=0) .* cos(Y/r)));
        percent1=sum(sum((C>=0.1*O2_full) .* cos(Y/r)))/Area*1000;
        percent2=sum(sum((C>=0.3*O2_full) .* cos(Y/r)))/Area*1000;
        percent3=sum(sum((C>=0.5*O2_full) .* cos(Y/r)))/Area*1000;
        percent4=sum(sum((C>=0.7*O2_full) .* cos(Y/r)))/Area*1000;
        percent5=sum(sum((C>=0.9*O2_full) .* cos(Y/r)))/Area*1000;
        %annual average the percentage number is coming from the WOA isopycnal output
        per0_diff = (1000-percent1) - (1000 - 986.3256);
        per1_diff = percent1 - percent2 - (986.3256-902.0042);%986.3256 902.0042 672.3308 458.4773 180.6421
        per2_diff = percent2 - percent3 - (902.0042-672.3308);
        per3_diff = percent3 - percent4 - (672.3308-458.4773);
        per4_diff = percent4 - percent5 - (458.4773-180.6421);
        per5_diff = percent5 - 180.6421;
        

        L1 = abs(per1_diff) + abs(per2_diff) + abs(per3_diff) + abs(per4_diff);
        L2 = per1_diff^2 + per2_diff^2 + per3_diff^2 + per4_diff^2;
        if isempty(ans)
            ans = [J_bulk,J_eq, A_run, per1_diff, per2_diff, per3_diff, per4_diff, L1, L2, per0_diff, per5_diff];
        else
            cat(1, ans, [J_bulk,J_eq, A_run, per1_diff, per2_diff, per3_diff, per4_diff, L1, L2, per0_diff, per5_diff]);
        end
        end
    end
end


new_ans = sortrows(ans,[9,8]);
writematrix(new_ans,'costGen.csv');