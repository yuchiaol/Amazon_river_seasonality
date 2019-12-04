% Final plot for trends analysis

clear;clc

load('final_trends.mat')
%%

% Plume Area
x = [1,1.5,1.75];
y1 = [trend_area_ctrl_c1,trend_area_1p50_c1,trend_area_1p75_c1];
y2 = [trend_area_ctrl_c2,trend_area_1p50_c2,trend_area_1p75_c2];
y3 = [trend_area_ctrl_c3,trend_area_1p50_c3,trend_area_1p75_c3];
y4 = [trend_area_ctrl_c4,trend_area_1p50_c4,trend_area_1p75_c4];
y5 = [trend_area_ctrl_c5,trend_area_1p50_c5,trend_area_1p75_c5];
ya1 = [y1(1),y2(1),y3(1),y4(1),y5(1)]; ya1 = mean(ya1);
ya2 = [y1(2),y2(2),y3(2),y4(2),y5(2)]; ya2 = mean(ya2);
ya3 = [y1(3),y2(3),y3(3),y4(3),y5(3)]; ya3 = mean(ya3);
ya = [ya1,ya2,ya3];

%clf(1)
figure(1)

%subplot(3,2,1)
plot(x,ya,'marker','o','color','k','linewidth',2)
hold on
scatter(x,y1,'ks')
scatter(x,y2,'k^')
scatter(x,y3,'kp')
scatter(x,y4,'kv')
scatter(x,y5,'kd')
box on 
grid on
xlim([0.5 2.0])
%ylim([0.1 0.2])
ylabel('(km^{2} per decade)')
ax = gca;
ax.FontSize=12;
ax.XTick = x;
ax.XTickLabel = {'Control','x1.5','x1.75'};
legend('average','cycle 1','cycle 2','cycle 3','cycle 4','cycle 5','Location','NorthWest')
title('(a) APR Area')
%
path = '/nara/data4/jsteffen/data/amazon_select/figs/';
name = [path 'plume_area_season_trend'];
print('-dpng',name,'-r300')
%}

%%

% Salinity
y1 = [trend_salt_ctrl_c1,trend_salt_1p50_c1,trend_salt_1p75_c1];
y2 = [trend_salt_ctrl_c2,trend_salt_1p50_c2,trend_salt_1p75_c2];
y3 = [trend_salt_ctrl_c3,trend_salt_1p50_c3,trend_salt_1p75_c3];
y4 = [trend_salt_ctrl_c4,trend_salt_1p50_c4,trend_salt_1p75_c4];
y5 = [trend_salt_ctrl_c5,trend_salt_1p50_c5,trend_salt_1p75_c5];
ya1 = [y1(1),y2(1),y3(1),y4(1),y5(1)]; ya1 = mean(ya1);
ya2 = [y1(2),y2(2),y3(2),y4(2),y5(2)]; ya2 = mean(ya2);
ya3 = [y1(3),y2(3),y3(3),y4(3),y5(3)]; ya3 = mean(ya3);
ya = [ya1,ya2,ya3];

%clf(2)
figure(2)

%subplot(3,2,2)
plot(x,ya,'marker','o','color','b','linewidth',2)
hold on
scatter(x,y1,'bs')
scatter(x,y2,'b^')
scatter(x,y3,'bp')
scatter(x,y4,'bv')
scatter(x,y5,'bd')
box on 
grid on
xlim([0.5 2.0])
%ylim([0.1 0.2])
ylabel('(PSU per decade)')
ax = gca;
ax.FontSize=12;
ax.XTick = x;
ax.XTickLabel = {'Control','x1.5','x1.75'};
legend('average','cycle 1','cycle 2','cycle 3','cycle 4','cycle 5','Location','NorthWest')
title('(b) APR 5 m Salinity')
%
path = '/nara/data4/jsteffen/data/amazon_select/figs/';
name = [path 'plume_salt_season_trend'];
print('-dpng',name,'-r300')
%}

%%
% Temp
y1 = [trend_temp_ctrl_c1,trend_temp_1p50_c1,trend_temp_1p75_c1];
y2 = [trend_temp_ctrl_c2,trend_temp_1p50_c2,trend_temp_1p75_c2];
y3 = [trend_temp_ctrl_c3,trend_temp_1p50_c3,trend_temp_1p75_c3];
y4 = [trend_temp_ctrl_c4,trend_temp_1p50_c4,trend_temp_1p75_c4];
y5 = [trend_temp_ctrl_c5,trend_temp_1p50_c5,trend_temp_1p75_c5];
ya1 = [y1(1),y2(1),y3(1),y4(1),y5(1)]; ya1 = mean(ya1);
ya2 = [y1(2),y2(2),y3(2),y4(2),y5(2)]; ya2 = mean(ya2);
ya3 = [y1(3),y2(3),y3(3),y4(3),y5(3)]; ya3 = mean(ya3);
ya = [ya1,ya2,ya3];

%clf(3)
figure(3)
%subplot(3,2,3)
plot(x,ya,'marker','o','color','r','linewidth',2)
hold on
scatter(x,y1,'rs')
scatter(x,y2,'r^')
scatter(x,y3,'rp')
scatter(x,y4,'rv')
scatter(x,y5,'rd')
box on 
grid on
xlim([0.5 2.0])
%ylim([0.1 0.2])
ylabel(['(' char(176) 'C per decade)'])
ax = gca;
ax.FontSize=12;
ax.XTick = x;
ax.XTickLabel = {'Control','x1.5','x1.75'};
legend('average','cycle 1','cycle 2','cycle 3','cycle 4','cycle 5','Location','NorthWest')
title('   (c) APR 5 m Temperature')
%
path = '/nara/data4/jsteffen/data/amazon_select/figs/';
name = [path 'plume_temp_season_trend'];
print('-dpng',name,'-r300')
%}

%%

% Buoyancy frequency squared
y1 = [trend_nsqr_ctrl_c1,trend_nsqr_1p50_c1,trend_nsqr_1p75_c1];
y2 = [trend_nsqr_ctrl_c2,trend_nsqr_1p50_c2,trend_nsqr_1p75_c2];
y3 = [trend_nsqr_ctrl_c3,trend_nsqr_1p50_c3,trend_nsqr_1p75_c3];
y4 = [trend_nsqr_ctrl_c4,trend_nsqr_1p50_c4,trend_nsqr_1p75_c4];
y5 = [trend_nsqr_ctrl_c5,trend_nsqr_1p50_c5,trend_nsqr_1p75_c5];
ya1 = [y1(1),y2(1),y3(1),y4(1),y5(1)]; ya1 = mean(ya1);
ya2 = [y1(2),y2(2),y3(2),y4(2),y5(2)]; ya2 = mean(ya2);
ya3 = [y1(3),y2(3),y3(3),y4(3),y5(3)]; ya3 = mean(ya3);
ya = [ya1,ya2,ya3];

%clf(4)
figure(4)
%subplot(3,2,4)
plot(x,ya,'marker','o','color',[0.466 0.674 0.188],'linewidth',2)
hold on
scatter(x,y1,[],[0.466 0.674 0.188],'s')
scatter(x,y2,[],[0.466 0.674 0.188],'^')
scatter(x,y3,[],[0.466 0.674 0.188],'p')
scatter(x,y4,[],[0.466 0.674 0.188],'v')
scatter(x,y5,[],[0.466 0.674 0.188],'d')
box on 
grid on
xlim([0.5 2.0])
%ylim([0.1 0.2])
ylabel('(s^{-2} per decade)')
ax = gca;
ax.FontSize=12;
ax.XTick = x;
ax.XTickLabel = {'Control','x1.5','x1.75'};
legend('average','cycle 1','cycle 2','cycle 3','cycle 4','cycle 5','Location','NorthWest')
title('(d) APR N^{2}')
%
path = '/nara/data4/jsteffen/data/amazon_select/figs/';
name = [path 'plume_nsqr_season_trend'];
print('-dpng',name,'-r300')
%}

%%
% BLT
y1 = [trend_blt_ctrl_c1,trend_blt_1p50_c1,trend_blt_1p75_c1];
y2 = [trend_blt_ctrl_c2,trend_blt_1p50_c2,trend_blt_1p75_c2];
y3 = [trend_blt_ctrl_c3,trend_blt_1p50_c3,trend_blt_1p75_c3];
y4 = [trend_blt_ctrl_c4,trend_blt_1p50_c4,trend_blt_1p75_c4];
y5 = [trend_blt_ctrl_c5,trend_blt_1p50_c5,trend_blt_1p75_c5];
ya1 = [y1(1),y2(1),y3(1),y4(1),y5(1)]; ya1 = mean(ya1);
ya2 = [y1(2),y2(2),y3(2),y4(2),y5(2)]; ya2 = mean(ya2);
ya3 = [y1(3),y2(3),y3(3),y4(3),y5(3)]; ya3 = mean(ya3);
ya = [ya1,ya2,ya3];

%clf(5)
figure(5)
%subplot(3,2,5)
plot(x,ya,'marker','o','color','k','linewidth',2)
hold on
scatter(x,y1,'ks')
scatter(x,y2,'k^')
scatter(x,y3,'kp')
scatter(x,y4,'kv')
scatter(x,y5,'kd')
box on 
grid on
xlim([0.5 2.0])
%ylim([0.1 0.2])
ylabel('(m per decade)')
ax = gca;
ax.FontSize=12;
ax.XTick = x;
ax.XTickLabel = {'Control','x1.5','x1.75'};
ax.YTick = [1.4,1.5,1.6,1.7,1.8,1.9];
legend('average','cycle 1','cycle 2','cycle 3','cycle 4','cycle 5','Location','NorthWest')
title('(e) APR Barrier Layer Thickness')
%
path = '/nara/data4/jsteffen/data/amazon_select/figs/';
name = [path 'plume_blt_season_trend'];
print('-dpng',name,'-r300')
%}

%%
% BLPE
y1 = [trend_blpe_ctrl_c1,trend_blpe_1p50_c1,trend_blpe_1p75_c1];
y2 = [trend_blpe_ctrl_c2,trend_blpe_1p50_c2,trend_blpe_1p75_c2];
y3 = [trend_blpe_ctrl_c3,trend_blpe_1p50_c3,trend_blpe_1p75_c3];
y4 = [trend_blpe_ctrl_c4,trend_blpe_1p50_c4,trend_blpe_1p75_c4];
y5 = [trend_blpe_ctrl_c5,trend_blpe_1p50_c5,trend_blpe_1p75_c5];
ya1 = [y1(1),y2(1),y3(1),y4(1),y5(1)]; ya1 = mean(ya1);
ya2 = [y1(2),y2(2),y3(2),y4(2),y5(2)]; ya2 = mean(ya2);
ya3 = [y1(3),y2(3),y3(3),y4(3),y5(3)]; ya3 = mean(ya3);
ya = [ya1,ya2,ya3];

%clf(6)
figure(6)
%subplot(3,2,6)
plot(x,ya,'marker','o','color','k','linewidth',2)
hold on
scatter(x,y1,'ks')
scatter(x,y2,'k^')
scatter(x,y3,'kp')
scatter(x,y4,'kv')
scatter(x,y4,'kd')
box on 
grid on
xlim([0.5 2.0])
%ylim([0.1 0.2])
ylabel('(J m^{-2} per decade)')
ax = gca;
ax.FontSize=12;
ax.XTick = x;
ax.XTickLabel = {'Control','x1.5','x1.75'};
legend('average','cycle 1','cycle 2','cycle 3','cycle 4','cycle 5','Location','NorthWest')
title('(f) APR Barrier Layer Potential Energy')
%
path = '/nara/data4/jsteffen/data/amazon_select/figs/';
name = [path 'plume_blpe_season_trend'];
print('-dpng',name,'-r300')
%}


%{
path = '/nara/data4/jsteffen/data/amazon_select/figs/';
name = [path 'plume_final_season_trend'];
print('-dpng',name,'-r300')
%}



