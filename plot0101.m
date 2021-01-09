%UC visualization
%cooperate with milp0101
%{
20210101
PV 1000MW

20201230
output the figure for 全国大会

20201216
UC_result=[timeline;fire;hydro;nuclear;P_PV;P_WT;-P_EV;rest_demand;surplus;P_demand;Ncon]
%}
clear;
load('UC_result_20210101.mat');


tiledlayout(2,1);
ax1 = nexttile;


area(ax1,timeline',[hydro;nuclear;fire;P_WT;P_PV]');

ax1.XLim = [0 23.5];
xlabel(ax1,'Time[hour]');
ylabel(ax1,'Power[MW]');
xticks(ax1,[0:5:23.5]);
hold on;
plot(timeline',P_demand','--','LineWidth',2);
hold on;
plot(timeline',P_EV','--','LineWidth',2);
hold on
legend('Hydro','Nuclear','Thermal','WT','PV','Load','EV');

ax2 = nexttile;
plot(ax2,timeline,Ncon);
xlabel(ax2,'Time[hour]');
ylabel(ax2,'Number');
xticks(ax2,[0:5:23.5]);








