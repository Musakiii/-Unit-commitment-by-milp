%EVSOC visualization
%cooperate with UCaggregator_20210101
%{
20210101ver2


20210101
output the figure for 全国大会

%}
clear;
load('EVSOC20210101ver2.mat');

ax1=nexttile;
timeline=[24/86400:24/86400:24];

plot(timeline,EVSOC(2,:));
hold on;
plot(timeline,EVSOC(20,:));
hold on;
plot(timeline,EVSOC(50,:));
hold on;
plot(timeline,EVSOC(33,:));
hold on;
plot(timeline,EVSOC(90,:));
hold on;
plot(timeline,EVSOC(300,:));
hold on;
plot(timeline,EVSOC(401,:));
hold on;
plot(timeline,EVSOC(402,:));
hold on;
plot(timeline,EVSOC(502,:));
hold on;
plot(timeline,EVSOC(600,:));
hold on;
plot(timeline,EVSOC(777,:));
hold on;
plot(timeline,EVSOC(888,:));
hold on;
plot(timeline,EVSOC(999,:));
hold on;
plot(timeline,EVSOC(1200,:));
hold on;
plot(timeline,EVSOC(1111,:));
hold on;


ax1.XLim=[0 23.5];
xlabel(ax1,'Time[hour]');
ylabel(ax1,'SOC');
xticks(ax1,[0:5:23.5]);
yticks(ax1,[0.3:0.05:1]);


