clear;
load('P_load20201229.mat');

ax1=nexttile;
x=[1:48]/2;
y=testdata';
x1=[1:0.01:48]/2;
y1=interp1(x,y,x1,'spline');
plot(x1,y1);

ax1.XLim = [0 24];
xlabel(ax1,'Time[hour]');
ylabel(ax1,'Power[MW]');
xticks(ax1,[0:5:24]);