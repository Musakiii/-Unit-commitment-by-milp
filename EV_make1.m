%EV走行計画作成用プログラム
%{
20201210
fun_forciblycharge & fun_travelconsumptionを作った
とりあえずSOCの初期値とNconをできた
%}
%{
20201217
use rng function in case of unforeseen problems
%}

clear;
%一つのアグリゲータが1500台EV
EVnum=1500;
P_EVmax=6;
B_EV=60;
time=24*60*60;
%0時のSOC初期値を作る
%SOC初期値は30%~90%で平均分布
rng(1,'philox');
S_SOC_initial=rng;
SOC_initial=(0.9-0.3)*rand(EVnum,1)+0.3;
EVSOC=zeros(EVnum,time);
EVSOC(:,1)=SOC_initial;


%travel時刻を正規分布に作る
b1=8;
a1=1.5;
rng(2,'philox');
S_t_travel1=rng;
t_travel1=a1*randn(EVnum,1)+b1;
b2=18;
a2=2;
rng(3,'philox');
S_t_travel2=rng;
t_travel2=a2*randn(EVnum,1)+b2;

for ii=1:EVnum
    if t_travel2(ii)>24
        t_travel2(ii)=24; 
    end
end

%EVの状態表を作る
%0:controllable,1:forcibly,2:traveling

EVstate=zeros(EVnum,time);
t_travel1_second=fix(t_travel1*3600);
t_travel2_second=fix(t_travel2*3600);

for ii=1:EVnum
    EVstate(ii,t_travel1_second(ii))=2;
    EVstate(ii,t_travel2_second(ii))=2;
end

%EVの走行時間を作る
%0.4~2時間で平均分布
travel_duration=(2-0.4)*rand(EVnum,1)*60*60;
travel_duration=fix(travel_duration);

%travel時間更新（t_travel2_second(ii)+travel_duration(ii)<=timeのため）
for ii=1:EVnum
    if t_travel2_second(ii)+travel_duration(ii)>time
        travel_duration(ii)=time-t_travel2_second(ii); 
    end
end


%EVstate更新
%走行開始～終わり　２になる
%走行開始1時間前 1になる（forcibly）
for ii=1:EVnum
    EVstate(ii,t_travel1_second(ii):t_travel1_second(ii)+travel_duration(ii))=2;
    EVstate(ii,t_travel2_second(ii):t_travel2_second(ii)+travel_duration(ii))=2;
    EVstate(ii,t_travel1_second(ii)-3600:t_travel1_second(ii)-1)=1;
    EVstate(ii,t_travel2_second(ii)-3600:t_travel2_second(ii)-1)=1;
end

%test:平均走行時間
%testNUM=length(find(EVstate==2));
%testNUM=testNUM/3600/EVnum;

%EVSOC更新(初期値から)
%SOC初期値から一定の時間帯は把握できる
SOCtarget=0.75;
[EVSOC,initial_duration]=fun_forciblycharge(EVSOC,EVnum,ones(EVnum,1),P_EVmax,B_EV,SOCtarget);

%test2:SOCが75%になる時点はtravel時点超えるか
%{
test2vector=zeros(EVnum,1);
for ii=1:EVnum
    for jj=1:time
        if EVSOC(ii,jj)==0
            test2vector(ii)=jj;
            break;
        end
    end
end
compare=zeros(EVnum,1);
for ii=1:EVnum
    if test2vector(ii)>t_travel1_second
        compare(ii)=1; 
    end
end
test2NUM=length(find(compare==1));
%}


%強制充放電終わり時刻
%{
t_forcibly_end=zeros(EVnum,1);
for ii=1:EVnum
    for jj=1:time
        if EVSOC(ii,jj)==0
            t_forcibly_end(ii)=jj;
            break;
        end
    end
end
%}

%EVstate更新
for ii=1:EVnum
    EVstate(ii,1:initial_duration(ii))=1;
end


%走行消耗SOC更新
%EV走行時刻のSOCはとりあえず85%で
for ii=1:EVnum
    EVSOC(ii,t_travel1_second(ii))=0.85;
    EVSOC(ii,t_travel2_second(ii))=0.85;
end

%v=40km/h,cost=8km/kWh
velocity=40;
cost=8;
t_end1=t_travel1_second+travel_duration;
t_end2=t_travel2_second+travel_duration;
[EVSOC] = fun_travelconsumption(EVSOC,t_travel1_second,t_end1,velocity,cost,B_EV,EVnum);
[EVSOC] = fun_travelconsumption(EVSOC,t_travel2_second,t_end2,velocity,cost,B_EV,EVnum);

%EV到着後強制充放電
[EVSOC,travel1_duration]=fun_forciblycharge(EVSOC,EVnum,t_end1,P_EVmax,B_EV,SOCtarget);
[EVSOC,travel2_duration]=fun_forciblycharge(EVSOC,EVnum,t_end2,P_EVmax,B_EV,SOCtarget);

%EVstate更新(二回の走行後の強制充放電)
for ii=1:EVnum
    EVstate(ii,t_end1(ii)+1:t_end1(ii)+travel1_duration(ii))=1;
    EVstate(ii,t_end2(ii)+1:t_end2(ii)+travel2_duration(ii))=1; 
end

%各時間断面で制御可能EV台数を把握
%deltaT=0.5h
deltaT=0.5;
Ncon=zeros(24/deltaT,1);

for ii=1:EVnum
    for jj=1:24/deltaT
        if EVstate(ii,(jj-1)*deltaT*60*60+1)==0
            Ncon(jj)=Ncon(jj)+1; 
        end
    end
end



%% 値を保存する

%20201210
%とりあえずSOCの初期値とNconをできた、今後はこれらの値でUCをやってみる
%{
save('SOC_initial_20201210.mat','SOC_initial');
save('Ncon_20201210.mat','Ncon');
%}


%20201217
%rng関数加入、新データを作る
save('SOC_initial_20201217.mat','SOC_initial');
save('Ncon_20201217.mat','Ncon');
save('EV_make_data_20201217.mat');













