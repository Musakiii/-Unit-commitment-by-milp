%% UC by intlinprog
%{
20201214
cplex cannot operate on matlab2020, so try to code a intlinprog version
%}


clear

%% data load
load('PVdata_20201210.mat');
load('load_weekday_20201211.mat');
load('Ncon_20201210.mat');


%% constant
E_oneEV=0.06;
P_oneEV=0.004;
agg_EVnum=1500;
Time=48;
width_num=4;
agg_num=20;
P_WTrate=200;
TG_LFC=0.2;%LFC capacity?

%% generator data
%PV=630,WT=200,thermal=2300,nuclear=800,hydro=200,NaS=300
%demand=3200
PVdata=PVdata*100;
load_weekday(2,:)=load_weekday(2,:)*0.22;
P_nuclear=600;
P_hydro=200;
P_load=load_weekday(2,0.25*36000:0.5*36000:864001)';

%% thermal data
%{
build a thermal matrix
Pmax,Pmin,a,b,c,startupcost
%}

thermal_data=[
60 80 304 300 620 591 350;
30 40 152 150 310 295.5 175;
0.0646 0.025 0.0533 0.0224 0.0067 0.0081 0.0032;
16.8899 43.5 9.2374 19.7128 10.2202 20.2909 10.102;
49.7108 200 164.3492 287.4982 207.1786 378.5918 350.6398;
181.6416 13.356 111.4208 1511.8992 1777.1544 2070.18 8331.9264
];

%quadratic 2 linear
dT=0.5;
[alpha,beta,P_bound] = func_linearized(thermal_data,width_num,dT);







%% PV WT TG 
P_PVmax=PVdata(15:30:1440)';
P_WTmax=ones(48,1)*P_WTrate;
P_TGmax=thermal_data(1,:)';
P_PVmin=zeros(48,1);
P_WTmin=zeros(48,1);
P_TGmin=thermal_data(2,:)';
S0=thermal_data(6,:)';%thermal startup cost
dPre0=P_TGmax*0.05;%LFC 5%
Gen_N=size(thermal_data,2);




%% variable
%TG
clearTG_level = zeros(Time*width_num,Gen_N);%閬嬭虎銉儥銉倰瑜囨暟鐢ㄦ剰銇嬨亼銈嬫檪闁撳箙0or1
clearP_TG     = zeros(Time*width_num,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
clearS        = zeros(Time,Gen_N);%璧峰嫊銉曘儵銉冦偘0or1

%NaS
%20201211 NaS銇俱仛灏庡叆銇椼仾銇?
%{
clearNaS_plus    = zeros(Time,1);%鏀鹃浕銇仺銇嶃伀锛?
clearNaS_minus   = zeros(Time,1);%鍏呴浕銇仺銇嶃伀锛?
clearP_NaS_plus  = zeros(Time,1);%鏀鹃浕銇仺銇嶃伄鍑哄姏
clearP_NaS_minus = zeros(Time,1);%鍏呴浕銇檪銇嚭鍔?
%}

%HP
%20201211 HP銇俱仛灏庡叆銇椼仾銇?
%{
clearHP_on = zeros(Time,HP_Group);%銇濄伄鏅傞枔甯伀璧峰嫊銇欍倠銇倝锛戯紝銈般儷銉笺儣鏁般伄鏁般仩銇戠敤鎰?0or1
%}

%TG銇甃FC澶夋暟
clearTG_LFCup = zeros(Time,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
clearTG_LFCdown = zeros(Time,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
%clearNaS_LFCup = zeros(Time,1);
%clearNaS_LFCdown = zeros(Time,1);

%PV銇╓T闁總銇鏁?
%level:1銇?-1銇?
%20201211 level銇勩倝銇亜
%{
clearPV_level = zeros(Time*4,1);
clearWT_level = zeros(Time*4,1);
%}
clearP_PV = zeros(Time,1);
clearP_WT = zeros(Time,1);


%20201211鎻愭鎵嬫硶銇枹銇欍倠澶夋暟
clearP_EV=zeros(Time,1);
clearE_EV=zeros(Time,1);



%% lower?
%TG
lb_TG_level = zeros(Time*width_num,Gen_N);%閬嬭虎銉儥銉倰瑜囨暟鐢ㄦ剰銇嬨亼銈嬫檪闁撳箙0or1
lb_P_TG     = zeros(Time*width_num,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
lb_S        = zeros(Time,Gen_N);%璧峰嫊銉曘儵銉冦偘0or1


%TG LFC
lb_TG_LFCup = zeros(Time,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
lb_TG_LFCdown = zeros(Time,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
%lb_NaS_LFCup = zeros(Time,1);
%lb_NaS_LFCdown = zeros(Time,1);

%PV WT

%{
lb_PV_level = zeros(Time*4,1);
lb_WT_level = zeros(Time*4,1);
%}
lb_P_PV = zeros(Time,1);
lb_P_WT = zeros(Time,1);


% EV P and E
lb_P_EV=-1*P_oneEV*Ncon*agg_num*agg_EVnum;
lb_E_EV=0.65*E_oneEV*Ncon*agg_num*agg_EVnum;

lb=[lb_TG_level(:);lb_P_TG(:);lb_S(:);lb_TG_LFCup(:);lb_TG_LFCdown(:);lb_P_PV(:);lb_P_WT(:);lb_P_EV(:);lb_E_EV(:)];
lb(3841)=0;


%% upper
%TG
ub_TG_level = ones(Time*width_num,Gen_N);%閬嬭虎銉儥銉倰瑜囨暟鐢ㄦ剰銇嬨亼銈嬫檪闁撳箙0or1
ub_P_TG     = inf(Time*width_num,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
ub_S        = ones(Time,Gen_N);%璧峰嫊銉曘儵銉冦偘0or1


%TG LFC
ub_TG_LFCup = inf(Time,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
ub_TG_LFCdown = inf(Time,Gen_N);%鍚勯亱杌€儸銉欍儷銇銇椼仸鍑哄姏
%lb_NaS_LFCup = zeros(Time,1);
%lb_NaS_LFCdown = zeros(Time,1);

%PV
%{
lb_PV_level = zeros(Time*4,1);
lb_WT_level = zeros(Time*4,1);
%}
ub_P_PV = P_PVmax;
ub_P_WT = ones(Time,1)*P_WTrate;


%
ub_P_EV=P_oneEV*Ncon*agg_num*agg_EVnum;
ub_E_EV=0.85*E_oneEV*Ncon*agg_num*agg_EVnum;

ub=[ub_TG_level(:);ub_P_TG(:);ub_S(:);ub_TG_LFCup(:);ub_TG_LFCdown(:);ub_P_PV(:);ub_P_WT(:);ub_P_EV(:);ub_E_EV(:)];
ub(3841)=Inf;

%% Ax<=b


%% TG lower boundary in each level
A=spalloc(Time*width_num*Gen_N,length(lb),Time*width_num*Gen_N*2);
counter=1;
low_P_box = repmat((P_bound(:,1:width_num)'),Time,1);
for loop1=1:Time
    for loop2=1:width_num
        for loop3=1:Gen_N
            tempTG_level=clearTG_level;
            tempP_TG=clearTG_level;
            if (loop2==1)
                tempTG_level((loop1-1)*width_num+loop2,loop3)= low_P_box((loop1-1)*width_num+loop2,loop3)+dPre0(loop3,1);
            else
                tempTG_level((loop1-1)*width_num+loop2,loop3)= low_P_box((loop1-1)*width_num+loop2,loop3);
            end
            tempP_TG((loop1-1)*width_num+loop2,loop3)=-1;
            addrow=[tempTG_level(:);tempP_TG(:);clearS(:);clearTG_LFCup(:);clearTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)]';
            A(counter,:)=sparse(addrow);
            conuter=counter+1;
        end
    end
end

b=zeros(Time*width_num*Gen_N,1);

%% TG higher boundary in each level
tempA = spalloc(Time*width_num*Gen_N,length(lb),Time*width_num*2*Gen_N);
counter = 1;
high_P_box = repmat((P_bound(:,2:width_num+1)'),Time,1);
for loop1 = 1:Time
    for loop2 = 1:width_num
        for loop3 = 1:Gen_N
            tempTG_level = clearTG_level;
            tempP_TG = clearP_TG;
            if loop2==width_num
                tempTG_level((loop1-1)*width_num+loop2,loop3) = -( high_P_box((loop1-1)*width_num+loop2,loop3) - dPre0(loop3,1) );
            else
                tempTG_level((loop1-1)*width_num+loop2,loop3) = -high_P_box((loop1-1)*width_num+loop2,loop3);
            end
            tempP_TG((loop1-1)*width_num+loop2,loop3) = 1;
            addrow=[tempTG_level(:);tempP_TG(:);clearS(:);clearTG_LFCup(:);clearTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)]';
            tempA(counter,:) = sparse(addrow);
            counter = counter + 1;
        end
    end
end
tempb = zeros(Time*width_num*Gen_N,1);
A = [A;tempA];
b = [b;tempb];


%% TG LEC(up)        
tempA = spalloc(Time*Gen_N,length(lb),Time*Gen_N*(width_num*2+1));
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        tempTG_level = clearTG_level;
        tempP_TG = clearP_TG;
        tempTG_LFCup = clearTG_LFCup;
        tempTG_level( (loop1-1)*width_num + 1: loop1*width_num ,loop2) = -P_TGmax(loop2,1);
        tempP_TG( (loop1-1)*width_num + 1: loop1*width_num ,loop2) = 1;
        tempTG_LFCup(loop1,loop2) = 1;
        addrow=[tempTG_level(:);tempP_TG(:);clearS(:);tempTG_LFCup(:);clearTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)]';
        tempA(counter,:) = sparse(addrow);
        counter = counter + 1;
    end
end
tempb = zeros(Time*Gen_N,1);
A = [A;tempA];
b = [b;tempb];

%% in 10%
tempA = spalloc(Time*Gen_N,length(lb),Time*Gen_N*(width_num+1));
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        tempTG_level = clearTG_level;
        tempTG_LFCup = clearTG_LFCup;
        tempTG_level( (loop1-1)*width_num + 1: loop1*width_num ,loop2) = - (P_TGmax(loop2,1).*TG_LFC);
        tempDG_LFCup(loop1,loop2) = 1;
        addrow=[tempTG_level(:);clearP_TG(:);clearS(:);tempTG_LFCup(:);clearTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)]';
        tempA(counter,:) = sparse(addrow);
        counter = counter + 1;
    end
end
tempb = zeros(Time*Gen_N,1);
A = [A;tempA];
b = [b;tempb];
%20201211銇亾銇撱伨銇?

%% TG LEC(down)
tempA = spalloc(Time*Gen_N,length(lb),Time*Gen_N*(width_num*2+1));
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        tempTG_level = clearTG_level;
        tempP_TG = clearP_TG;
        tempTG_LFCdown = clearTG_LFCdown;
        tempTG_level( (loop1-1)*width_num + 1: loop1*width_num ,loop2) = P_TGmin(loop2,1);
        tempP_TG( (loop1-1)*width_num + 1: loop1*width_num ,loop2) = -1;
        tempTG_LFCdown(loop1,loop2) = 1;
        addrow=[tempTG_level(:);clearP_TG(:);clearS(:);clearTG_LFCup(:);tempTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)]';
        tempA(counter,:) = sparse(addrow);
        counter = counter + 1;
    end
end
tempb = zeros(Time*Gen_N,1);
A = [A;tempA];
b = [b;tempb];


%% in 10%?
tempA = spalloc(Time*Gen_N,length(lb),Time*Gen_N*(width_num+1));
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        tempTG_level = clearTG_level;
        tempTG_LFCdown = clearTG_LFCdown;
        tempTG_level( (loop1-1)*width_num + 1: loop1*width_num ,loop2) = - (P_TGmax(loop2,1).*TG_LFC);
        tempTG_LFCdown(loop1,loop2) = 1;
        addrow=[tempTG_level(:);clearP_TG(:);clearS(:);clearTG_LFCup(:);tempTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)]';
        tempA(counter,:) = sparse(addrow);
        counter = counter + 1;
    end
end
tempb = zeros(Time*Gen_N,1);
A = [A;tempA];
b = [b;tempb];


%% TG start up
tempA = spalloc((Time)*Gen_N,length(lb),(width_num*2+1)*Gen_N*Time);
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        tempTG_level = clearTG_level;
        tempS = clearS;
        if loop1 == 1
            %tempDG_level((Time)*width_num-(width_num-1):(Time)*width_num,loop2) = -1;%绲傜鏅傚埢銇ㄣ伄璧峰嫊鍙版暟銇樊銇犮亼璧峰嫊璨荤敤銈掕冻銇?
            tempTG_level((loop1)*width_num-(width_num-1):(loop1)*width_num,loop2) = 1;%闁嬪鏅傚埢銇捣鍕曘仐銇︺亜銈嬪垎銈掕拷鍔?()
        else
            tempTG_level((loop1-1)*width_num-(width_num-1):(loop1-1)*width_num,loop2) = -1;
            tempTG_level((loop1)*width_num-(width_num-1):(loop1)*width_num,loop2) = 1;
        end
        tempS(loop1,loop2) = -1;
        addrow=[tempTG_level(:);clearP_TG(:);tempS(:);clearTG_LFCup(:);clearTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)]';
        tempA(counter,:) = sparse(addrow);
        counter = counter + 1;
    end
end
tempb = zeros((Time)*Gen_N,1);
A = [A;tempA];
b = [b;tempb];

on_Time = 2;
tempA = spalloc((Time)*Gen_N*on_Time,length(lb),(width_num*2+1)*Gen_N*Time);
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        for loop3 = 1:on_Time
            tempTG_level = clearTG_level;
            tempS = clearS;
            tempS(loop1,loop2) = 1;
            tempTG_level( ( (loop1+loop3-1) - 48*( (loop1+loop3-1)>48 ) )*width_num-(width_num-1):( (loop1+loop3-1) - 48*( (loop1+loop3-1)>48 ) )*width_num,loop2) = -1;
            addrow=[tempTG_level(:);clearP_TG(:);tempS(:);clearTG_LFCup(:);clearTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)]';
            tempA(counter,:) = sparse(addrow);
            counter = counter + 1;
        end
    end
end
tempb = zeros((Time)*Gen_N*on_Time,1);
A = [A;tempA];
b = [b;tempb];

%% Aeq*x=beq

%% EV
%equation between E and P
Aeq=spalloc(Time-1,length(lb),(Time-1)*3);
counter=1;
for loop1=1:Time-1
    tempE_EV=clearE_EV;
    tempP_EV=clearP_EV;
    tempE_EV(loop1+1)=1;
    tempE_EV(loop1)=-1;
    tempP_EV(loop1)=-0.5;
    addrow=[clearTG_level(:);clearP_TG(:);clearS(:);clearTG_LFCup(:);clearTG_LFCdown(:);clearP_PV(:);clearP_WT(:);tempP_EV(:);tempE_EV(:)]';
    Aeq(counter,:)=sparse(addrow);
    counter=counter+1;
end
beq=[];
for loop1=1:Time-1
    tempbeq=0.75*E_oneEV*(Ncon(loop1+1)-Ncon(loop1))*agg_num*agg_EVnum;
    beq=[beq;tempbeq];
end


%% demand and supply
tempAeq=spalloc(Time,length(lb),Time*(Gen_N*width_num+3));
counter=1;
for loop1=1:Time
    tempP_TG=clearP_TG;
    tempP_PV=clearP_PV;
    tempP_WT=clearP_WT;
    tempP_EV=clearP_EV;
    tempP_TG(1+(loop1-1)*width_num:loop1*width_num,:)=1;
    tempP_PV(loop1)=1;
    tempP_WT(loop1)=1;
    tempP_EV(loop1)=-1;
    adrrow=[clearTG_level(:);tempP_TG(:);clearS(:);clearTG_LFCup(:);clearTG_LFCdown(:);tempP_PV(:);tempP_WT(:);tempP_EV(:);clearE_EV(:)]';
    tempAeq(counter,:)=sparse(addrow);
    counter=counter+1;
end
tempbeq=P_load-P_nuclear-P_hydro;
Aeq=[Aeq;tempAeq];
beq=[beq;tempbeq];


%% object function
tempTG_level = clearTG_level;
tempP_TG = clearP_TG;
tempS = clearS;
for loop2 = 1:Gen_N
    tempTG_level(:,loop2) = repmat(beta(loop2,:)',Time,1);
    tempP_TG(:,loop2) = repmat(alpha(loop2,:)',Time,1);
    tempS(:,loop2) = repmat(S0(loop2,1),Time,1);
end
f =  [tempTG_level(:);tempP_TG(:);tempS(:);clearTG_LFCup(:);clearTG_LFCdown(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)];



 %% 澶夋暟銈裤偆銉楁寚瀹?
intcon=[];
for i = 1:length(lb)
    if i <= size(clearTG_level(:),1)
        intcon=[intcon,i];
    elseif i <= size([clearTG_level(:);clearP_TG(:)],1)
        continue;
        %consistent
    elseif i <= size([clearTG_level(:);clearP_TG(:);clearS(:)],1)
        intcon=[intcon,i];;%binary
    elseif i <= size([clearTG_level(:);clearP_TG(:);clearS(:);clearTG_LFCup(:);clearTG_LFCdown(:)],1)
        continue;
        %consistent
    elseif i <= size([clearTG_level(:);clearP_TG(:);clearS(:);clearTG_LFCup(:);clearTG_LFCdown(:);tempP_PV(:);tempP_WT(:);tempP_EV(:);clearE_EV(:)],1)
        continue;
    else
        intcon=[intcon,i];
    end
end


%% MILP

x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
%f =  [tempTG_level(:);tempP_TG(:);tempS(:);clearTG_LFCup(:);clearTG_LFCdown(:);tempP_PV(:);tempP_WT(:);tempP_EV(:);clearE_EV(:)]


count = 1;
TG_level = [];
for loop1 = 1:Gen_N
    for loop2 = 1:Time
        for loop3 = 1:width_num
            TG_level(loop3,loop2,loop1) = x(count);
            count = count + 1;
        end
    end
end
P_TG = [];
for loop1 = 1:Gen_N
    for loop2 = 1:Time
        for loop3 = 1:width_num
            P_TG(loop3,loop2,loop1) = x(count);
            count = count + 1;
        end
    end
end
S = [];
for loop1 = 1:Gen_N
    for loop2 = 1:Time
        S(loop1,loop2) = x(count);
        count = count + 1;
    end
end
TG_LFCup = [];
for loop1 = 1:Gen_N
    for loop2 = 1:Time
        TG_LFCup(loop1,loop2) = x(count);
        count = count + 1;
    end
end
TG_LFCdown = [];
for loop1 = 1:Gen_N
    for loop2 = 1:Time
        TG_LFCdown(loop1,loop2) = x(count);
        count = count + 1;
    end
end
P_PV = [];
for loop1 = 1:Time
    P_PV(loop1) = x(count);
    count = count + 1;
end
P_WT = [];
for loop1 = 1:Time
    P_WT(loop1) = x(count);
    count = count + 1;
end
P_EV=[];
for loop1=1:Time
    P_EV(loop1)=x(count);
    count=count+1;
end
E_EV=[];
for loop1=1:Time
    E_EV(loop1)=x(count);
    count=count+1;
end

%% TG
for loop1 = 1:Gen_N
    P_fire(loop1,:) = sum(P_TG(:,:,loop1));
end
TOTAL_P_TG = sum(P_fire);

%20201213銇亾銇撱伨銇с仩
%% 
Load=P_load'
surplus=TOTAL_P_TG+P_PV+P_WT+P_nuclear+P_hydro-Load;

%% 
timeline=0:0.5:24;
Hydro=repmat(P_hydro,1,48);
Nuclear=repmat(P_nuclear,1,48);
Thermal=TOTAL_P_TG;
PV=P_PV;
WT=P_WT;
Surplus=surplus;
EV=P_EV;

UC_result=[time1;Hydro;Nuclear;Thermal;PV;WT;Surplus;Load-EV];
  