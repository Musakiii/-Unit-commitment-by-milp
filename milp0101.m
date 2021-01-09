%{
20210101 Release notes
according to Professor Yokoyama's comment
rated PV from 630 to 1000
the end of data name is 20210101

20201230 Release notes
全国大会用データの出力
the end of data name is 20201230

20201222 Release notes
introduce up&down group to restrict a EV charge or discharge constantly in 
more than 3 dT (one deltaT EV could charge 4*0.5=2kWh equal to 2/60 SOC, 
and if constantly charge more than 3 dT, it will up to 10% SOC which is 
equal to EVSOC constraints 75-85 or 75-65)
%}

clear;
load('PVdata_20201210.mat','PVdata');
load('load_weekday_20201211.mat','load_weekday');
load('Ncon_20201217','Ncon');
load('wind_data20201229.mat');
load('P_load20201229.mat');%load the load data for 全国大会
%% constant
EVUCrate=4/6;
P_EVrate=0.006;
E_oneEV=0.06;
P_oneEV=P_EVrate*EVUCrate;
agg_EVnum=1500;
Time=48;
width_num=4;
agg_num=40;
P_WTrate=200;

PVdata=PVdata*100*1.5873;%20210101 make rated P_PV=1000
load_weekday(2,:)=load_weekday(2,:)*0.22;
P_nuclear=800;
P_hydro=200;
%P_load=load_weekday(2,0.25*36000:0.5*36000:864001)';
P_load=testdata;%20201230 input the 全国大会 load curve,this curve's max is 3500 that 300 more than 3200

%Ncon(1)=Ncon(2);

thermal_data=[
60 80 304 300 620 591 350;
30 40 152 150 310 295.5 175;
0.0646 0.025 0.0533 0.0224 0.0067 0.0081 0.0032;
16.8899 43.5 9.2374 19.7128 10.2202 20.2909 10.102;
49.7108 200 164.3492 287.4982 207.1786 378.5918 350.6398;
181.6416 13.356 111.4208 1511.8992 1777.1544 2070.18 8331.9264
];


dT=0.5;
[alpha,beta,P_bound] = func_linearized(thermal_data,width_num,dT);


P_PVmax=PVdata(15:30:1440)';
P_WTmax=ones(48,1)*P_WTrate;
P_TGmax=thermal_data(1,:)';
P_PVmin=zeros(48,1);
P_WTmin=zeros(48,1);
P_TGmin=thermal_data(2,:)';
S0=thermal_data(6,:)';%thermal startup cost
dPre0=P_TGmax*0.05;%LFC 5%
Gen_N=size(thermal_data,2);


clearDG_level = zeros(Time*width_num,Gen_N);
clearP_DG     = zeros(Time*width_num,Gen_N);
clearS        = zeros(Time,Gen_N);
clearP_PV = zeros(Time,1);
clearP_WT = zeros(Time,1);
clearP_EV=zeros(Time,1);
clearE_EV=zeros(Time,1);
%20201222 number of group_plus+group_minus=Ncon
clearEVgroup_plus=zeros(Time,1);
clearEVgroup_minus=zeros(Time,1);




lb_DG_level = zeros(Time*width_num,Gen_N);
lb_P_DG     = zeros(Time*width_num,Gen_N);
lb_S        = zeros(Time,Gen_N);
lb_P_PV = zeros(Time,1);
lb_P_WT = zeros(Time,1);
lb_P_EV=-1*P_oneEV*Ncon*agg_num;
lb_E_EV=0.65*E_oneEV*Ncon*agg_num;
lb_EVgroup_plus=zeros(Time,1);
lb_EVgroup_minus=zeros(Time,1);
% lb_EVgroup_plus=-inf(Time,1);
% lb_EVgroup_minus=-inf(Time,1);
lb=[lb_DG_level(:);lb_P_DG(:);lb_S(:);lb_P_PV(:);lb_P_WT(:);lb_P_EV(:);lb_E_EV(:);lb_EVgroup_plus(:);lb_EVgroup_minus(:)];


ub_DG_level = ones(Time*width_num,Gen_N);
ub_P_DG     = repmat(P_TGmax',Time*width_num,1);
ub_S        = ones(Time,Gen_N);
ub_P_PV = P_PVmax;
ub_P_WT = ones(Time,1)*P_WTrate;
ub_P_EV=P_oneEV*Ncon*agg_num;
ub_E_EV=0.85*E_oneEV*Ncon*agg_num;
ub_EVgroup_plus=Ncon*agg_num;
ub_EVgroup_minus=Ncon*agg_num;
% ub_EVgroup_plus=inf(Time,1);
% ub_EVgroup_minus=inf(Time,1);
ub=[ub_DG_level(:);ub_P_DG(:);ub_S(:);ub_P_PV(:);ub_P_WT(:);ub_P_EV(:);ub_E_EV(:);ub_EVgroup_plus(:);ub_EVgroup_minus(:)];




%% DGの運転レベルごとの最低出力制約
A = spalloc(Time*width_num*Gen_N,length(lb),Time*width_num*2*Gen_N);
counter = 1;
low_P_box = repmat((P_bound(:,1:width_num)'),Time,1);
for loop1 = 1:Time
    for loop2 = 1:width_num
        for loop3 = 1:Gen_N
            tempDG_level = clearDG_level;
            tempP_DG = clearP_DG;
            if (loop2==1)
%                 tempDG_level((loop1-1)*width_num+loop2,loop3) = low_P_box((loop1-1)*width_num+loop2,loop3) + dPre0(loop3,1);
                tempDG_level((loop1-1)*width_num+loop2,loop3) = 0;
            else
%                 tempDG_level((loop1-1)*width_num+loop2,loop3) = low_P_box((loop1-1)*width_num+loop2,loop3);
                tempDG_level((loop1-1)*width_num+loop2,loop3) = 0;
            end
            tempP_DG((loop1-1)*width_num+loop2,loop3) = -1;
            addrow = [tempDG_level(:);tempP_DG(:);clearS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);clearEVgroup_plus(:);clearEVgroup_minus(:)]';
            A(counter,:) = sparse(addrow);
            counter = counter + 1;
        end
    end
end
b = zeros(Time*width_num*Gen_N,1);

%% DGの運転レベルごとの最高出力
tempA = spalloc(Time*width_num*Gen_N,length(lb),Time*width_num*2*Gen_N);
counter = 1;
high_P_box = repmat((P_bound(:,2:width_num+1)'),Time,1);
for loop1 = 1:Time
    for loop2 = 1:width_num
        for loop3 = 1:Gen_N
            tempDG_level = clearDG_level;
            tempP_DG = clearP_DG;
            if loop2==width_num
                tempDG_level((loop1-1)*width_num+loop2,loop3) = -( high_P_box((loop1-1)*width_num+loop2,loop3) - dPre0(loop3,1) );
            else
                tempDG_level((loop1-1)*width_num+loop2,loop3) = -high_P_box((loop1-1)*width_num+loop2,loop3);
            end
            tempP_DG((loop1-1)*width_num+loop2,loop3) = 1;
            addrow = [tempDG_level(:);tempP_DG(:);clearS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);clearEVgroup_plus(:);clearEVgroup_minus(:)]';
            tempA(counter,:) = sparse(addrow);
            counter = counter + 1;
        end
    end
end
tempb = zeros(Time*width_num*Gen_N,1);
A = [A;tempA];
b = [b;tempb];

% level
tempA = spalloc(Time*Gen_N,length(lb),Time*Gen_N*width_num);
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        tempDG_level = clearDG_level;
        tempDG_level(loop1*width_num-(width_num-1):loop1*width_num,loop2) = -1;
        addrow = [tempDG_level(:);clearP_DG(:);clearS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);clearEVgroup_plus(:);clearEVgroup_minus(:)]';
        tempA(counter,:) = sparse(addrow);
        counter = counter + 1;
    end
end
tempb = zeros(Time*Gen_N,1);
A = [A;tempA];
b = [b;tempb];
%sum(z)<=1
tempA = spalloc(Time*Gen_N,length(lb),Time*Gen_N*width_num);
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        tempDG_level = clearDG_level;
        tempDG_level(loop1*width_num-(width_num-1):loop1*width_num,loop2) = 1;
        addrow = [tempDG_level(:);clearP_DG(:);clearS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);clearEVgroup_plus(:);clearEVgroup_minus(:)]';
        tempA(counter,:) = sparse(addrow);
        counter = counter + 1;
    end
end
tempb = ones(Time*Gen_N,1);
A = [A;tempA];
b = [b;tempb];


%% DG起動制約式（起動費に関するところ）
% s >= sum(z(t)) - sum(z(t-1))で0 or 1
% 勿論0のほうが安いから右辺が0ならsも当然0をせんたくする
tempA = spalloc((Time)*Gen_N,length(lb),(width_num*2+1)*Gen_N*Time);
counter = 1;
for loop1 = 1:Time
    for loop2 = 1:Gen_N
        tempDG_level = clearDG_level;
        tempS = clearS;
        if loop1 == 1
            %tempDG_level((Time)*width_num-(width_num-1):(Time)*width_num,loop2) = -1;%終端時刻との起動台数の差だけ起動費用を足す
            tempDG_level((loop1)*width_num-(width_num-1):(loop1)*width_num,loop2) = 1;%開始時刻に起動している分を追加()
        else
            tempDG_level((loop1-1)*width_num-(width_num-1):(loop1-1)*width_num,loop2) = -1;
            tempDG_level((loop1)*width_num-(width_num-1):(loop1)*width_num,loop2) = 1;
        end
        tempS(loop1,loop2) = -1;
        addrow = [tempDG_level(:);clearP_DG(:);tempS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);clearEVgroup_plus(:);clearEVgroup_minus(:)]';
        tempA(counter,:) = sparse(addrow);
        counter = counter + 1;
    end
end
tempb = zeros((Time)*Gen_N,1);
A = [A;tempA];
b = [b;tempb];


%2020122 EVgroup inequality constraints
%minus2+minus3+minus4>=plus1
tempA = spalloc(Time-3,length(lb),4*(Time-3));
counter = 1;
for loop1 = 1:Time-3
    tempEVgroup_plus=clearEVgroup_plus;
    tempEVgroup_minus=clearEVgroup_minus;
    tempEVgroup_plus(loop1)=1;
    tempEVgroup_minus(loop1+1:loop1+3)=-1;
    addrow = [tempDG_level(:);clearP_DG(:);tempS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);tempEVgroup_plus(:);tempEVgroup_minus(:)]';
    tempA(counter,:) = sparse(addrow);
    counter = counter + 1;
end
tempb = zeros(Time-3,1);
A = [A;tempA];
b = [b;tempb];


%plus2+plus3+plus4>=minus1
tempA = spalloc(Time-3,length(lb),4*(Time-3));
counter = 1;
for loop1 = 1:Time-3
    tempEVgroup_plus=clearEVgroup_plus;
    tempEVgroup_minus=clearEVgroup_minus;
    tempEVgroup_minus(loop1)=1;
    tempEVgroup_plus(loop1+1:loop1+3)=-1;
    addrow = [tempDG_level(:);clearP_DG(:);tempS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);tempEVgroup_plus(:);tempEVgroup_minus(:)]';
    tempA(counter,:) = sparse(addrow);
    counter = counter + 1;
end
tempb = zeros(Time-3,1);
A = [A;tempA];
b = [b;tempb];


%% 需給一致制約
Aeq = spalloc(Time,length(lb),Time*(width_num*Gen_N+3));
counter = 1;
beq = zeros(48,1);
for loop1 = 1:Time
    tempP_DG = clearP_DG;
    tempP_PV = clearP_PV;
    tempP_WT = clearP_WT;
    tempP_EV = clearP_EV;
    tempP_DG(1+(loop1-1)*width_num:loop1*width_num,:) = 1;
    tempP_PV(loop1,1) = 1;
    tempP_WT(loop1,1) = 1;
    tempP_EV(loop1,1) = -1; 
    addrow = [clearDG_level(:);tempP_DG(:);clearS(:);tempP_PV(:);tempP_WT(:);tempP_EV(:);clearE_EV(:);clearEVgroup_plus(:);clearEVgroup_minus(:)]';
    Aeq(counter,:) = sparse(addrow);
    counter = counter + 1;
    beq(loop1,1) = P_load(loop1,1) - P_hydro-P_nuclear;% + EVout_P;
end

tempAeq=spalloc(Time-1,length(lb),(Time-1)*3);
counter=1;
for loop1=1:Time-1
    tempE_EV=clearE_EV;
    tempP_EV=clearP_EV;
    tempE_EV(loop1+1)=1;
    tempE_EV(loop1)=-1;
    tempP_EV(loop1)=-0.5;
    addrow=[clearDG_level(:);clearP_DG(:);clearS(:);clearP_PV(:);clearP_WT(:);tempP_EV(:);tempE_EV(:);clearEVgroup_plus(:);clearEVgroup_minus(:)]';
    tempAeq(counter,:)=sparse(addrow);
    counter=counter+1;
end
Aeq=[Aeq;tempAeq];
for loop1=1:Time-1
    tempbeq=0.75*E_oneEV*(Ncon(loop1+1)-Ncon(loop1))*agg_num;
    beq=[beq;tempbeq];
end

%2020122 EVgroup equality constraints
%plus+minus=Ncon
tempAeq=spalloc(Time,length(lb),Time*2);
counter=1;
for loop1=1:Time
    tempEVgroup_plus=clearEVgroup_plus;
    tempEVgroup_minus=clearEVgroup_minus;
    tempEVgroup_plus(loop1)=1;
    tempEVgroup_minus(loop1)=1;
    addrow=[clearDG_level(:);clearP_DG(:);clearS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);tempEVgroup_plus(:);tempEVgroup_minus(:)]';
    tempAeq(counter,:)=sparse(addrow);
    counter=counter+1;
end
Aeq=[Aeq;tempAeq];
for loop1=1:Time
    tempbeq=Ncon(loop1)*agg_num;
    beq=[beq;tempbeq];
end
% P_oneEV*plus-P_oneEVminus=P
tempAeq=spalloc(Time,length(lb),Time*3);
counter=1;
for loop1=1:Time
    tempEVgroup_plus=clearEVgroup_plus;
    tempEVgroup_minus=clearEVgroup_minus;
    tempP_EV=clearP_EV;
    tempEVgroup_plus(loop1)=P_oneEV;
    tempEVgroup_minus(loop1)=-1*P_oneEV;
    tempP_EV(loop1)=-1;
    addrow=[clearDG_level(:);clearP_DG(:);clearS(:);clearP_PV(:);clearP_WT(:);tempP_EV(:);clearE_EV(:);tempEVgroup_plus(:);tempEVgroup_minus(:)]';
    tempAeq(counter,:)=sparse(addrow);
    counter=counter+1;
end
Aeq=[Aeq;tempAeq];
for loop1=1:Time
    tempbeq=0;
    beq=[beq;tempbeq];
end






% object function       

tempDG_level=clearDG_level;
tempP_DG = clearP_DG;
tempS=clearS;
for loop2 = 1:Gen_N
    tempDG_level(:,loop2) = repmat(beta(loop2,:)',Time,1);
    tempP_DG(:,loop2) = repmat(alpha(loop2,:)',Time,1);
    tempS(:,loop2) = repmat(S0(loop2,1),Time,1);
end
f =  [tempDG_level(:);tempP_DG(:);tempS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:);clearEVgroup_plus(:);clearEVgroup_minus(:)];


% TGub=repmat(sum(thermal_data(1,:)),Time,1);
% TGlb=repmat(0*sum(thermal_data(1,:)),Time,1);
% supplyub=TGub+P_PVmax+P_WTmax+P_nuclear+P_hydro;
% supplylb=TGlb+P_PVmin+P_WTmin+P_nuclear+P_hydro;
% 
% compare(:,1)=supplylb;
% compare(:,2)=P_load;
% compare(:,3)=supplyub;
% 
% y=0.5:0.5:24.1;
% plot(y,compare(:,1)');
% hold on;
% plot(y,compare(:,2)');
% hold on;
% plot(y,compare(:,3)');
% hold on;

intcon=[];
for i=1:length(lb)
    if i<=size([clearDG_level(:)],1)
        intcon=[intcon,i];
    elseif i<=size([clearDG_level(:);clearP_DG(:)],1)
        continue;            
    elseif i<=size([clearDG_level(:);clearP_DG(:);clearS(:)],1)
        intcon=[intcon,i];
    elseif i<=size([clearDG_level(:);clearP_DG(:);clearS(:);clearP_PV(:);clearP_WT(:);clearP_EV(:);clearE_EV(:)],1)
        continue;
    end
end


[x,fval,exitflag,output]  = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);

%% allocate variables
count = 1;
DG_level = [];
for loop1 = 1:Gen_N
    for loop2 = 1:Time
        for loop3 = 1:width_num
            DG_level(loop3,loop2,loop1) = x(count);
            count = count + 1;
        end
    end
end
P_DG = [];
for loop1 = 1:Gen_N
    for loop2 = 1:Time
        for loop3 = 1:width_num
            P_DG(loop3,loop2,loop1) = x(count);
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
P_PV=[];
for loop=1:Time
    P_PV(1,loop)=x(count);
    count = count + 1;
end
P_WT=[];
for loop=1:Time
    P_WT(1,loop)=x(count);
    count = count + 1;
end
P_EV=[];
for loop=1:Time
    P_EV(1,loop)=x(count);
    count = count + 1;
end
E_EV=[];
for loop=1:Time
    E_EV(1,loop)=x(count);
    count = count + 1;
end
EVgroup_plus=[];
for loop=1:Time
    EVgroup_plus(1,loop)=x(count);
    count = count + 1;
end
EVgroup_minus=[];
for loop=1:Time
    EVgroup_minus(1,loop)=x(count);
    count = count + 1;
end



tempfire=[];
for loop1 = 1:Gen_N
    tempfire(loop1,:) = sum(P_DG(:,:,loop1));
end
fire = sum(tempfire);

rest_demand=P_load'-P_PV-P_WT-P_hydro-P_nuclear;
surplus=fire-rest_demand;

%% 需給データ
timeline = 0:0.5:23.5;
hydro = repmat(P_hydro,1,48);
nuclear=repmat(P_nuclear,1,48);
P_demand=P_load';
P_EV=P_EV;
Ncon=Ncon';
UC_result = [timeline;fire;hydro;nuclear;P_PV;P_WT;P_EV;rest_demand;surplus;P_demand;Ncon];
% save('EV_result_20201222.mat','E_oneEV','P_EVrate','EVUCrate','P_oneEV','agg_EVnum','agg_num');
% save('UC_result_20201222.mat','timeline','fire','hydro','nuclear','P_PV','P_WT','P_EV','rest_demand','surplus','P_demand','Ncon');
save('EV_result_20210101.mat','E_oneEV','P_EVrate','EVUCrate','P_oneEV','agg_EVnum','agg_num');%20201230
save('UC_result_20210101.mat','timeline','fire','hydro','nuclear','P_PV','P_WT','P_EV','rest_demand','surplus','P_demand','Ncon');%20201230












