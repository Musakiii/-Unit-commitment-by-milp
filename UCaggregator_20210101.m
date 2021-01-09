% allocate P to EVs according to result of UC
%cooperate with milp0101
%{
20210101
for 全国大会　
rated P_PV=1000MW

20201230
for 全国大会

20201222
test with milp1222
EVSOC seems better than milp1220 and the front part
but it still will be more than 85 when the time near to 24h
%}
%{
20201220
a new algorithm to complete EVSOC
because of the difference between EV patterns, build a P table to each EV
and each EV's SOC can be caculated by this P table.
the point is how to build the P table
according to newEVstate:
0 in control state, so need to get data from P_agg allocation loop
1 in forcibly charging/discharging state, P is decided by the 1st 1 start
time and the SOC at that time
2 in travel state, it's easy
3 in none state, equal to SOC_t-1
%}
%{
20201219
update EVSOC
load newEVstate (with state 3)
%}
%{
20201218
operate successfully (allocate P)(Ncon=size(EVlist,2))
%}
%{
20201216
('EV_result.mat','E_oneEV','P_EVrate','EVUCrate','P_oneEV','agg_EVnum','agg_num')
('UC_result.mat','timeline','fire','hydro','nuclear','P_PV','P_WT','P_EV','rest_demand','surplus','P_demand','Ncon');
sort controallable EVs by charging or discharging in a power of 4kW or -4kW(UC rate power)
%}






clear;
%import bioma.data.*;
load('UC_result_20210101.mat');
load('EV_result_20210101.mat');
load('EV_make_data_20201217.mat');
%EV state
%0:controllable,1:forcibly,2:traveling,3:none state



% MW 2 kW
P_EVrate=P_EVrate*1000;
P_agg=P_EV'*1000/agg_num;
P_oneEV=P_EVrate*EVUCrate;
E_oneEV=E_oneEV*1000;
P_SOC=P_oneEV/B_EV;
%
SOC_incontrol=SOCtarget;




%20201219
%update EVSOC
% newEVSOC=zeros(EVnum,86400);
% newEVSOC(:,1)=SOC_initial;
% [newEVSOC,initial_duration]=fun_forciblycharge(newEVSOC,EVnum,ones(EVnum,1),P_EVmax,B_EV,SOCtarget);
% timepoint1=[];
% for loop1=1:EVnum
%     for loop2=1:86400
%         if newEVSOC(loop1,loop2)==0
%             timepoint1=[timepoint1;loop2]; 
%             break;
%         end
%     end  
% end
% for loop1=1:EVnum
%     EVSOC(loop1,timepoint1(loop1):1800*(fix((timepoint1(loop1)-1)/1800)+1)+1)=newEVSOC(loop1,timepoint1(loop1)-1);
% end
% for loop1=1:EVnum
%     EVSOC(loop1,fix((t_travel1_second(loop1)-3600-1)/1800)*1800+1:t_travel1_second(loop1))= linspace(EVSOC(loop1,fix((t_travel1_second(loop1)-3600-1)/1800)*1800),0.85,t_travel1_second(loop1)-(fix((t_travel1_second(loop1)-3600-1)/1800)*1800+1)+1); 
%     EVSOC(loop1,t_travel1_second(loop1):t_travel1_second(loop1)+travel1_duration(loop1))=linspace(EVSOC(loop1,t_travel1_second(loop1)),EVSOC(loop1,t_travel1_second(loop1)+travel1_duration(loop1)),travel1_duration(loop1)+1);
% end  


%20201220
%build a P_table
P_table=zeros(EVnum,86400);
%2 travel state
P_travel=-1*velocity/cost;
for loop1=1:EVnum
    P_table(loop1,t_travel1_second(loop1):t_end1(loop1))=P_travel;
    P_table(loop1,t_travel2_second(loop1):t_end2(loop1))=P_travel; 
end
%3 none state
for loop1=1:EVnum
    for loop2=1:86400
        if newEVstate(loop1,loop2)==3
            P_table(loop1,loop2)=0; 
        end
    end
end

test=[];
upgroup=[];
downgroup=[];
%0 controllable state
for t=1:deltaT*60*60:24*60*60
    %update data
    controlinEVlist=[];
    if t<=1800
        EVlist=[];
        continue;
    else
        for loop1=1:EVnum
            if (EVstate(loop1,t)==0) & (EVstate(loop1,t-deltaT*60*60)~=0) & (EVstate(loop1,min(t+deltaT*3600-1,86400))==0)
                controlinEVlist=[controlinEVlist,[loop1;0.75]]; 
            end
        end
    end
    EVlist=fun_updateEVlist(EVlist,controlinEVlist);
    test=[test;size(EVlist,2)];
    
    EVlist=fun_sort(EVlist);
    
    %allocate signal to EVs
    signal=P_agg(fix(t/1800)+1);
    [up,down]=fun_allocate2EV(signal,EVlist);
    upgroup=[upgroup,up];
    downgroup=[downgroup,down];
    %20201220
    %update P_table
    count=1;
    for loop1=1:size(EVlist,2)
        if count<=down
            P_table(EVlist(1,loop1),t:t+1800-1)=-1*P_oneEV;
            count=count+1;
        else
            P_table(EVlist(1,loop1),t:t+1800-1)=P_oneEV;
            count=count+1;
        end
        
    end 
   
    
    
    %after deltaT
    %update EV data
    if down>0
        EVlist(2,1:down)=EVlist(2,1:down)-P_SOC*deltaT;
    end
    if down<size(EVlist,2)
        EVlist(2,down+1:end)=EVlist(2,down+1:end)+P_SOC*deltaT;
    end
    
    %sort the EVlist
    %{
    EVlist = DataMatrix(EVlist);
    EVlist=sortcols(EVlist,2,'descend');
    %}
    EVlist=fun_sort(EVlist);
    
    

    
    
    
    %update EVSOC
    count=0;
    for loop1=EVlist(1,:)
        count=count+1;
        EVSOC(loop1,t:min(t+deltaT*3600,24*3600))=linspace(EVSOC(loop1,t),EVlist(2,count),min(t+deltaT*3600,24*3600)-t+1); 
    end
    
    %delete uncontrollable EV at next t in EVlist
    
    colnum=0;
    deletelist=[];
    for loop1=EVlist(1,:)
        colnum=colnum+1;
        if (t+deltaT*3600<=86400) & (EVstate(loop1,min(86400,t+deltaT*3600*2-1))~=0)
            deletelist=[deletelist,colnum];
        end
    end
    EVlist(:,deletelist)=[];
end

%20201219
%update EVSOC
%20201220
%it's a too bad algorithm
%{
for loop1=1:EVnum
    for loop2=1:86400
        if newEVstate(loop1,loop2)==3
            EVSOC(loop1,loop2)=EVSOC(loop1,loop2-1); 
        end
    end
end


for loop1=1:EVnum
    EVSOC(loop1,fix((t_travel1_second(loop1)-3600-1)/1800)*1800+1:t_travel1_second(loop1))= linspace(EVSOC(loop1,fix((t_travel1_second(loop1)-3600-1)/1800)*1800),0.85,t_travel1_second(loop1)-(fix((t_travel1_second(loop1)-3600-1)/1800)*1800+1)+1);
    EVSOC(loop1,fix((t_travel2_second(loop1)-3600-1)/1800)*1800+1:t_travel2_second(loop1))= linspace(EVSOC(loop1,fix((t_travel2_second(loop1)-3600-1)/1800)*1800),0.85,t_travel2_second(loop1)-(fix((t_travel2_second(loop1)-3600-1)/1800)*1800+1)+1); 
    EVSOC(loop1,t_travel1_second(loop1):t_travel1_second(loop1)+travel1_duration(loop1))=linspace(EVSOC(loop1,t_travel1_second(loop1)),EVSOC(loop1,t_travel1_second(loop1)+travel1_duration(loop1)),travel1_duration(loop1)+1);
    EVSOC(loop1,t_travel2_second(loop1):t_travel2_second(loop1)+travel2_duration(loop1))=linspace(EVSOC(loop1,t_travel2_second(loop1)),EVSOC(loop1,t_travel2_second(loop1)+travel2_duration(loop1)),travel2_duration(loop1)+1);
end

%caculate the 3rd and 5th last state3 time
complete=zeros(EVnum,6);
for loop1=1:EVnum
    count=1;
    for loop2=1:86400
        if (newEVstate(loop1,loop2)==3) & (newEVstate(loop1,loop2-1)==1)
            complete(loop1,count)=loop2;
            count=count+1;
        end
        if (newEVstate(loop1,loop2)==3) & (newEVstate(loop1,min(loop2+1,86400))==0)
            complete(loop1,count)=loop2;
            count=count+1;
        end        
    end
end


for loop1=1:EVnum
    if complete(loop1,5)~=complete(loop1,6)
        EVSOC(loop1,fix((complete(loop1,3)-1)/1800)*1800:complete(loop1,4))=linspace(EVSOC(loop1,fix((complete(loop1,3)-1)/1800)*1800),EVSOC(loop1,complete(loop1,4)),complete(loop1,4)-fix((complete(loop1,3)-1)/1800)*1800+1);
        EVSOC(loop1,fix((complete(loop1,5)-1)/1800)*1800:complete(loop1,6))=linspace(EVSOC(loop1,fix((complete(loop1,5)-1)/1800)*1800),EVSOC(loop1,complete(loop1,6)),complete(loop1,6)-fix((complete(loop1,5)-1)/1800)*1800+1);
    else
        EVSOC(loop1,max(fix((complete(loop1,1)-1)/1800)*1800,1):complete(loop1,2))=linspace(EVSOC(loop1,max(fix((complete(loop1,1)-1)/1800)*1800,1)),EVSOC(loop1,complete(loop1,2)),complete(loop1,2)-max(fix((complete(loop1,1)-1)/1800)*1800,1)+1);
        EVSOC(loop1,fix((complete(loop1,3)-1)/1800)*1800:complete(loop1,4))=linspace(EVSOC(loop1,fix((complete(loop1,3)-1)/1800)*1800),EVSOC(loop1,complete(loop1,4)),complete(loop1,4)-fix((complete(loop1,3)-1)/1800)*1800+1);
    end
end
%}

%20201220
%all EV are in 1 state and has been construct in now EVSOC
%in 3 state SOC don't change
for loop1=1:EVnum
    for loop2=1:86400
        if newEVstate(loop1,loop2)==3
            EVSOC(loop1,loop2)=EVSOC(loop1,loop2-1); 
        end
    end   
end
%in 0 state
for loop1=1:EVnum
    for loop2=2:86400
        if newEVstate(loop1,loop2)==0
            EVSOC(loop1,loop2)=EVSOC(loop1,loop2-1)+P_table(loop1,loop2-1)/E_oneEV/3600;
        end
    end   
end

%find each EV's now SOC=0's start time and end time
for loop1=1:EVnum
    count=1;
    for loop2=1:86400-1
        if (EVSOC(loop1,loop2)~=0) & (EVSOC(loop1,loop2+1)==0)
            find(loop1,count)=loop2;
            count=count+1;
        end
        if (EVSOC(loop1,loop2)==0) & (EVSOC(loop1,loop2+1)~=0)
            find(loop1,count)=loop2+1;
            count=count+1;
        end
    end
    
end

%complete the last part(after contorl's forcibly part)
for loop1=1:EVnum
    EVSOC(loop1,find(loop1,1):find(loop1,2))=linspace(EVSOC(loop1,find(loop1,1)),EVSOC(loop1,find(loop1,2)),find(loop1,2)-find(loop1,1)+1);
    EVSOC(loop1,find(loop1,3):find(loop1,4))=linspace(EVSOC(loop1,find(loop1,3)),EVSOC(loop1,find(loop1,4)),find(loop1,4)-find(loop1,3)+1); 
end


%check the highest and lowest SOC in 0 state
highest=zeros(EVnum,1);
lowest=ones(EVnum,1);
for loop1=1:EVnum
    for loop2=1:86400
        if newEVstate(loop1,loop2)==0
            if EVSOC(loop1,loop2)>highest(loop1)
               highest(loop1)=EVSOC(loop1,loop2); 
            end
            if EVSOC(loop1,loop2)<lowest(loop1)
                lowest(loop1)=EVSOC(loop1,loop2); 
            end
        end
    end
end

 save('EVSOC20210101ver2.mat','EVSOC');


