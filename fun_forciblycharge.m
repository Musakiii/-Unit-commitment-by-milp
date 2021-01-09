function [newSOC,output_duration] = fun_forciblycharge(EVSOC,EVnum,time,P,B,SOCtarget)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
tempEVSOC=EVSOC;
P_SOC=P/B/3600;
timeMAX=24*60*60;
output_duration=zeros(EVnum,1);
for ii=1:EVnum
    if tempEVSOC(ii,time(ii))<SOCtarget
        duration=fix((SOCtarget-tempEVSOC(ii,time(ii)))/P_SOC)+1;
    else
        duration=fix((tempEVSOC(ii,time(ii))-SOCtarget)/P_SOC)+1;
    end
    if duration+time(ii)<timeMAX
        tempEVSOC(ii,time(ii):time(ii)+duration)=linspace(tempEVSOC(ii,time(ii)),SOCtarget,duration+1);
        output_duration(ii)=duration;
    else
        delta=timeMAX-time(ii);
        output_duration(ii)=delta;
        tempEVSOC(ii,time(ii):timeMAX)=linspace(tempEVSOC(ii,time(ii)),tempEVSOC(ii,time(ii))+P_SOC*delta,delta+1);
    end
end

newSOC=tempEVSOC;

end

