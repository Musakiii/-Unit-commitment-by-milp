function [newSOC] = fun_travelconsumption(EVSOC,t_start,t_end,velocity,cost,B,EVnum)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
tempSOC=EVSOC;
duration=t_end-t_start;
consumption=(t_end-t_start)*(velocity/3600)/cost/B;
for ii=1:EVnum
    tempSOC(ii,t_start(ii):t_end(ii))=linspace(EVSOC(ii,t_start(ii)),EVSOC(ii,t_start(ii))-consumption(ii),duration(ii)+1); 
end

newSOC=tempSOC;

end

