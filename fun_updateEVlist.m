function [newEVlist] = fun_updateEVlist(oldEVlist,controlinEVlist)
%UNTITLED10 此处显示有关此函数的摘要
%   
templist=oldEVlist;
count=0;
if size(oldEVlist,2)<1
    newEVlist=controlinEVlist;
else
    for loop1=1:size(oldEVlist,2)
        count=count+1;
        if oldEVlist(2,loop1)<0.75
            newEVlist=[templist(:,1:loop1-1),controlinEVlist,templist(:,loop1:end)];
            break;
        else
            continue; 
        end
    end
if (count==size(oldEVlist,2)) & (oldEVlist(2,count)>=0.75)
    newEVlist=[templist(:,1:end),controlinEVlist]; 
end    
    
end



end

