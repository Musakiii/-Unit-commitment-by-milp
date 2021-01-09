function [newEVlist] = fun_sort(oldEVlist)
%UNTITLED 此处显示有关此函数的摘要
%   descend

% result=[];
% newlist=oldEVlist;
%     for loop1=1:size(oldEVlist,2)
%         [tempmax,newlist]=fun_findmax(oldEVlist(:,loop1:end));
%         oldEVlist=newlist;
%         result=[result,tempmax];
%     end
% 
%     function [max,newlist]=fun_findmax(list)
%         if size(list,2)==0
%             max=[];
%             newlist=[];
%         else
%             note=1;
%             max=list(:,1); 
%             for i=1:size(list,2)
%                 if list(2,i)>max(2,1)
%                     max=list(:,i);
%                     note=i;
%                 end
%             end
%             newlist=list;
%             newlist(:,[1 note])=newlist(:,[note 1]);
%         end     
%     end
% 
%     newEVlist=result;


result=[];
tempoldlist=oldEVlist;
if size(oldEVlist,2)==0
    result=[];
elseif size(oldEVlist,2)==1
    result=oldEVlist;
else
    num=size(oldEVlist,2);
    for ii=1:num-1
        [tempmax,note]=fun_findmax(tempoldlist);
        result=[result,tempmax];
        tempoldlist(:,note)=[];
    end
    result=[result,tempoldlist];
end
newEVlist=result;

    function [max,note]=fun_findmax(list)
        if size(list,2)==0
            max=[];
        else
            note=1;
            max=list(:,1); 
            for i=1:size(list,2)
                if list(2,i)>max(2,1)
                    max=list(:,i);
                    note=i;
                end
            end  
        end     
    end





end

