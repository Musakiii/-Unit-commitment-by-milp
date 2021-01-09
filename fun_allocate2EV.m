function [up,down] = fun_allocate2EV(signal,EVlist)
%UNTITLED11 此处显示有关此函数的摘要
%
num=size(EVlist,2);
if signal>=0
    temp=fix(signal/4)*4;       
else
    temp=-1*fix(abs(signal)/4)*4;    
end
up=floor((temp+4*num)/8);
down=num-up; 

end

