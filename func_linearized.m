function [alpha,beta,P_bound] = func_linearized(Diesel_data,width_num,dT)
    %P_bound��P_min,i,k^DG    �R��cost=alpha*�o��+beta
    %% ���d�@��LFC�����͂��������㉺���̐ݒ�
    Gen_N = size(Diesel_data,2);
    Pmax0 = Diesel_data(1,:)';
    Pmin0 = Diesel_data(2,:)';
    A0 = Diesel_data(3,:)';
    %A0 = abs(A0);
    B0 = Diesel_data(4,:)';
    C0 = Diesel_data(5,:)';
    dPre0 = Pmax0.*0;
    width = ( ( Pmax0 - dPre0 ) - ( Pmin0 + dPre0 ) ) ./width_num;%�e���d�@�̐��`��ԕ� width��deltaP�͂�
    for loop1 = 1:Gen_N,
       for loop2 = 1:width_num,
           low_cost   = C0(loop1,1).*dT + B0(loop1,1).*dT*( Pmin0(loop1,1) + dPre0(loop1,1) + (loop2-1)*width(loop1,1) ) + A0(loop1,1).*dT*( Pmin0(loop1,1) + dPre0(loop1,1) + (loop2-1)*width(loop1,1) )^2;
           high_cost  = C0(loop1,1).*dT + B0(loop1,1).*dT*( Pmin0(loop1,1) + dPre0(loop1,1) + (loop2)*width(loop1,1) ) + A0(loop1,1).*dT*( Pmin0(loop1,1) + dPre0(loop1,1) + (loop2)*width(loop1,1) )^2;
           alpha(loop1,loop2) = (high_cost-low_cost)/width(loop1,1);
           beta(loop1,loop2) = low_cost - alpha(loop1,loop2)*( Pmin0(loop1,1) + dPre0(loop1,1) + (loop2-1)*width(loop1,1) );
           P_bound(loop1,loop2) =  Pmin0(loop1,1) + dPre0(loop1,1) + (loop2-1)*width(loop1,1);
       end
       P_bound(loop1,width_num+1) = Pmin0(loop1,1) + dPre0(loop1,1) + (width_num)*width(loop1,1);
    end
end

%{
���`��Ԋ֐�
�i�R���R�X�g�j = �A���t�@������o�͂����ׁ[��
%}