function [alpha,beta,P_bound] = func_linearized(Diesel_data,width_num,dT)
    %P_boundはP_min,i,k^DG    燃料cost=alpha*出力+beta
    %% 発電機のLFC調整力を除いた上下限の設定
    Gen_N = size(Diesel_data,2);
    Pmax0 = Diesel_data(1,:)';
    Pmin0 = Diesel_data(2,:)';
    A0 = Diesel_data(3,:)';
    %A0 = abs(A0);
    B0 = Diesel_data(4,:)';
    C0 = Diesel_data(5,:)';
    dPre0 = Pmax0.*0;
    width = ( ( Pmax0 - dPre0 ) - ( Pmin0 + dPre0 ) ) ./width_num;%各発電機の線形補間幅 widthはdeltaPはず
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
線形補間関数
（燃料コスト） = アルファかける出力たすべーた
%}