function [K1,K2] = RK_Iteration(B,C,B1,B2,b,precision)
%% B*K1 = b+B1*K2
%% B*K2 = b+B2*K1
K1 = zeros(size(b)); K1_temp = 100*ones(size(b));
K2 = zeros(size(b)); K2_temp = 100*ones(size(b));

while ( norm(K1_temp-K1,inf) > precision ) || (norm(K2_temp-K2,inf)>precision)     
    K1_temp = K1;
    K2_temp = K2;
    K1 = inv(B)*(b+B1*K2);
    K2 = inv(B)*(b+B2*K1);
end


end