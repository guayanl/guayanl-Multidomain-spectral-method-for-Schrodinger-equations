function [T] = Coef_matrix(L,N,order)
%% L = [L_0,L_1,...L_(M-1)]
%% T_k(x) = cos(k arccos(x))
%% Order means different derivatives of T_k(x),the representation of following matrix omits the order 
%  T_0(L_0)     T_1(L_0)  T_2(L_0)  ... T_(N-1)(L_0)
%  T_0(L_1)     T_1(L_1)  T_2(L_1)  ... T_(N-1)(L_1)
%  T_0(L_2)     T_1(L_2)  T_2(L_2)  ... T_(N-1)(L_2)
%  ...
%  T_0(L_(M-1)) T_1(L_(M-1))  T_2(L_(M-1))  ... T_(N-1)(L_(M-1))
M = length(L);
T = zeros(M,N);
for k =1:N
    T(:,k) = Derivative_of_Tn(k-1,L,order);
end
end