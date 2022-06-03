function [A,U,X,Uxx] = Spatial_Discretization(x_L,x_R,N1,N2,N3,interval_number)
t_start = 0;
%% i U_t + A*U ( + V*U) = 0
%% A！！Modified Chebyshev differentiation matrix        U！！Values of function u at x
%  X！！Discrete points on the whole real line                 Uxx！！Second derivative of function u at x

%% Interval:(-infinity,x_1)
[D_1 ,L_1] = Chebyshev_Differentiation_Matrix(N1);
X_1 = 2*x_L./(1-L_1);
% Note: If u->0 when x-> \infinity,then u(inf)=0 in numerical calculation
[u_1,Uxx_1] = Equation(X_1,t_start);
A_1 = ((1-L_1).^4/4/(x_L)^2).*(D_1*D_1)-(1-L_1).^3/2/(x_L)^2.*D_1;

%% Interval:(x_l,x_r)    
[D_2,L_2] = Chebyshev_Differentiation_Matrix(N2);
X_2 = x_L*(1+L_2)/2+x_R*(1-L_2)/2;
[u_2,Uxx_2] = Equation(X_2,t_start);
A_2 = 4/(x_L-x_R).^(2)*D_2*D_2;
%% Uxx_cal_2 = D_2*u_2*2/(x_L-x_R);    formula of first derivative

%% Interval:(x_r,+infinity)
[D_3,L_3] = Chebyshev_Differentiation_Matrix(N3);
X_3 = 2*x_R./(1+L_3);
[u_3,Uxx_3] = Equation(X_3,t_start);
A_3 = ((1+L_3).^4/4/x_R^2).*D_3*D_3 +(1+L_3).^3/2/x_R^2.*D_3;

%% Combination
if interval_number == 3
    A = blkdiag(A_1,A_2,A_3);
    U = [u_1;u_2;u_3];
    X = [X_1;X_2;X_3];
    Uxx = [Uxx_1;Uxx_2;Uxx_3];
end

if interval_number == 1
    A = A_2;
    U = u_2;
    X = X_2;
    Uxx = Uxx_2;
end

end