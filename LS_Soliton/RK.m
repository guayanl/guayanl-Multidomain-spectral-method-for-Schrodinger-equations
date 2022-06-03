function [Result] = RK(A,U,t,h,N1,N2,N3,x_L,x_R)
%% Fourth order implicit Runge-Kutta method
%% U(t_(n+1)) = U(t_n) + h*(0.5)*(K1+K2)
%% B*K1 = C*U(t_n)+B1*K2
%% B*K2 = C*U(t_n)+B2*K1

%% Initialization
[D1,~] = Chebyshev_Differentiation_Matrix(N1);
[D2,~] = Chebyshev_Differentiation_Matrix(N2);
[D3,~] = Chebyshev_Differentiation_Matrix(N3);

%% - - - - - Matching conditions and boundary - - - - -
B = eye(size(A))-1i*h/4*A;
B(N1+1,:) = zeros(N1+N2+N3+3,1) ;B(N1+1,N1+1) = 1 ;B(N1+1,N1+2) = -1;
B(N1+2,1:N1+N2+2) = [2/x_L*D1(end,:) , -2/(x_L-x_R)*D2(1,:)];
B(N1+N2+2,:) = zeros(N1+N2+N3+3,1) ;B(N1+N2+2,N1+N2+2) = 1;B(N1+N2+2,N1+N2+3) = -1;
B(N1+N2+3,N1+2:end) = [2/(x_L-x_R)*D2(end,:) ,2/x_R*D3(1,:)];

A(N1+1,:) = zeros(N1+N2+N3+3,1);
A(N1+2,:) = zeros(N1+N2+N3+3,1);
A(N1+N2+2,:) = zeros(N1+N2+N3+3,1);
A(N1+N2+3,:) = zeros(N1+N2+N3+3,1);
C = 1i*A;
B1 = 1i*h*(1/4-sqrt(3)/6)*A;
B2 = 1i*h*(1/4+sqrt(3)/6)*A;

%% Solve Problem by Simplified Netwon Method
Result = 1000*ones(length(U),length(t));
Result(:,1) = U;

for k = 1:(length(t)-1)
    b = C*Result(:,k);
    [K1,K2] = RK_Iteration(B,C,B1,B2,b,1e-8);
    x = Result(:,k) + h*0.5*(K1+K2);
    Result(:,k+1) = x;
end

end