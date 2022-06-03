function [Result] = RK_PML(A,U,t,h,N1,N2,N3,x_L,x_R,delta)
%% Fourth order implicit Runge-Kutta method
%% U(t_(n+1)) = U(t_n) + h*(0.5)*(K1+K2)
%% B*K1 = C*U(t_n)+B1*K2
%% B*K2 = C*U(t_n)+B2*K1
%% (K1+K2) = B^(-1)*C*U(t_n) + ()*
%% - - - - - Matching conditions and boundary - - - - -
[D1,~] = Chebyshev_Differentiation_Matrix(N1);
[D2,~] = Chebyshev_Differentiation_Matrix(N2);
[D3,~] = Chebyshev_Differentiation_Matrix(N3);

B = eye(size(A))-1i*h/4*A;
B(N1+1,:) = zeros(N1+N2+N3+3,1) ;B(N1+1,N1+1) = 1 ;B(N1+1,N1+2) = -1;
B(N1+2,1:N1+N2+2) = [2/(-delta)*D1(end,:) , -2/(x_L-x_R)*D2(1,:)];
B(N1+N2+2,:) = zeros(N1+N2+N3+3,1) ;B(N1+N2+2,N1+N2+2) = 1;B(N1+N2+2,N1+N2+3) = -1;
B(N1+N2+3,N1+2:end) = [2/(x_L-x_R)*D2(end,:) ,2/delta*D3(1,:)];

B(1,:) = zeros(N1+N2+N3+3,1);
B(end,:) = zeros(N1+N2+N3+3,1);
B(1,1) = 1;
B(end,end) = 1;

%% Solve Problem by Simplified Netwon Method
Result = 1000*ones(length(U),length(t));
Result(:,1) = U;
for k = 1:(length(t)-1)
    b = Result(:,k);
    % K = K1+K2
    [K1,K2] = RK_Iteration(B,A,b,h,N1,N2,N3,1e-8);
    x = Result(:,k) + h*0.5*(K1+K2);
    Result(:,k+1) = x;
end


end