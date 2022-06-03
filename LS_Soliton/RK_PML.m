function [Result] = RK_PML(A,U,X,t,h,N1,N2,N3,x_L,x_R,interval_number)
%% Fourth order implicit Runge-Kutta method
%% U(t_(n+1)) = U(t_n) + h*(0.5)*(K1+K2)
%% B*K1 = C*U(t_n)+B1*K2
%% B*K2 = C*U(t_n)+B2*K1
%% (K1+K2) = B^(-1)*C*U(t_n) + ()*
%% - - - - - Matching conditions and boundary - - - - -
delta = 0.5 ;
[D1,x1] = Chebyshev_Differentiation_Matrix(N1);
[D2,x2] = Chebyshev_Differentiation_Matrix(N2);
[D3,x3] = Chebyshev_Differentiation_Matrix(N3);

B = eye(size(A))-1i*h/4*A;
B(N1+1,:) = zeros(N1+N2+N3+3,1) ;B(N1+1,N1+1) = 1 ;B(N1+1,N1+2) = -1;
B(N1+2,1:N1+N2+2) = [2/(-delta)*D1(end,:) , -2/(x_L-x_R)*D2(1,:)];
B(N1+N2+2,:) = zeros(N1+N2+N3+3,1) ;B(N1+N2+2,N1+N2+2) = 1;B(N1+N2+2,N1+N2+3) = -1;
B(N1+N2+3,N1+2:end) = [2/(x_L-x_R)*D2(end,:) ,2/delta*D3(1,:)];

%% X_LP X_RP

B(1,:) = zeros(N1+N2+N3+3,1);
B(end,:) = zeros(N1+N2+N3+3,1);
B(1,1) = 1;
B(end,end) = 1;
A(1,:) = zeros(N1+N2+N3+3,1);
A(end,:) = zeros(N1+N2+N3+3,1);

%%
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
    %% x = BB*b;
    %% x = Matrix_Inversion(B,b,1e-10);
    CMY = [B,-B1;-B2,B];
    LGS = [b;b];
    K = CMY\LGS;
    x = Result(:,k) + h*0.5*(K(1:N1+N2+N3+3)+K(N1+N2+N3+4:end));
    Result(:,k+1) = x;
end




end