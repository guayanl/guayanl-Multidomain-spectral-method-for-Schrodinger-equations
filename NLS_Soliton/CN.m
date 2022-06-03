function [Result] = CN(A,U,t,h,N1,N2,N3,x_L,x_R)
%% Crank-Nicolson   B * U(t_n+1) = C * U(t_n) = b£»
%% Boundary Conditions
[D1,~] = Chebyshev_Differentiation_Matrix(N1);
[D2,~] = Chebyshev_Differentiation_Matrix(N2);
[D3,~] = Chebyshev_Differentiation_Matrix(N3);
B = eye(size(A))-1i*h/2*A;
C = eye(size(A))+1i*h/2*A;


B(N1+1,:) = zeros(N1+N2+N3+3,1) ;B(N1+1,N1+1) = 1 ;B(N1+1,N1+2) = -1;
B(N1+2,1:N1+N2+2) = [2/x_L*D1(end,:) , -2/(x_L-x_R)*D2(1,:)];
B(N1+N2+2,:) = zeros(N1+N2+N3+3,1) ;B(N1+N2+2,N1+N2+2) = 1;B(N1+N2+2,N1+N2+3) = -1;
B(N1+N2+3,N1+2:end) = [2/(x_L-x_R)*D2(end,:) ,2/x_R*D3(1,:)];

C(N1+1,:) = zeros(N1+N2+N3+3,1);
C(N1+2,:) = zeros(N1+N2+N3+3,1);
C(N1+N2+2,:) = zeros(N1+N2+N3+3,1);
C(N1+N2+3,:) = zeros(N1+N2+N3+3,1);


%% Solve Problem by inverting the matrix
Result = 1000*ones(length(U),length(t));
Result(:,1) = U;
for k = 1:(length(t)-1)
    b = Result(:,k);
    x = CN_Iteration(B,C,b,h,N1,N2,N3,1e-8);
    Result(:,k+1) = x;
end


end