

for ccc = [1,10,100]
    x_L = -5 ;x_R = 5 ;delta = 0.5 ;R = exp(1i*pi/4) ;sigma0 = 50;
    N1 = 20; N2 = 120; N3 = 50;
    %% Construct the spatial matrix
    [D1 ,L1] = Chebyshev_Differentiation_Matrix(N1);
    [D2 ,L2] = Chebyshev_Differentiation_Matrix(N2);
    [D3 ,L3] = Chebyshev_Differentiation_Matrix(N3);

    X1 = (x_L-delta)*(1+L1)/2+x_L*(1-L1)/2;
    X2 = x_L*(1+L2)/2+x_R*(1-L2)/2;
    X3 = x_R*(1+L3)/2+(x_R+delta)*(1-L3)/2;

    A1 = 4/(delta).^(2)*1./(1+R*sigma0*(X1-x_L).^2).^2.*D1*D1-2/(-delta)*( R*2*sigma0*(X1-x_L)./( 1+R*sigma0*(X1-x_L).^2).^3 ).*D1;
    A2 = D2*D2*4/(x_L-x_R).^(2);
    A3 = 4/(delta).^(2)*1./(1+R*sigma0*(X3-x_R).^2).^2.*D3*D3-2/(-delta)*( R*2*sigma0*(X3-x_R)./( 1+R*sigma0*(X3-x_R).^2).^3 ).*D3;

    X1_hat = X1+R*sigma0*(X1-x_L).^3/3;
    X2_hat = X2;
    X3_hat = X3+R*sigma0*(X3-x_R).^3/3;
    [U1,q] = Equation(X1_hat,0);
    [U2,q] = Equation(X2_hat,0);
    [U3,q] = Equation(X3_hat,0);

    A = blkdiag(A1,A2,A3);
    U = [U1;U2;U3];
    X = [X1;X2;X3];
    %% Time scheme
    t_start = 0 ;t_end = 0.5 ;Nt = 1e3*ccc;t = linspace(t_start,t_end,Nt)';h = (t_end-t_start)/(Nt-1) ;interval_number = 3;
    Result = CN_PML(A,U,t,h,N1,N2,N3,x_L,x_R);


    %% Calculate the true result
    error = zeros(length(t),1);
    True_result = ones(length(X2),length(t));

    for k = 1:length(t)
        [True_result(:,k),q] = Equation(X2,t(k));
    end
    for k = 1:length(t)
        temp = Result(N1+2:N1+N2+2,k) - True_result(:,k);
        error(k) = norm(temp,2)/norm(True_result(:,k),2);
    end

%% Plot
    semilogy(t,error(:,1));hold on;
end
% legend('40','50','60')
legend('1e3','1e4','1e5')
xlabel('t');




