for index = [1]
    %% - - - - - Initialize Data - - - - -
    %% £¨-infinity ,x_l) , (x_l, x_r) , (x_r,+infinity)     x_l<0 & x_r>0
    %%        N_1              N_2             N_3
    x_L = -25 ;x_R = 25 ;
    N1 = 20 ;N2 = 700 ;N3 = 500 ;
    t_start = 0 ;t_end = 2 ;Nt = 1000*index ;t = linspace(t_start,t_end,Nt+1)';h = (t_end-t_start)/(Nt);
    interval_number = 3;

    %% - - - - - Spatial Dependence - - - -   -1 <= L <= 1
    %% i U_t + A*U ( + V*U) = 0
    [A,U,X,Uxx] = Spatial_Discretization(x_L,x_R,N1,N2,N3,interval_number);
    U(1) = 0;U(end)=0;
    Uxx(1) = 0;Uxx(end)=0;

    %% - - - - - Boundary and matching conditions & Time integration- - - - -
    %% Crank-Nicolson   B * U(t_n+1) = C * U(t_n) = b£»
    Result = CN(A,U,t,h,N1,N2,N3,x_L,x_R);
    %% Fourth order implicit Runge-Kutta method
    %% Result = RK(A,U,X,t,h,N1,N2,N3,x_L,x_R,interval_number);

    %% - - - - - - - - Calculation - - - - - - - -
    error = zeros(length(t),1);
    True_result = ones(length(U),length(t));
    
    for k = 1:length(t)
        [True_result(:,k),q] = Equation(X,t(k));
        True_result(1,k) = 0 ;True_result(end,k) = 0;
        temp = Result(:,k) - True_result(:,k);
        error(k) = norm(temp,2)/norm(True_result(:,k),2);
    end

    %% Plot
    semilogy(t,error(:,1));hold on;
end
legend('Nt = 1e3','Nt = 1e4')
xlabel('t');


