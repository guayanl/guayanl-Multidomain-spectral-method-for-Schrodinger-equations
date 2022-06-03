for index = [1,2]
    %% - - - - - Initialize Data - - - - -
    %% ги-infinity ,x_l) , (x_l, x_r) , (x_r,+infinity)     x_l<0 & x_r>0
    %%        N_1              N_2             N_3
    x_L = -5 ;x_R = 5 ;
    N1 = 200 ;N2 = 400 ;N3 = 200 ;
    t_start = 0 ;t_end = 1 ;Nt = 1000*index ;t = linspace(t_start,t_end,Nt)';h = (t_end-t_start)/(Nt-1);
    interval_number = 3;

    %% - - - - - Spatial Dependence - - - -   -1 <= L <= 1
    %% i U_t + A*U ( + V*U) = 0
    [A,U,X,Uxx] = Spatial_Discretization(x_L,x_R,N1,N2,N3,interval_number);

    %% - - - - - Boundary and matching conditions & Time integration- - - - -
    %% Fourth order implicit Runge-Kutta method
    Result = RK(A,U,t,h,N1,N2,N3,x_L,x_R);

    %% - - - - - - - - Calculation - - - - - - - -
    error = zeros(length(t),1);
    True_result = ones(length(U),length(t));
    
    for k = 1:length(t)
        [True_result(:,k),q] = Equation(X,t(k));
        True_result(1,k) = 0 ;True_result(end,k) = 0;
        temp = Result(:,k) - True_result(:,k);
        error(k) = norm(temp,inf)/norm(True_result(:,k),inf);
    end

    %% Plot
    semilogy(t,error);hold on;
end
legend('Nt = 1e3','Nt = 1e4')
xlabel('t');


