%% Equation:  i u_t(x) + u_xx(x) = 0          u: R¡ÁR -> C
%%   Method:  Compactified exterior domains(CED) 
%%  Spatial:  Chebyshev collocation methods
%%     Time:  CN/RK Method
%%     Date: 2022.2.12 by Francis
%% -----------------------------------------------------------------------------------------------------------------

for index = [1,10]
    %% - - - - - Initialize Data - - - - -
    %% £¨-infinity ,x_l) , (x_l, x_r) , (x_r,+infinity)     x_l<0 & x_r>0
    %%        N_1              N_2             N_3
    x_L = -5 ;x_R = 5 ;
    N1 = 20 ;N2 = 120 ;N3 = 600 ;
    t_start = 0 ;t_end = 0.5 ;Nt = 1000*index ;t = linspace(t_start,t_end,Nt+1)';h = (t_end-t_start)/(Nt);
    interval_number = 3;

    %% - - - - - Spatial Dependence - - - -   -1 <= L <= 1
    %% i U_t + A*U ( + V*U) = 0
    [A,U,X,Uxx] = Spatial_Discretization(x_L,x_R,N1,N2,N3,interval_number);
    
    %% - - - - - Boundary and matching conditions & Time integration- - - - -
    %% Crank-Nicolson   B * U(t_n+1) = C * U(t_n) = b£»
     %Result = CN(A,U,t,h,N1,N2,N3,x_L,x_R);

    %% Fourth order implicit Runge-Kutta method
     Result = RK(A,U,t,h,N1,N2,N3,x_L,x_R);
     
    %% - - - - - - - - Calculation - - - - - - - -
    error = zeros(length(t),1);
    True_result = ones(length(U),length(t));
    % generate the true_result matrix
    for k = 1:length(t)
        [True_result(:,k),q] = Equation(X,t(k));
        True_result(1,k) = 0 ;True_result(end,k) = 0;
    end
    % calculate the error
    for k = 1:length(t)
        temp = Result(:,k) - True_result(:,k);
        error(k) = norm(temp,2)/norm(True_result(:,k),2);
    end
    %% Plot the image
    semilogy(t,error(:,1));hold on;
end
legend('Nt = 1e3','Nt = 1e4')
xlabel('t');

