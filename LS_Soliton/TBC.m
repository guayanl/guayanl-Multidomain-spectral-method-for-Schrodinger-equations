 for index = [1,10]
    %% - - - - - - - - - - - - - - - - - - -Initialization- - - - - - - - - - - - - - - - - - -
    x_L = -5;x_R = 5 ;N1 = 0 ;N2 = 120 ;N3 = 0 ;
    t_start = 0 ;t_end = 0.5 ;Nt = 1000*index  ;t = linspace(t_start,t_end,Nt+1)';h = (t_end-t_start)/(Nt);
    interval_number = 1;
    [A,U,X,Uxx] = Spatial_Discretization(x_L,x_R,N1,N2,N3,interval_number);

    %% Crank-Nicolson method with Transparent Boundary Conditions
    % B*U(t_(n+1)) = C*U(t_n)
    B = eye(size(A))-1i*h/2*A ;B(1,:) = [1 zeros(1,N2)] ;B(end,:) = [zeros(1,N2) 1];
    C = eye(size(A))+1i*h/2*A ;
    [D,~] = Chebyshev_Differentiation_Matrix(N2);
    % F and Beta
    F = -exp(-1i*pi/4)*sqrt(2/h);
    beta = [1,-1];
    for k =1:Nt
        temp = beta(2*k-1)*(1-1/2/k);
        beta = [beta temp -temp]; 
    end
    
    BB = inv(B);
    Rrr = 2/(x_L-x_R)*( D(end,:)*BB(:,end) ) ;
    Rrl = 2/(x_L-x_R)*( D(end,:)*BB(:,1) ) ;
    Rlr = 2/(x_L-x_R)*( D(1,:)*BB(:,end) ) ;
    Rll = 2/(x_L-x_R)*( D(1,:)*BB(:,1) ) ;
    E = [Rlr,Rll+F*beta(1);Rrr-F*beta(1),Rrl] ;
    
    Result = 1000*ones(length(U),length(t));
    Result(:,1) = U;
    for k = 1:length(t)-1
        %% Construct the Dirichlet to Neumann map
        temp1 = Result(1,1:k) ;temp2 = Result(end,1:k);
        Vl = F*beta(2:1+k)*(fliplr(temp1)).';
        Vr = F*beta(2:1+k)*(fliplr(temp2)).';
        U_hat = C*Result(:,k);
        Ul = 2/(x_L-x_R)*D(1,:)*BB(:,2:end-1)*U_hat(2:end-1);
        Ur = 2/(x_L-x_R)*D(end,:)*BB(:,2:end-1)*U_hat(2:end-1);
        H = [-Ul-Vl;Vr-Ur];
        % g_0r g_0l
        g = E\H;
        U_hat(1) = g(2);
        U_hat(end) = g(1);
        %% Calculate the result
        x = B\U_hat;
        Result(:,k+1) = x;
    end



    %% - - - - - - - - Calculation - - - - - - - -
    error = zeros(1,length(t));
    True_result = ones(length(U),length(t));
    for k = 1:length(t)
        [True_result(:,k),useless] = Equation(X,t(k));
        temp = Result(:,k) - True_result(:,k);
        error(1,k) = norm(temp,2)/norm(True_result(:,k),2);
    end

    %% - - - - - Plot - - - - -   
    semilogy(t',error);hold on
 end
 legend('Nt=1e3','Nt=1e4','Nt=1e5');
 xlabel('t');

 