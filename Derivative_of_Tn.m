function [derivative] = Derivative_of_Tn(n,L,order)
%% Calculate the different derivatives of Chebyshev polynomial
%% L = [1,....,-1]
%% Tn(L)   /   Tn'(L)   /   Tn''(L)
derivative = zeros(size(L));
if order == 0
    derivative = cos(n*acos(L));
elseif L(1) == 1 && L(end) == -1 && order == 2
    L = L(2:end-1);
    derivative(1) = (n^4-n^2)/3;
    derivative(end) = (-1)^n*(n^4-n^2)/3;
    derivative(2:end-1,1) = n*((n+1)*cos(n*acos(L))-sin((n+1)*acos(L))./sin(acos(L)))./(L.^2-1);
    
elseif L(1) == 1 && L(end) ~= -1 && order == 2
    L = L(2:end);
    derivative(1) = (n^4-n^2)/3;
    derivative(2:end,1) = n*((n+1)*cos(n*acos(L))-sin((n+1)*acos(L))./sin(acos(L)))./(L.^2-1);
    
    
elseif L(1) ~=1 && L(end) == -1 && order == 2
    L = L(1:end-1);
    derivative(end) = (-1)^n*(n^4-n^2)/3;
    derivative(1:end-1,1) = n*((n+1)*cos(n*acos(L))-sin((n+1)*acos(L))./sin(acos(L)))./(L.^2-1);
    
elseif L(1) ~=1 && L(end) ~= -1 && order == 2
    derivative = n*((n+1)*cos(n*acos(L))-sin((n+1)*acos(L))./sin(acos(L)))./(L.^2-1);
    
elseif L(1) == 1 && L(end) == -1 && order == 1
    L = L(2:end-1);
    derivative(1) = n^2;
    derivative(end) = n^2*(-1)^(n-1);
    derivative(2:end-1,1) = n*sin(n*acos(L))./sin(acos(L));
elseif L(1) == 1 && L(end) ~= -1 && order == 1
    L = L(2:end);
    derivative(1) = n^2;
    derivative(2:end,1) = n*sin(n*acos(L))./sin(acos(L));
elseif L(1) ~= 1 && L(end) == -1 && order == 1
    L = L(1:end-1);
    derivative(end) =  n^2*(-1)^(n-1);
    derivative(1:end-1,1) = n*sin(n*acos(L))./sin(acos(L));
    
else
    derivative = n*sin(n*acos(L))./sin(acos(L));
end

end
