%% Initialization
a = -5 ;b = 5 ;N1 = 100 ;N2 = 120 ;N3 = 5000;t = 0.5;
% Interpolation Nodes
L1 = cos((0:N1)/(N1)*pi)' ;L2 = cos((0:N2)/(N2)*pi)' ;L3 = cos((0:N3)/(N3)*pi)';
X1 = 2*a./(1-L1) ;X2 = a*(1+L2)/2+b*(1-L2)/2 ;X3 = 2*b./(1+L3);
[U1,~] = Equation(X1,t) ;[U2,~] = Equation(X2,t) ;[U3,~] = Equation(X3,t);
% values of function f at infinity are zero
U1(1) = 0;U3(end)=0;

%% Calculate coefficient of chebyshev expansion of polynomials
coef1 = Coef_matrix(L1,N1+1,0)\U1;coef2 = Coef_matrix(L2,N2+1,0)\U2 ;coef3 = Coef_matrix(L3,N3+1,0)\U3;

%% Plot the image
semilogy(abs(coef1),'.');hold on; semilogy(abs(coef2),'g*');hold on; semilogy(abs(coef3),'r--');
xlim([0,600]);ylim([1e-20,1]);
handle = legend('Domain $\uppercase\expandafter{\romannumeral1}$','Domain $\uppercase\expandafter{\romannumeral2}$'...
,'$Domain \uppercase\expandafter{\romannumeral3}$');
set(handle,'Interpreter','latex')
ylabel([1e-16,0])