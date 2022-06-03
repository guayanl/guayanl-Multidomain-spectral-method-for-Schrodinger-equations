
a = -10 ;b = 10 ;
N1 = 50 ;N2 = 700 ;N3 = 50;t = 1;
% Interpolation Nodes
L1 = cos((0:N1)/(N1)*pi)' ;L2 = cos((0:N2)/(N2)*pi)' ;L3 = cos((0:N3)/(N3)*pi)';
X1 = 2*a./(1-L1) ;X2 = a*(1+L2)/2+b*(1-L2)/2 ;X3 = 2*b./(1+L3);
U1 = Equation(X1,t) ;U2 = Equation(X2,t) ;U3 = Equation(X3,t);

%% Calculate coefficient of chebyshev expansion of polynomials
coef1 = Coef_matrix(L1,N1+1,0)\U1;coef2 = Coef_matrix(L2,N2+1,0)\U2 ;coef3 = Coef_matrix(L3,N3+1,0)\U3;

semilogy(abs(coef1),'b.');hold on; semilogy(abs(coef2),'g*');hold on; semilogy(abs(coef3),'rd');
xlim([0,600]);ylim([1e-20,1]);
handle = legend('Domain $\uppercase\expandafter{\romannumeral1}$','Domain $\uppercase\expandafter{\romannumeral2}$'...
,'$Domain \uppercase\expandafter{\romannumeral3}$');
set(handle,'Interpreter','latex')
xlabel('N');ylabel('log_{10}^{|a_{n}|}');