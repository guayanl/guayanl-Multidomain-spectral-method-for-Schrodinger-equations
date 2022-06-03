%% interpolation Node: cos(j*pi/N)   j=0,1,...,N
%% Chebyshev differentiation matrix
%% Cite  Trefethen L N. Spectral methods in MATLAB[M]. Society for industrial and applied mathematics, 2000.
function [D,x] = Chebyshev_Differentiation_Matrix(N)
if N == 0 ,D = 0 ;x = 1 ;return,end
x = cos((0:N)/N*pi)';
c = [2 ;ones(N-1,1) ;2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1)));
D = D-diag(sum(D'));
end
