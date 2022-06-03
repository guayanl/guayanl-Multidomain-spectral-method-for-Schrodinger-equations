function [u,uxx] = Equation(x,t)
%% Values of function at the nodes x     u(x,0) 
%% Note: t is number ,not vector
% t = t*ones(size(x));
a = 2;
c = 15;
u = sqrt(a)*sech(sqrt(a)*(x-c*t)).*exp( 1i*( c/2*x+(a-c^2/4)*t ) );

uxx = -1/4*sqrt(a)*exp( 1/4*1i*(4*a*t+c*(2*x-c*t)) ).*sech(sqrt(a)*(x-c*t)).*( 4*a*(sech(sqrt(a)*(x-c*t))).^2  + (c+2i*sqrt(a)*tanh( sqrt(a)*(x-c*t) )).^2    )    ;

% u = (1./((1+4*1i*t).^(0.5))).*exp((-x.^2+8*1i*x-64*1i*t)./(1+4*1i*t));
% uxx = (1./((1+4*1i*t).^(2.5))).*exp((-x.^2+8*1i*x-64*1i*t)./(1+4*1i*t)).*(4*x.^2-32*1i*x-8*1i*t-66);


end