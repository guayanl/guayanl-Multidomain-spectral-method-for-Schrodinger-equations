function [u,uxx] = Equation(x,t)
% x   Values of function at the nodes x     u(x,0) 
% t   time
% u   value of function at points x
% uxx second derivative of function u at x
t = t*ones(size(x));
u = (1./((1+4*1i*t).^(0.5))).*exp((-x.^2+8*1i*x-64*1i*t)./(1+4*1i*t));
uxx = (1./((1+4*1i*t).^(2.5))).*exp((-x.^2+8*1i*x-64*1i*t)./(1+4*1i*t)).*(4*x.^2-32*1i*x-8*1i*t-66);

end