function [U,uxx] = Equation(x,t)

U = ( 1 - 4*(1+4*1i*t)./( 1+4*x.^2+16*t.^2 ) ).*exp(2*1i*t);
uxx=1;
if t==0
    U = ( 1 - 4*(1+4*1i*t)./( 1+4*x.^2+16*t.^2 ) ).*exp(2*1i*t)+0.1*exp(-x.^2);
end
end