function [x] = CN_Iteration(B,C,b,h,N1,N2,N3,precision)

temp = 2*abs(b).^2 ;V = diag(temp);
% V(1,:) = 0;
% V(end,:) = 0;
V(N1+1,:) = zeros(N1+N2+N3+3,1);
V(N1+2,:) = zeros(N1+N2+N3+3,1);
V(N1+N2+2,:) = zeros(N1+N2+N3+3,1);
V(N1+N2+3,:) = zeros(N1+N2+N3+3,1);

x = b;
x_temp = 1e8*ones(size(b));
while norm(x-x_temp,inf)>precision
    x_temp = x;
    temp = 2*abs(x).^2 ;
    V_temp = diag(temp) ;
% It doesn't matter
%     V_temp(1,:) = 0;
%     V_temp(end,:) = 0;
    V_temp(N1+1,:) = zeros(1,N1+N2+N3+3);
    V_temp(N1+2,:) = zeros(1,N1+N2+N3+3);
    V_temp(N1+N2+2,:) = zeros(1,N1+N2+N3+3);
    V_temp(N1+N2+3,:) = zeros(1,N1+N2+N3+3);
    x = inv(B)*C*b+inv(B)*1i*h/2*( V_temp*x+V*b );
    
end

end