function [K1,K2] = RK_Iteration(B,A,b,h,N1,N2,N3,precision)
%% b = U(t_n)
%% B = E-ih*a11*L
%% A = L
K1 = zeros(size(b)); K1_temp = 100*ones(size(b));
K2 = zeros(size(b)); K2_temp = 100*ones(size(b));
a12 = 1/4 - sqrt(3)/6; a21 = 1/4 + sqrt(3)/6;
a11 = 1/4; a22 = 1/4;

while norm( [K1;K2]-[K1_temp;K2_temp] , inf )>precision
    K1_temp = K1; K2_temp = K2;
    V1 = 2*(( abs( b+h*a11*K1+h*a12*K2 ) ).^2); V2 = 2*( ( abs(b+h*a21*K1+h*a22*K2) ).^2);
    temp1 = 1i*A*b + 1i*h*a12*A*K2_temp + 1i*diag(V1)*(b + h*a11*K1_temp + h*a12*K2_temp);
    temp2 = 1i*A*b + 1i*h*a21*A*K1_temp + 1i*diag(V2)*(b + h*a21*K1_temp + h*a22*K2_temp);
    temp1(N1+1,:) = 0;temp1(N1+2,:) = 0;
    temp1(N1+N2+2,:) = 0;temp1(N1+N2+3,:) = 0;
    temp2(N1+1,:) = 0;temp2(N1+2,:) = 0;
    temp2(N1+N2+2,:) = 0;temp2(N1+N2+3,:) = 0;
    K1 = inv(B)*temp1;
    K2 = inv(B)*temp2;
end

end
