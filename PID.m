clear variables
close all

% zmienne i macierze regulatora
load('stepRes')
D=length(s);

deltaup=zeros(1,D-1);

% dane
n = 10000;
tau = 50;
U0 = 54;
D0 = 10;
Y0 = 16;
dY = 5;
start = 60;

U = U0*ones(1,n);
D = D0*ones(1,n);
Y = Y0*ones(1,n);
Yz = Y;
%Yz(start:n) = Y0+dY;
Yz(start:2000) = Y0+dY*2;
Yz(2000:5000) = Y0-dY;
Yz(5000:n) = Y0;
%Yz(start:n) = Y0+1;
e = zeros(1,n);
Tp = 1;
kp = 0.7;
Ti = 37;
Td = 5;
r1 = kp*(1+Tp/(2*Ti)+Td/Tp);
r2 = kp*(Tp/(2*Ti)-2*Td/Tp-1);
r3 = kp*Td/Tp;

resetObj();
for k = start:n
    % symulacja
    Y(k) = obj(U(k-1-tau), D(k-1));
    % uchyb
    e(k) = Yz(k) - Y(k);
    
    U(k)=U(k-1)+r1*e(k)+r2*e(k-1)+r3*e(k-2);
    
end

plot(Yz, 'r')
hold on;
plot(Y, 'b')
hold off;
figure
plot(U)