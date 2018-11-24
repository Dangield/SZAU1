clear variables
close all

n = 600;
tau = 50;
U0 = 54;
D0 = 10;
Y0 = 16;
start = 100;

% U = U0*ones(1,n);
% D = D0*ones(1,n);
% Y = Y0*ones(1,n);

% step = 10;
% U(1,start:n) = U0 + step;
% resetObj();
% for i = start+1:n
%     Y(i) = objLin(U(i-1-tau), D(i-1));
% end
% s = (Y(start+1:start+225)-Y0)/step;

% U = U0*ones(1,n);
% D = D0*ones(1,n);
% Y = Y0*ones(1,n);

step = 1;
% D(1,start:n) = D0 + step;
% resetObj();
% for i = start+1:n
%     Y(i) = objLin(U(i-1-tau), D(i-1));
% end
% z = (Y(start+1:start+225)-Y0)/step;

S = zeros(61, n-start);
Z = zeros(61, n-start);
load('params')
for U0 = 24:84
    h10 = ((U0 + D0)/a1)^2;
    V10 = C1*h10^2;
    h20 = h10;
    V20 = V10;
    U = U0*ones(1,n);
    D = D0*ones(1,n);
    U(1,start:n) = U0 + step;
    Y = h20*ones(1,n);
    V1 = V10;
    V2 = V20;
    h1 = h10;
    
    for i = start+1:n
        V1 = V1 + (U(i-1-tau) - U0) + (D(i-1) - D0) - a1/2*h10^-0.5 * (h1 - h10);
        V2 = V2 + a1/2*h10^-0.5 * (h1 - h10) - a1/2*h20^-0.5 * (Y(i-1) - h20);
        h1 = h10 + 1/2*(C1*V10)^-0.5 * (V1 - V10);
        Y(i) = h20 + 1/2*(C2*V20)^-0.5 * (V2 - V20);
    end
    S(U0-23,:) = (Y(start+1:n)-h20)/step;
end
for U0 = 24:84
    h10 = ((U0 + D0)/a1)^2;
    V10 = C1*h10^2;
    h20 = h10;
    V20 = V10;
    U = U0*ones(1,n);
    D = D0*ones(1,n);
    D(1,start:n) = D0 + step;
    Y = h20*ones(1,n);
    V1 = V10;
    V2 = V20;
    h1 = h10;
    
    for i = start+1:n
        V1 = V1 + (U(i-1-tau) - U0) + (D(i-1) - D0) - a1/2*h10^-0.5 * (h1 - h10);
        V2 = V2 + a1/2*h10^-0.5 * (h1 - h10) - a1/2*h20^-0.5 * (Y(i-1) - h20);
        h1 = h10 + 1/2*(C1*V10)^-0.5 * (V1 - V10);
        Y(i) = h20 + 1/2*(C2*V20)^-0.5 * (V2 - V20);
    end
    Z(U0-23,:) = (Y(start+1:n)-h20)/step;
end
s = S(31,1:225);
z = Z(31,1:225);
subplot(2,1,1)
plot(s)
subplot(2,1,2)
plot(z)
figure
hold on
plot(S(31,:),'b')
plot(S(1,:),'r')
plot(S(61,:),'g')
figure
hold on
plot(Z(31,:),'b')
plot(Z(1,:),'r')
plot(Z(61,:),'g')
save('stepRes.mat', 's', 'z', 'S', 'Z')