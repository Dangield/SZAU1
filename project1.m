clear variables
close all

n = 300;
C1 = 0.95;
C2 = 0.95;
a1 = 16;
a2 = 16;
tau = 50;

F10 = 54;
FD0 = 10;
h10 = 16;
h20 = 16;
V10 = C1*h10*h10;
V20 = C2*h20*h20;

figure
hold on
for i = 24:10:84
    F1in = i * ones(1,n);
    F1 = F10 * ones(1,n);
    FD = FD0 * ones(1,n);
    h1 = h10 * ones(2,n);
    h2 = h20 * ones(2,n);
    V1 = C1*h10*h10 * ones(2,n);
    V2 = C2*h20*h20 * ones(2,n);
    for t = tau+1 : n
        F1(t) = F1in(t-tau);
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
        V2(1,t) = V2(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
        h1(1,t) = (V1(1,t)/C1)^0.5;
        h2(1,t) = (V2(1,t)/C2)^0.5;
        
        V1(2,t) = V1(2,t-1) + (F1(t-1) - F10) + (FD(t-1) - FD0) - a1/2*h10^-0.5 * (h1(2,t-1) - h10);
        V2(2,t) = V2(2,t-1) + a1/2*h10^-0.5 * (h1(2,t-1) - h10) - a1/2*h20^-0.5 * (h2(2,t-1) - h20);
        h1(2,t) = h10 + 1/2*(C1*V10)^-0.5 * (V1(2,t) - V10);
        h2(2,t) = h20 + 1/2*(C2*V20)^-0.5 * (V2(2,t) - V20);
    end
    plot(h2(1,:), 'b')
    plot(h2(2,:), 'm')
end

figure
hold on
for i = 5:2.5:15
    F1in = 54 * ones(1,n);
    F1 = F10 * ones(1,n);
    FD = FD0 * ones(1,n);
    FD(1,50:n) = i;
    h1 = h10 * ones(2,n);
    h2 = h20 * ones(2,n);
    V1 = C1*h10*h10 * ones(2,n);
    V2 = C2*h20*h20 * ones(2,n);
    for t = tau+1 : n
        F1(t) = F1in(t-tau);
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - a1*h1(1,t-1)^0.5;
        V2(1,t) = V1(1,t-1) + a1*h1(1,t-1)^0.5 - a2*h2(1,t-1)^0.5;
        h1(1,t) = (V1(1,t)/C1)^0.5;
        h2(1,t) = (V2(1,t)/C2)^0.5;
        
        V1(2,t) = V1(2,t-1) + (F1(t-1) - F10) + (FD(t-1) - FD0) - a1/2*h10^-0.5 * (h1(2,t-1) - h10);
        V2(2,t) = V2(2,t-1) + a1/2*h10^-0.5 * (h1(2,t-1) - h10) - a1/2*h20^-0.5 * (h2(2,t-1) - h20);
        h1(2,t) = h10 + 1/2*(C1*V10)^-0.5 * (V1(2,t) - V10);
        h2(2,t) = h20 + 1/2*(C2*V20)^-0.5 * (V2(2,t) - V20);
    end
    plot(h2(1,:), 'b')
    plot(h2(2,:), 'm')
end