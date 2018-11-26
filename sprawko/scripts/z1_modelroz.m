% function [] = z1_modelroz(il)
clear variables
close all
    
    n = 400;
    C1 = 0.95;
    C2 = 0.95;
    a1 = 16;
    a2 = 16;
    tau = 50;
    start = 100;

    il = 5;
    ymax = 34.5;
    ymin = 4.5;
    dy = (ymax-ymin)/(2*il-1);
    a = ymin*ones(1,il);
    b = ymin*ones(1,il);
    c = ymax*ones(1,il);
    d = ymax*ones(1,il);
    for r = 2:il
        a(r) = ymin+dy*(2*r-3);
        b(r) = ymin+dy*(2*r-2);
        c(r-1) = ymin+dy*(2*r-3);
        d(r-1) = ymin+dy*(2*r-2);
    end
%     a = [ymin    5.2656   10.5977   16.6914   23.9102];
%     b = [ymin   11.2656   16.5977   22.6914   29.9102];
%     c = [5.2656   10.5977   16.6914   23.9102    ymax];
%     d = [11.2656   16.5977   22.6914   29.9102       ymax];
    figure
    hold on
    for r = 1:il
        if r == 1
            plot([ymin c(r) d(r)],[1 1 0],'b')
        elseif r == il
            plot([a(r) b(r) ymax],[0 1 1],'b')
        else
            plot([a(r) b(r) c(r) d(r)],[0 1 1 0],'b')
        end
    end
    hr0 = (b+c)./2;
    hr0 = (a+d)./2;
    Vr0 = C1*hr0.^2;
    Fr0 = a1*hr0.^0.5-10;
%     Fr0 = [29 43 54 67 79];
%     hr0 = ((Fr0+10)/a1).^2;
%     Vr0 = C1*hr0.^2;
    a(1) = -Inf;
    b(1) = a(1);
    c(il) = Inf;
    d(il) = c(il);

    F10 = 54;
    FD0 = 10;
    h10 = 16;
    h20 = 16;
    V10 = C1*h10*h10;
    V20 = C2*h20*h20;

    figure
    title('Przebiegi wyjœcia dla skoku wartoœci sterowania w chwili 60s.')
    xlabel('czas[t]')
    ylabel('h2[cm]')
    hold on
    for i = 84:-10:24
        F1in = F10 * ones(1,n);
        F1in(start:n) = i;
        F1 = F10 * ones(1,n);
        FD = FD0 * ones(1,n);
        h1 = h10 * ones(2+il+1,n);
        h2 = h20 * ones(2+il+1,n);
        V1 = C1*h10*h10 * ones(2+il+1,n);
        V2 = C2*h20*h20 * ones(2+il+1,n);
%         htemp = ones(1,il);
        w = ones(1,il);
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
            
            for r = 1:il
%                 htemp(r) = hr0(r) + 1/2*(C2*Vr0(r))^-0.5 * (V1(3,t-1) + (F1(t-1) - Fr0(r)) +...
%                     (FD(t-1) - FD0) - a1/2*hr0(r)^-0.5 * (h1(3,t-1) - hr0(r)) - Vr0(r));
                V1(2+r,t) = V1(2+il+1,t-1) + (F1(t-1) - Fr0(r)) + (FD(t-1) - FD0) - a1/2*hr0(r)^-0.5 * (h1(2+il+1,t-1) - hr0(r));
                V2(2+r,t) = V2(2+il+1,t-1) + a1/2*hr0(r)^-0.5 * (h1(2+il+1,t-1) - hr0(r)) - a1/2*hr0(r)^-0.5 * (h2(2+il+1,t-1) - hr0(r));
                h1(2+r,t) = hr0(r) + 1/2*(C1*Vr0(r))^-0.5 * (V1(2+r,t) - Vr0(r));
                h2(2+r,t) = hr0(r) + 1/2*(C2*Vr0(r))^-0.5 * (V2(2+r,t) - Vr0(r));
                if h2(2+il+1,t-1) <= a(r)
                    w(r) = 0;
                elseif h2(2+il+1,t-1) > a(r) && h2(2+il+1,t-1) < b(r)
                    w(r) = (h2(2+il+1,t-1)-a(r))/(b(r)-a(r));
                elseif h2(2+il+1,t-1) >= b(r) && h2(2+il+1,t-1) <= c(r)
                    w(r) = 1;
                elseif h2(2+il+1,t-1) > c(r) && h2(2+il+1,t-1) < d(r)
                    w(r) = (d(r)-h2(2+il+1,t))/(d(r)-c(r));
                else
                    w(r) = 0;
                end
            end
%             h2(2+il+1,t) = w*htemp'/sum(w);
%             V2(2+il+1,t) = C2*h2(2+il+1,t)^2;
%             V1(2+il+1,t) = V1(3,t-1) + F1(t-1)+ FD(t-1) - a1*h1(3,t-1)^0.5;
%             h1(3,t) = (V1(3,t)/C1)^0.5;
            h2(2+il+1,t) = w*h2(3:2+il, t)/sum(w);
            h1(2+il+1,t) = w*h1(3:2+il, t)/sum(w);
            V1(2+il+1,t) = w*V1(3:2+il, t)/sum(w);
            V2(2+il+1,t) = w*V2(3:2+il, t)/sum(w);
        end
        plot(start:n,h2(1,start:end),'b')
        plot(start:n,h2(2,start:end),'m')
        plot(start:n,h2(2+il+1,start:end),'g')
    end
% end

