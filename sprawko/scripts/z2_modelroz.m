function [a, c, hr0] = z2_modelroz(il, draw)
    n = 400;
    C1 = 0.95;
    C2 = 0.95;
    a1 = 16;
    a2 = 16;
    tau = 50;
    start = 100;

    ymax = 34.5;
    ymin = 4.5;
    dy = (ymax-ymin)/il;
    a = 3;
    c = ymin+dy:dy:ymax-dy;
    hr0 = ones(1,il);
    hr0(1) = (c(1)+ymin)/2-1;
    hr0(il) = min((ymax+c(il-1))/2+1, ymax);
    if il > 2
        hr0(2:il-1) = (c(2:il-1)+c(1:il-2))./2;
    end
    Vr0 = C1*hr0.^2;
    Fr0 = a1*hr0.^0.5-10;
    if draw
        figure
        hold on
        for r = 1:il
            if r == 1
                plot(ymin:0.1:ymax,1-1./(1+exp(-a*([ymin:0.1:ymax]-c(1)))))
            elseif r == il
                plot(ymin:0.1:ymax,1./(1+exp(-a*([ymin:0.1:ymax]-c(il-1)))))
            else
                plot(ymin:0.1:ymax,1./(1+exp(-a*([ymin:0.1:ymax]-c(r-1))))-1./(1+exp(-a*([ymin:0.1:ymax]-c(r)))))
            end
        end
        plot(hr0, ones(1,il), 'ko')
        xlabel('h2')
        ylabel('przynale¿noœæ')
        title('Funkcje przynale¿noœci regulatorów lokalnych w regulacji rozmytej wzglêdem wartoœci wyjœcia')
    end
    F10 = 54;
    FD0 = 10;
    h10 = 16;
    h20 = 16;
    V10 = C1*h10*h10;
    V20 = C2*h20*h20;

    if draw
        figure
        title('Przebiegi wyjœcia dla skoku wartoœci sterowania w chwili 60s.')
        xlabel('czas[t]')
        ylabel('h2[cm]')
        hold on
    end
    for i = 84:-10:24
        F1in = F10 * ones(1,n);
        F1in(start:n) = i;
        F1 = F10 * ones(1,n);
        FD = FD0 * ones(1,n);
        h1 = h10 * ones(2+il+1,n);
        h2 = h20 * ones(2+il+1,n);
        V1 = C1*h10*h10 * ones(2+il+1,n);
        V2 = C2*h20*h20 * ones(2+il+1,n);
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
                V1(2+r,t) = V1(2+il+1,t-1) + (F1(t-1) - Fr0(r)) + (FD(t-1) - FD0) - a1/2*hr0(r)^-0.5 * (h1(2+il+1,t-1) - hr0(r));
                V2(2+r,t) = V2(2+il+1,t-1) + a1/2*hr0(r)^-0.5 * (h1(2+il+1,t-1) - hr0(r)) - a1/2*hr0(r)^-0.5 * (h2(2+il+1,t-1) - hr0(r));
                h1(2+r,t) = hr0(r) + 1/2*(C1*Vr0(r))^-0.5 * (V1(2+r,t) - Vr0(r));
                h2(2+r,t) = hr0(r) + 1/2*(C2*Vr0(r))^-0.5 * (V2(2+r,t) - Vr0(r));
                if r == 1
                    w(r) = 1-1/(1+exp(-a*(h2(2+il+1,t-1)-c(1))));
                elseif r == il
                    w(r) = 1/(1+exp(-a*(h2(2+il+1,t-1)-c(il-1))));
                else
                    w(r) = 1/(1+exp(-a*(h2(2+il+1,t-1)-c(r-1)))) - 1/(1+exp(-a*(h2(2+il+1,t-1)-c(r))));
                end
            end
            h2(2+il+1,t) = w*h2(3:2+il, t)/sum(w);
            h1(2+il+1,t) = w*h1(3:2+il, t)/sum(w);
            V1(2+il+1,t) = w*V1(3:2+il, t)/sum(w);
            V2(2+il+1,t) = w*V2(3:2+il, t)/sum(w);
        end
        if draw
            plot(start:n,h2(1,start:end),'b')
            plot(start:n,h2(2,start:end),'m')
            plot(start:n,h2(2+il+1,start:end),'g')
        end
    end
end

