function[E] = z3_npl(D, N, Nu, DZ, lambda, draw)
    close all
    il = 5;

    ymin = 4.5;
    ymax = 34.5;
    a = 10;
    c = [8.5000   17.5000   24.5000   33.5000];
    hr0(1) = (c(1)+ymin)/2-1;
    hr0(il) = min((ymax+c(il-1))/2+1, ymax);
    if il > 2
        hr0(2:il-1) = (c(2:il-1)+c(1:il-2))./2;
    end

    % zmienne i macierze regulatora
    D = min(D, 400);
    DZ = min(DZ, 400);
    N = min(min(N,D),DZ);
    Nu = min(Nu, N);
    lambda = lambda*ones(1,il);

    M=zeros(N,Nu,il);
    Snr = zeros(D,il);
    for r = 1:il
        [s, z] = z1_step(hr0(r), false);
        Snr(:,r)=s(1:D);

        for i=1:N
           for j=1:Nu
              if (i>=j)
                 M(i,j,r)=s(i-j+1);
              end
           end
        end
    end
    deltaup=zeros(1,D-1)';
    deltazp=zeros(1,DZ-1)';
    w = zeros(1,il);
    duk = ones(1,Nu)';

    % dane
    C1 = 0.95;
    C2 = 0.95;
    a1 = 16;
    a2 = 16;
    tau = 50;
    U0 = 54;
    D0 = 10;
    Y0 = 16;
    V1 = C1*Y0^2;
    h1 = Y0;
    V2 = V1;
    start = 100;
    n = 12100;
    dY = [start 34; start+3000 4.5; start+6000 24; start+9000 14];
    dZ = [start+1500 15; start+4500 5; start+7500 12.5; start+10500 7.5];

    U = U0*ones(1,n);
    Dist = D0*ones(1,n);
    Y = Y0*ones(1,n);
    Yz = Y;
    for i = 1:length(dY)
        Yz(dY(i,1):n) = dY(i,2);
        Dist(dZ(i,1):n) = dZ(i,2);
    end
    e = zeros(1,n);

    hold on
    yzad = zeros(1,N)';
    yk = zeros(1,N)';
    A = [tril(ones(Nu));tril(ones(Nu))*-1];
    B = zeros(2*Nu,1);
    dk = 0;
    for k = start:n
        % symulacja
        V1 = V1 + U(k-1-tau) + Dist(k-1) - a1*h1^0.5;
        V2 = V2 + a1*h1^0.5 - a2*Y(k-1)^0.5;
        h1 = (V1/C1)^0.5;
        Y(k) = (V2/C2)^0.5;
        % uchyb
        e(k) = Yz(k) - Y(k);

        %uwzglêdnianie zak³ócenia
        for i = DZ:-1:2
           deltazp(i) = deltazp(i-1);
        end
        deltazp(1) = Dist(k) - Dist(k-1);

        % Prawo regulacji
        for i = 1:il
            if i == 1
                w(i) = 1-1/(1+exp(-a*(Y(k)-c(1))));
            elseif i == il
                w(i) = 1/(1+exp(-a*(Y(k)-c(il-1))));
            else
                w(i) = 1/(1+exp(-a*(Y(k)-c(i-1)))) - 1/(1+exp(-a*(Y(k)-c(i))));
            end
        end
        Mr=zeros(N,Nu);
        Sr = zeros(D,1);
        for i = 1:il
            Mr = Mr + w(i)*M(:,:,i)/sum(w);
            Sr = Sr + w(i)*Snr(:,i)/sum(w);
        end
        lambdar = w*lambda'/sum(w);
        y0 = zeros(N,1);
        dk = Y(k);
        if k<=D
            dk = dk-Sr(D)*U(1);
        else
            dk = dk-Sr(D)*U(k-D);
        end
        for i = 1:D-1
            if i<k-1
                dk=dk-Sr(i)*(U(k-i)-U(k-i-1));
            end
        end
        for i=1:N
            if k-1<D-i
                y0(i)=Sr(D)*U(1)+dk;
            else
                y0(i)=Sr(D)*U(k-D+i)+dk;
            end
            for j = i+1:D-1
                if j-i<k-1
                    y0(i)=y0(i)+Sr(j)*(U(k-j+i)-U(k-j+i-1));
                end
            end
        end
        B(1:Nu)=(90-U(k-1));
        B(Nu+1:end) = (U(k-1)-20);
        yzad(1:end)=Yz(k);
        duk = fmincon(@(duk)(yzad-y0-Mr*duk)'*(yzad-y0-Mr*duk)+lambdar*duk'*duk,duk,A,B);
        DELTAuk = duk(1);
        %duk = circshift(duk,-1);
        k

        for i = D-1:-1:2
          deltaup(i) = deltaup(i-1);
        end
        deltaup(1) = DELTAuk;
        U(k) = U(k-1)+deltaup(1);
        if draw
            plot(Yz, 'r')
            hold on
            plot(Y, 'b')
            plot(U,'g')
            drawnow;
        end
    end

    if draw
        figure
        subplot(3,1,1)
        plot(1:n, Yz, 'r')
        xlabel('t[s]')
        ylabel('h2[cm]')
        hold on
        plot(1:n, Y, 'b')
        subplot(3,1,2)
        plot(1:n, U)
        xlabel('t[s]')
        ylabel('F1in[cm3/s]')
        subplot(3,1,3)
        plot(1:n, Dist)
        xlabel('t[s]')
        ylabel('FD[cm3/s]')
    end
    E = sum(e.^2)/n;
end