clear variables
close all
resetObj();
load('state.mat');
load('params.mat');
load('stepRes');
% C1 = 0.95;
% C2 = 0.95;
% a1 = 16;
% a2 = 16;
il = 5;

h2_stat = zeros (100,1);
for u = 24:84
    h2_stat(u)= ((u+FD0)/a1)^2;
end
if il == 5
    F1r = [29 43 54 67 79];
elseif il == 4
    F1r = [32 45 61 77];
elseif il == 3
    F1r = [35 54 75];
elseif il == 2
    F1r = [40 70];
end
a = zeros(1,il);
b = zeros(1,il);
c = Inf*ones(1,il);
d = Inf*ones(1,il);
for j = 1:il-1
    a(j+1) = h2_stat(round((F1r(j+1) + F1r(j))/2)) - 3;
    b(j+1) = h2_stat(round((F1r(j+1) + F1r(j))/2)) + 3;
    c(j) = h2_stat(round((F1r(j+1) + F1r(j))/2)) - 3;
    d(j) = h2_stat(round((F1r(j+1) + F1r(j))/2)) + 3;
end

% zmienne i macierze regulatora
D=length(S(1,:));
N=D;
N=220;
Nu=D;
Nu=70;
DZ = length(Z);
lambda = 5000:40000/(il-1):45000;
if il == 2
    lambda = [1000 7000];
elseif il == 3
    lambda = [10 1000 3000];
elseif il == 4
    lambda = [10 100 100 1000];
elseif il == 5
    lambda = [10 100 100 1000 100];
end
    
M=zeros(N,Nu,il);
MP=zeros(N,D-1,il);
MZP=zeros(N,DZ,il);
for r = 1:il
    s = S(F1r(r)-23,:);
    z = Z(F1r(r)-23,:);
    
   
    for i=1:N
       for j=1:Nu
          if (i>=j)
             M(i,j,r)=s(i-j+1);
          end
       end
    end

    
    for i=1:N
       for j=1:D-1
          if i+j<=D
             MP(i,j,r)=s(i+j)-s(j);
          else
             MP(i,j,r)=s(D)-s(j);
          end    
       end
    end

    
    for i=1:N
        MZP(i,1,r) = z(i);
       for j=2:DZ
          if i+j-1<=DZ
             MZP(i,j,r)=z(i+j-1)-z(j);
          else
             MZP(i,j,r)=z(DZ)-z(j);
          end      
       end
    end
end
deltaup=zeros(1,D-1)';
deltazp=zeros(1,DZ-1)';
w = zeros(1,il);
duk = ones(1,Nu)';

% dane
tau = 50;
U0 = 54;
D0 = 10;
Y0 = 16;
start = 100;
dY = [start 34; start+2000 4.5; start+4000 24; start+6000 14];
n = length(dY)*2000+100;
n = 800;

U = U0*ones(1,n);
Dist = D0*ones(1,n);
Y = Y0*ones(1,n);
Yz = Y;
for i = 1:length(dY)
    Yz(dY(i,1):n) = dY(i,2);
end
e = zeros(1,n);

hold on
yzad = zeros(1,N)';
yk = zeros(1,N)';
A = [tril(ones(Nu));tril(ones(Nu))*-1];
b = zeros(2*Nu,1);
for k = start:n
    % symulacja
    V1 = V1 + U(k-1-tau) + Dist(k-1) - a1*h1^0.5;
    V2 = V2 + a1*h1^0.5 - a2*h2^0.5;
    h1 = (V1/C1)^0.5;
    h2 = (V2/C2)^0.5;
    Y(k) = h2;
    % uchyb
    e(k) = Yz(k) - Y(k);
    
    %uwzglêdnianie zak³ócenia
    for i = DZ:-1:2
       deltazp(i) = deltazp(i-1);
    end
    deltazp(1) = Dist(k) - Dist(k-1);

    % Prawo regulacji
    for i = 1:il
        if Y(k) <= a(i)
            w(i) = 0;
        elseif Y(k) > a(i) && Y(k) < b(i)
            w(i) = (Y(k)-a(i))/(b(i)-a(i));
        elseif Y(k) >= b(i) && Y(k) <= c(i)
            w(i) = 1;
        elseif Y(k) > c(i) && Y(k) < d(i)
            w(i) = (d(i)-Y(k))/(d(i)-c(i));
        else
            w(i) = 0;
        end
    end
    Mr=zeros(N,Nu);
    MPr=zeros(N,D-1);
    MZPr=zeros(N,DZ);
    lambdar = 0;
    for i = 1:il
        Mr = Mr + w(i)*M(:,:,i)/sum(w);
        MPr = MPr + w(i)*MP(:,:,i)/sum(w);
        MZPr = MZPr + w(i)*MZP(:,:,i)/sum(w);
    end
    lambdar = w*lambda'/sum(w);
    b(1:Nu)=(84-U(k-1));
    b(Nu+1:end) = (U(k-1)-24);
    yzad(1:end)=Yz(k);
    yk(1:end)=Y(k);
    duk = fmincon(@(duk)(yzad-yk+MPr*deltaup+MZPr*deltazp-Mr*duk)'*(yzad-yk+MPr*deltaup+MZPr*deltazp-Mr*duk)+lambdar*duk'*duk,duk,A,b,[],[],ones(Nu,1)*-60,ones(Nu,1)*60);
    DELTAuk = duk(1);
    %duk = circshift(duk,-1);
    k

    for i = D-1:-1:2
      deltaup(i) = deltaup(i-1);
    end
    deltaup(1) = DELTAuk;
    U(k) = U(k-1)+deltaup(1);
end

plot(Yz, 'r')
hold on
plot(Y, 'b')
% figure
% plot(U)