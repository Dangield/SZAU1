clear variables
close all

% zmienne i macierze regulatora
load('stepRes')
D=length(s);
N=D;
Nu=D;
DZ = length(z);
lambda = 1;

M=zeros(N,Nu);
for i=1:N
   for j=1:Nu
      if (i>=j)
         M(i,j)=s(i-j+1);
      end
   end
end

MP=zeros(N,D-1);
for i=1:N
   for j=1:D-1
      if i+j<=D
         MP(i,j)=s(i+j)-s(j);
      else
         MP(i,j)=s(D)-s(j);
      end    
   end
end

MZP=zeros(N,DZ);
for i=1:N
    MZP(i,1) = z(i);
   for j=2:DZ
      if i+j-1<=DZ
         MZP(i,j)=z(i+j-1)-z(j);
      else
         MZP(i,j)=z(DZ)-z(j);
      end      
   end
end

I=eye(Nu);
K=((M'*M+lambda*I)^-1)*M';
ku=K(1,:)*MP;
kz=K(1,:)*MZP;
ke=sum(K(1,:));
deltaup=zeros(1,D-1);
deltazp=zeros(1,DZ-1);

% dane
n = 400;
tau = 50;
U0 = 54;
D0 = 10;
Y0 = 16;
dY = 5;
start = 100;

U = U0*ones(1,n);
D = D0*ones(1,n);
Y = Y0*ones(1,n);
Yz = Y;
Yz(start:n) = Y0+dY;
e = zeros(1,n);

resetObj();
for k = start:n
    % symulacja
    Y(k) = obj(U(k-1-tau), D(k-1));
    % uchyb
    e(k) = Yz(k) - Y(k);
    
    %uwzglêdnianie zak³ócenia
    for i = DZ:-1:2
       deltazp(i) = deltazp(i-1);
    end
    deltazp(1) = D(k) - D(k-1);

    % Prawo regulacji
    deltauk = ke*e(k)-ku*deltaup'-kz*deltazp';

    for i = D-1:-1:2
      deltaup(i) = deltaup(i-1);
    end
    deltaup(1) = deltauk;
    U(k) = U(k-1)+deltaup(1);
end

plot(Yz, 'r')
hold on
plot(Y, 'b')