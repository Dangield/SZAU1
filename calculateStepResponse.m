clear variables
close all

n = 400;
tau = 50;
U0 = 54;
D0 = 10;
Y0 = 16;
start = 100;

U = U0*ones(1,n);
D = D0*ones(1,n);
Y = Y0*ones(1,n);

step = 10;
U(1,start:n) = U0 + step;
resetObj();
for i = start+1:n
    Y(i) = objLin(U(i-1-tau), D(i-1));
end
s = (Y(start+1:start+225)-Y0)/step;
subplot(2,1,1)
plot(s)

U = U0*ones(1,n);
D = D0*ones(1,n);
Y = Y0*ones(1,n);

step = 1;
D(1,start:n) = D0 + step;
resetObj();
for i = start+1:n
    Y(i) = objLin(U(i-1-tau), D(i-1));
end
z = (Y(start+1:start+225)-Y0)/step;
subplot(2,1,2)
plot(z)

save('stepRes.mat', 's', 'z')