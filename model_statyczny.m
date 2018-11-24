clear variables
close all;
a1 = 16;
il = 5;
if il == 5
%     F1_0 = 29;%54
%     F1_1 = 43;%90
%     F1_2 = 54;%10
%     F1_3 = 67;%73
%     F1_4 = 79;%21
    F1 = [29 43 54 67 79];
    a = [0 7 12 17 23];
    b = [0 11 16 21 27];
    c = [7 12 17 23 Inf];
    d = [11 16 21 27 Inf];
end
if il == 4
%     F1_0 = 30;%54
%     F1_1 = 45;%90
%     F1_2 = 60;%10
%     F1_3 = 75;%73
%     F1_4 = 0;%21
    F1 = [32 45 61 77];
    a = [0 6 12.5 20.5];
    b = [0 12 18.5 26.5];
    c = [6 12.5 20.5 Inf];
    d = [12 18.5 26.5 Inf];
end
if il == 3
%     F1_0 = 35;%54
%     F1_1 = 54;%90
%     F1_2 = 75;%10
%     F1_3 = 0;%73
%     F1_4 = 0;%21
    F1 = [35 54 75];
    a = [0 9.5 18];
    b = [0 15.5 23];
    c = [9.5 18 Inf];
    d = [15.5 23 Inf];
end
if il == 2
%     F1_0 = 40;%54
%     F1_1 = 70;%90
%     F1_2 = 0;%10
%     F1_3 = 0;%73
%     F1_4 = 0;%21
    F1 = [40 70];
    a = [0 10.5];
    b = [0 22.5];
    c = [10.5 Inf];
    d = [22.5 Inf];
end
F_D = 10;
h2 = zeros (100,1);
% h2_l0 = h2;
% h2_l1 = h2;
% h2_l2 = h2;
% h2_l3 = h2;
% h2_l4 = h2;
h2_1 = zeros(il, 100);
h2_0 = zeros(100);
for i = 24:84
    h2(i)= ((i+F_D)/a1)^2;
end
for i = 1:il-1
    a(i+1) = h2(round((F1(i+1) + F1(i))/2)) - 3;
    b(i+1) = h2(round((F1(i+1) + F1(i))/2)) + 3;
    c(i) = h2(round((F1(i+1) + F1(i))/2)) - 3;
    d(i) = h2(round((F1(i+1) + F1(i))/2)) + 3;
end
%     h2_l0(i)= ((F1_0+F_D)/a1)^2+2*((F1_0+F_D)/a1)*(i-F1_0)/a1;
%     h2_l1(i)= ((F1_1+F_D)/a1)^2+2*((F1_1+F_D)/a1)*(i-F1_1)/a1;
%     h2_l2(i)= ((F1_2+F_D)/a1)^2+2*((F1_2+F_D)/a1)*(i-F1_2)/a1;
%     h2_l3(i)= ((F1_3+F_D)/a1)^2+2*((F1_3+F_D)/a1)*(i-F1_3)/a1;
%     h2_l4(i)= ((F1_4+F_D)/a1)^2+2*((F1_4+F_D)/a1)*(i-F1_4)/a1;
m = zeros(1, il);
for i = 24:84
    for j = 1:il
        h2_1(j,i) = ((F1(j)+F_D)/a1)^2+2*((F1(j)+F_D)/a1)*(i-F1(j))/a1;
        if h2_1(j,i) <= a(j)
            m(j) = 0;
        elseif h2_1(j,i) > a(j) && h2_1(j,i) < b(j)
            m(j) = (h2_1(j,i)-a(j))/(b(j)-a(j));
        elseif h2_1(j,i) >= b(j) && h2_1(j,i) <= c(j)
            m(j) = 1;
        elseif h2_1(j,i) > c(j) && h2_1(j,i) < d(j)
            m(j) = (d(j)-h2_1(j,i))/(d(j)-c(j));
        else
            m(j) = 0;
        end
    end
    h2_0(i) = m*h2_1(:,i)/sum(m);
end
plot(24:84,h2(24:84),'b');
hold on;
for i = 1:il
    if i == 1
        plot(24:round((F1(1)+F1(2))/2), h2_1(i, 24:round((F1(1)+F1(2))/2)),'r')
    elseif i == il
        plot(round((F1(il-1)+F1(il))/2):84, h2_1(i, round((F1(il-1)+F1(il))/2):84),'r')
    else
        plot(round((F1(i-1)+F1(i))/2):round((F1(i)+F1(i+1))/2), h2_1(i, round((F1(i-1)+F1(i))/2):round((F1(i)+F1(i+1))/2)),'r')
    end
end
% if il == 5
%     plot(24:uint8((F1_0+F1_1)/2),h2_l0(24:uint8((F1_0+F1_1)/2)),'r');
%     plot(uint8((F1_0+F1_1)/2):uint8((F1_1+F1_2)/2),h2_l1(uint8((F1_0+F1_1)/2):uint8((F1_1+F1_2)/2)),'r');
%     plot(uint8((F1_1+F1_2)/2):uint8((F1_2+F1_3)/2),h2_l2(uint8((F1_1+F1_2)/2):uint8((F1_2+F1_3)/2)),'r');
%     plot(uint8((F1_2+F1_3)/2):uint8((F1_3+F1_4)/2),h2_l3(uint8((F1_2+F1_3)/2):uint8((F1_3+F1_4)/2)),'r');
%     plot(uint8((F1_3+F1_4)/2):84,h2_l4(uint8((F1_3+F1_4)/2):84),'r');
% end
% if il == 4
%     plot(24:uint8((F1_0+F1_1)/2),h2_l0(24:uint8((F1_0+F1_1)/2)),'r');
%     plot(uint8((F1_0+F1_1)/2):uint8((F1_1+F1_2)/2),h2_l1(uint8((F1_0+F1_1)/2):uint8((F1_1+F1_2)/2)),'r');
%     plot(uint8((F1_1+F1_2)/2):uint8((F1_2+F1_3)/2),h2_l2(uint8((F1_1+F1_2)/2):uint8((F1_2+F1_3)/2)),'r');
%     plot(uint8((F1_2+F1_3)/2):84,h2_l3(uint8((F1_2+F1_3)/2):84),'r');
%     
% end
% if il == 3
%     plot(24:uint8((F1_0+F1_1)/2),h2_l0(24:uint8((F1_0+F1_1)/2)),'r');
%     plot(uint8((F1_0+F1_1)/2):uint8((F1_1+F1_2)/2),h2_l1(uint8((F1_0+F1_1)/2):uint8((F1_1+F1_2)/2)),'r');
%     plot(uint8((F1_1+F1_2)/2):84,h2_l2(uint8((F1_1+F1_2)/2):84),'r');
% end
% if il == 2
%     plot(24:uint8((F1_0+F1_1)/2),h2_l0(24:uint8((F1_0+F1_1)/2)),'r');
%     plot(uint8((F1_0+F1_1)/2):84,h2_l1(uint8((F1_0+F1_1)/2):84),'r');
% end
plot(24:84, h2_0(24:84), 'g');