function [Y] = objLin(U, D)
    load('state.mat');
    load('params.mat');
    
    V1 = V1 + (U - F10) + (D - FD0) - a1/2*h10^-0.5 * (h1 - h10);
    V2 = V2 + a1/2*h10^-0.5 * (h1 - h10) - a1/2*h20^-0.5 * (h2 - h20);
    h1 = h10 + 1/2*(C1*V10)^-0.5 * (V1 - V10);
    h2 = h20 + 1/2*(C2*V20)^-0.5 * (V2 - V20);
    
    Y = h2;
    save('state.mat', 'V1', 'V2', 'h1', 'h2')
end