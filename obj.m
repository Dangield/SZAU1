function [Y] = obj(U, D)
    load('state.mat');
    load('params.mat');
    
    V1 = V1 + U + D - a1*h1^0.5;
    V2 = V2 + a1*h1^0.5 - a2*h2^0.5;
    h1 = (V1/C1)^0.5;
    h2 = (V2/C2)^0.5;
    
    Y = h2;
    save('state.mat', 'V1', 'V2', 'h1', 'h2')
end

