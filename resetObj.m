function [] = resetObj()
    C1 = 0.95;
    C2 = 0.95;
    a1 = 16;
    a2 = 16;
    
    h1 = 16;
    h2 = 16;
    V1 = C1*h1*h1;
    V2 = C2*h2*h2;
    
    h10 = h1;
    h20 = h2;
    V10 = V1;
    V20 = V2;
    F10 = 54;
    FD0 = 10;
    
    save('params.mat', 'C1', 'C2', 'a1', 'a2', 'h10', 'h20', 'V10', 'V20', 'F10', 'FD0')
    save('state.mat', 'h1', 'h2', 'V1', 'V2')
end

