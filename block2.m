function [k2_lox,k2_upx,k2_loy,k2_upy] = block2(k)
    m = 2;
    n = 2;
    kx = k;
    ky = k;
    k_1 = linspace(0,100,26);
    k1x = zeros(n,1);
    k2x = zeros(1,n);
    k1y = zeros(1,n);
    k2y = zeros(n,1);
    
for j = 1 : length(k_1)
    kx(1,end) = k_1(j); % sensitivity of light blue block in x-direction
    ky(1,end) = k_1(j); % sensitivity of light blue block in x-direction
    for i = 1 : n
        k1x(i) = m*prod(kx(i,1:m))/sum(kx(i,:)); % harmonic average permeability in x direction
        k2y(i) = sum(ky(i,1:m))/(m); % aritmathic average permeability in y direction
    end
    
    for i = 1 : m
        k2x(i) = sum(kx(1:n,i))/(n); % aritmathic average permeability in x direction
        k1y(i) = n*prod(ky(1:n,i))/sum(ky(:,i)); % harmonic average permeability in y direction
    end
    k2_lox(j) = sum(k1x)/(n); % lower bound permeability in x direction
    k2_upx(j) = n*prod(k2x)/sum(k2x); % upper bound permeability in x direction
    k2_loy(j) = sum(k1y)/(m); % lower bound permeability in y direction
    k2_upy(j) = m*prod(k2y)/sum(k2y); % upper bound permeability in y direction
end
