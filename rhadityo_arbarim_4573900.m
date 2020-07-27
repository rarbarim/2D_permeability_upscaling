% upscaling permeability
% R B Arbarim : 4573900

clear all
format short

% Part 1 :
disp('Part 1 Number 1 :')
disp(' ')
n = ceil(str2double(input('enter number of blocks in vertical arrangement = ','s')));
m = ceil(str2double(input('enter number of blocks in horizontal arrangement = ','s')));

for i = 1 : n
    for j = 1 : m
       fprintf('block %s,%s \n',num2str(i),num2str(j))
       kg(i,j) = str2double(input('enter permeability in mD = ','s'));
    end
end

kgx = kg;
kgy = kg;
kg1x = zeros(n,1);
kg2x = zeros(1,n);
kg1y = zeros(1,n);
kg2y = zeros(n,1);
    
for i = 1 : n
    kg1x(i) = m*prod(kgx(i,1:m))/sum(kgx(i,:)); % harmonic average permeability in x direction
    kg2y(i) = sum(kgy(i,1:m))/(m); % aritmathic average permeability in y direction
end
    
for i = 1 : m
    kg2x(i) = sum(kgx(1:n,i))/(n); % aritmathic average permeability in x direction
    kg1y(i) = n*prod(kgy(1:n,i))/sum(kgy(:,i)); % harmonic average permeability in y direction
end
kg_lox = sum(kg1x)/(n); % lower bound permeability in x direction
kg_upx = n*prod(kg2x)/sum(kg2x); % upper bound permeability in x direction
kg_loy = sum(kg1y)/(m); % lower bound permeability in y direction
kg_upy = m*prod(kg2y)/sum(kg2y); % upper bound permeability in y direction
disp(' ')
fprintf('with configuration %u x %u blocks \n',n,m)
fprintf('permeability value in x-direction lies between %f and %f mD. \n',kg_lox,kg_upx)
fprintf('permeability value in y-direction lies between %f and %f mD. \n',kg_loy,kg_upy)

disp(' ')
disp('Part 1 Number 2 :')
disp(' ')

q = 2;
p = 2;
k = [1 100 ; 10 0.001 ];
kx = k;
ky = k;
k_4 = linspace(k(end,end),100,11);
k_1 = linspace(0,100,26);
k1x = zeros(q,1);
k2x = zeros(1,q);
k1y = zeros(1,q);
k2y = zeros(q,1);

for j = 1 : length(k_4)
    kx(end,end) = k_4(j); % sensitivity of light blue block in x-direction
    ky(end,end) = k_4(j); % sensitivity of light blue block in x-direction
    for i = 1 : q
        k1x(i) = p*prod(kx(i,1:p))/sum(kx(i,:)); % harmonic average permeability in x direction
        k2y(i) = sum(ky(i,1:p))/(p); % aritmathic average permeability in y direction
    end
    
    for i = 1 : p
        k2x(i) = sum(kx(1:q,i))/(q); % aritmathic average permeability in x direction
        k1y(i) = q*prod(ky(1:q,i))/sum(ky(:,i)); % harmonic average permeability in y direction
    end
    k_lox(j) = sum(k1x)/(q); % lower bound permeability in x direction
    k_upx(j) = q*prod(k2x)/sum(k2x); % upper bound permeability in x direction
    k_loy(j) = sum(k1y)/(p); % lower bound permeability in y direction
    k_upy(j) = p*prod(k2y)/sum(k2y); % upper bound permeability in y direction
end

fprintf('permeability in x-direction lies between %f and %f mD \n',k_lox(1),k_upx(1))
fprintf('permeability in y-direction lies between %f and %f mD \n',k_loy(1),k_upy(1))

disp(' ')
disp('Part 1 Number 4 is shown in figure(1)')
figure(1)
subplot (2,1,1)
hold on
plot(k_4,k_upx,'r')
plot(k_4,k_lox,'b')
ylabel('k-eff in x-direction (mD)')
xlabel('k light-blue block (mD)')
legend('upper bound permeability','lower bound permeability')
hold off

subplot (2,1,2)
hold on
plot(k_4,k_upy,'r')
plot(k_4,k_loy,'b')
ylabel('k-eff in y-direction (mD)')
xlabel('k light-blue block(mD)')
legend('upper bound permeability','lower bound permeability')
hold off

% Part 2 Number 1 :
disp(' ')
disp('Part 2 Number 1 : ')
k1_lox = block1(k);
k1_upx = block1(k);
k1_loy = block1(k);
k1_upy = block1(k);

k2_lox = block2(k);
k2_upx = block2(k);
k2_loy = block2(k);
k2_upy = block2(k);

k3_lox = block3(k);
k3_upx = block3(k);
k3_loy = block3(k);
k3_upy = block3(k);

figure(2)
subplot (2,2,1)
plot(k_1,k1_upx,'r')
ylabel('k-eff in x-direction (mD)')
xlabel('k sensitivity(mD)')
legend('upper bound permeability in x-direction')

subplot (2,2,2)
plot(k_1,k2_upx,'r')
ylabel('k-eff in x-direction (mD)')
xlabel('k sensitivity(mD)')
legend('upper bound permeability in x-direction')

subplot (2,2,3)
plot(k_1,k3_upx,'r')
ylabel('k-eff in x-direction (mD)')
xlabel('k sensitivity(mD)')
legend('upper bound permeability in x-direction')

subplot (2,2,4)
plot(k_4,k_upx,'r')
ylabel('k-eff in x-direction (mD)')
xlabel('k sensitivity(mD)')
legend('upper bound permeability in x-direction')

disp('block 1,1 is chosen as fractured block due to its high value of effective permeability in x-direction')
disp('as shown in figure(2)')
disp(' ')
k_efx = 200; % k-effective of fractured block in x-direction
k_efy = ky; % k-effective of fractured block in y-direction
b = 10; % fract aperture in micron
a = b*10^-6; % fract aperture in meter
a2 = 10^-6.*linspace(b,20,21);

k1xf = zeros(q,1);
k2xf = zeros(1,q);
k1yf = zeros(1,q);
k2yf = zeros(q,1);

kxf = kx;
kyf = ky;
kxf(1,1) = k_efx; % new value of kx in first block

for j = 1 : length(a2);
kf(j) = a2(j)^2/(12*10^-12)*1000; 
nf(j) = ceil((k_efx-kx(1,1))*1/(a*(kf(j)-kx(1,1))));

for i = 1 : q
     k1xf(i) = p*prod(kxf(i,1:p))/sum(kxf(i,:)); % harmonic average permeability in x direction
     k2yf(i) = sum(kyf(i,1:p))/(p); % aritmathic average permeability in y direction
end
 
for i = 1 : p
    k2xf(i) = sum(kxf(1:q,i))/(q); % aritmathic average permeability in x direction
    k1yf(i) = q*prod(kyf(1:q,i))/sum(kyf(:,i)); % harmonic average permeability in y direction
end
end

k_loxf = sum(k1xf)/(q); % lower bound permeability in x direction
k_upxf = q*prod(k2xf)/sum(k2xf); % upper bound permeability in x direction
k_loyf = sum(k1yf)/(p); % lower bound permeability in y direction
k_upyf = p*prod(k2yf)/sum(k2yf); % upper bound permeability in y direction

figure(3)
plot(10^6.*a2,nf)
xlabel('aperture (micron)');
ylabel('number of fractures');

disp('Part 2 Number 2 :')
fprintf('if we have effective permeability along fractured block 1,1 = %1.2f mD. \n',k_efx)
fprintf('we need %g number of %g micron width uniform fracture. \n',nf(1),b)
disp(' ')
disp('Part 2 Number 3 :')
fprintf('effective permeability in x-direction after fracturing block 1,1 lies between %1.2f mD to %1.2f mD. \n',k_loxf,k_upxf)
fprintf('effective permeability in y-direction after fracturing block 1,1 lies between %1.2f mD to %1.2f mD. \n',k_loyf,k_upyf)


