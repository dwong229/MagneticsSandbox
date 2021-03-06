%% Potential Function Plot

close all
clear all

R = 0.01; % radius of coil, 1cm [m]
L = .1; % separation between 2 coils [m];
u = 0.05; %Am^2
u0 = 4*pi*10^-7; % Tm/A

%% Given x1 and x2 solve for I1 and I2
inputx1 = 0.03;
inputx2 = 0.05;

[I1,I2] = solveforcurrent(inputx1,inputx2);

fprintf('For a coil separation of L = %2.2f m with coils of radius R = %2.2f \n',L,R)
fprintf('To hold 2 magnets at x = %3.2f and %3.2f m\n',inputx1,inputx2)
fprintf('Current needed: I1 = %7.4f A, I2 = %7.4f A \n',I1,I2)

%% Compute Force
idx = 1;
x1vec = 0.01:0.001:inputx2-0.01;
x2 = inputx2;
Ivec = [I1;I2;1];
xlength = length(x1vec);
F = zeros(2,xlength);
for idx = 1:xlength
    x1 = x1vec(idx);
    C(1) = -3*u*x1/(R^2+x1^2)^(5/2)*u0/2*R^2;
    C(2) = -3*(x1-L)*u/(R^2 + (x1-L)^2)^(5/2)*u0/2*R^2;
    C(3) = u0*u*u/(4*pi*(x1-x2)^2);
    C(4) = -3*u*(x1-L)/(R^2+(x1-L)^2)^(5/2)*u0/2*R^2;
    C(5) = -3*(x2-L)*u/(R^2 + (x2-L)^2)^(5/2)*u0/2*R^2;
    C(6) = -u0*u*u/(4*pi*(x1-x2)^2);
    A = [C(1) C(2) C(3);C(4) C(5) C(6)];
    
    F(:,idx) = A*Ivec;
end

%% Compute U
Cperm = 1;
% Permanent magnet as electro magnets
%U = u*(u0*I1*R^2./(2*(x1vec.^2+R^2).^(3/2)) + u0*I2*R^2./(2*((x1vec-L).^2+R^2).^(3/2)) + Cperm/(((x1-x2)^2 + R^2).^(3/2)));

% Permanent magnets using classical mechanics and integrating
U = u*(-u0*I1*R^2./(2*(x1vec.^2+R^2).^(3/2)) - u0*I2*R^2./(2*((x1vec-L).^2+R^2).^(3/2)) + u0 * u./(4*pi.*(x1vec-x2)));

% stable
%U = u*(-u0*I1*R^2./(2*(x1vec.^2+R^2).^(3/2)) - u0*I2*R^2./(2*((x1vec-L).^2+R^2).^(3/2)) - u0 * u./(4*pi.*(x1vec-x2)));
figure
subplot(2,1,1)
plot(x1vec,U)
title('Potential Function')
ylabel('U')
xlabel('x_1')
subplot(2,1,2)
plot(x1vec,F(1,:),'-b');
title('Force')
ylabel('F')
xlabel('x_1')