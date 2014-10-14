function [I1,I2] = solveforcurrent(x1,x2)

% solve for current through 2 coils of a 1D electromagnet system with 2
% permanent magnet inbetween. 

global R L u

if isempty(R)
    R = 0.01; % radius of coil, 1cm [m]
    L = .1; % separation between 2 coils [m];
    u = 0.05; %Am^2
end

u0 = 4*pi*10^-7; % Tm/A

C(1) = -3*u*x1/(R^2+x1^2)^(5/2)*u0/2*R^2;
C(2) = -3*(x1-L)*u/(R^2 + (x1-L)^2)^(5/2)*u0/2*R^2;
C(3) = -u0*u*u/(4*pi*(x1-x2)^2);
C(4) = -3*u*(x1-L)/(R^2+(x1-L)^2)^(5/2)*u0/2*R^2;
C(5) = -3*(x2-L)*u/(R^2 + (x2-L)^2)^(5/2)*u0/2*R^2;
C(6) = u0*u*u/(4*pi*(x1-x2)^2);

% Ax = B

B = [C(3);C(6)];
A = [C(1) C(2);C(4) C(5)];

x = A\B;
I1 = x(1);
I2 = x(2);