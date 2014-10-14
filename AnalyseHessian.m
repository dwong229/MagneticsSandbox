%% Solve for eigenvalues of Hessian Matrix

global R L u %for solveforcurrent.m
 

%variables to set up system
u0 = 4*pi*10^(-7); %Tm/A
u1 = 1; %Am^2 % magnetic moment
u2 = 1; %Am^2 % magnetic moment
R1 = .01;
R2 = .01;
L = .1;
m1 = 1; %Am % magnetic magnitude
m2 = 1; %Am % magnetic magnitude

R = 0.01; % radius of coil, 1cm [m]
L = .1; % separation between 2 coils [m];
u = 0.05; %Am^2

% Plug in stable configuration
I1 = 1; 
I2 = 1;

x1 = 0.1*L;
x2 = 0.3*L;

[I1,I2] = solveforcurrent(x1,x2);


% Hessian for magnet 1
H(1,1) = -3*u0*u1*I1*R2^2/(2*((x1-L)^2+R2^2)^(5/2)) + 15*(x1-L)^2*u0*u1*I1*R2^2/(2*((x1-L)^2+R2^2)^(7/2)) + 15*x1^2*u0*u1*I1*R1^2/(2*(R1^2+x1^2)^(7/2)) - 3*u0*u1*I1*R1^2/(2*(R1^2+x1^2)^(5/2)) + u0*m1*m2/(2*pi*(x2-x1)^3);
H(1,2) = -u0*m1*m2/(2*pi*(x2-x1)^3);
H(2,1) = -u0*m1*m2/(2*pi*(x2-x1)^3);
H(2,2) = u0*m1*m2/(2*pi*(x2-x1)^3);

% Solve for eigenvalues
[v,d] = eig(H)
