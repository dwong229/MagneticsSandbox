% Visualize interaction between 2 magnets

% place one magnets at [0,0] and simulate other magnets at p away.  
% 0 : magnet fixed at origin
% 1 : magnet in space

% Define geometry:
p = [1 1]'; % vector from [0,0] to magnet-c
m0 = [0,1]';
m1 = [0,1]';

% magnetude of magnetic field
m0mag = norm(m0);
m1mag = norm(m1);

m0hat = m0/m0mag;
m1hat = m1/m1mag;

% constants
u0 = 4*pi*10^-7;

F = 3*u0*m0mag*m1mag/(4*pi*norm(p)^2)*(m0hat * m1hat' + m1hat * m0hat' + (m1hat'*(eye(2) - 5*p*p')*m0hat)*eye(2)) * p;

tempp = p;

% compute B-field of 
B = u0/(4*pi * norm(tempp)^3)*(3*tempp/norm(tempp) * tempp'/norm(tempp) - eye(2))*m0hat/m0mag %|m|/(4pi |p|^3) (3hat(p)hat(p)' - I) hat(m)
Bx = u0/(4*pi*norm(p)^3) * ((3*p(1)^2/-1) * m0(1) + 3*p(1)*p(2)*m0(2))/m0mag
By = u0/(4*pi*norm(p)^3) * ((3*p(2)^2-1) * m0(2) + 3*p(1)*p(2)*m0(1))/m0mag
% torque:
D = 3*p/norm(p) * p'/norm(p) - eye(2);

%t = cross([u0 * m0mag * m1mag/(4*pi*norm(p)^3) * m1hat;0],[D*m0hat;0]);
t = det([u0 * m0mag * m1mag/(4*pi*norm(p)^3) * m1hat,D*m0hat]);

%% Visualize interaction forces
plot([0,p(1)],[0,p(2)],'ok')
hold on
quiver([0,p(1)],[0,p(2)],[m0(1),m1(1)],[m0(2),m1(2)],'Color','b')
hold on
quiver(p(1),p(2),F(1),F(2),10e5,'Color','r')
legend({'m','F'})
xlabel('x')
ylabel('y')
axis equal