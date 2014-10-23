function statedot = compute_magnet_derivatives(t,state,magnet,coilarray)

% compute derivative of state vector for 

x = state(1);
v = state(2);

Cd = magnet.Cd;
A = magnet.A;
u1 = magnet.ui;
u0 = coilarray(1).u0;
I = [coilarray(:).I];
Rcoil = [coilarray(:).R];
xcoil = [coilarray(:).x];
mass = magnet.mass;

p = xcoil - x;
rho = 1000;
F_drag = (-1)*sign(v)* 1/2* rho * v^2 * Cd * A;
F_drag = 0;
Fcoil = u1 * 3 * u0 * I .*p.* Rcoil.^2./(2*(p.^2 + Rcoil.^2).^(5/2));

a = (sum(Fcoil)+F_drag)/mass;

statedot = [v,a]';
