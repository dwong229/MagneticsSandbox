function statedot = compute_magnet_derivatives(t,state)

% unpack states
x = state(1);
v = state(2);

% pack systemdef
global rho magnet coil
Cd = magnet.Cd;
A = magnet.A;
u1 = magnet.u1;
u0 = coil.u0;
I = coil.I;
Rcoil = coil.R;
xcoil = coil.x;
mass = magnet.mass;

F_drag = (-1)*sign(v)* 1/2* rho * v^2 * Cd * A;
% Force calculations
p = xcoil - x; % vector from magnet to the coil
F_B = u1 * 3 * p * u0 * I * Rcoil^2/(2*(p^2 + Rcoil^2)^(5/2));

% compute acceleration
F_total = F_drag + F_B;
a = F_total/mass;

% pack state derivatives
statedot = [v a]';