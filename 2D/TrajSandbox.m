v = zeros(4,1);
x = zeros(4,1);

x(1) = 0.06;
x(2) = 0;


xycoil = 0.11;
R = 0.028575;

coords = [-xycoil 0;0 xycoil;xycoil 0;0 -xycoil];
rotation = [pi/2 pi/2 -pi/2;...
            0 pi/2 -pi/2;...
            pi/2 -pi/2 pi/2;...
            0 -pi/2 pi/2]; 
% define current seq
I = [0 0 250 250];
coil = struct('R',R,'current',I,'coords',coords,'rot',rotation);


%mass = .0707e-3 * 4; %0.0707g per magnet, 4 magnets, 0.01 default
magnet = struct('Area',0.001,'u',0.004,'m',2.377,'mass',0707e-3 * 4);

% xmin xmax ymin ymax = -.1 .1 -.1 .1
params.u0 = 4*pi*10^-7;
params.Cd = 1;
F = computeF_4coil(coil,magnet,params,x,v)
