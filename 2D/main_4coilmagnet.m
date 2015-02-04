%% main_4coilmagnet
% Start with 1 magnet

close all
clear all

% initialize coil posn, magnet posn
v = zeros(2,1);
x = zeros(2,1);

x(1) = 0.02;
x(2) = -0.02;

xycoil = 0.11;
R = 0.028575;

coords = [-xycoil 0;0 xycoil;xycoil 0;0 -xycoil];
rotation = [pi/2 pi/2 -pi/2;...
            0 pi/2 -pi/2;...
            pi/2 -pi/2 pi/2;...
            0 -pi/2 pi/2]; 

% define current seq
I = [0 512 512 0];
coil = struct('R',R,'current',I,'coords',coords,'rot',rotation);

% xmin xmax ymin ymax = -.1 .1 -.1 .1
params.u0 = 4*pi*10^-7;
params.Cd = 1;

%mass = .0707e-3 * 4; %0.0707g per magnet, 4 magnets, 0.01 default
magnet = struct('Area',0.001,'u',0.004,'m',2.377,'mass',.0707e-3 * 4);

freq = .25; %update freq
endtime = 8;


% loop through current updates
%for i = 1:size(I,1)
tvector = 0;
Xmatrix = [x;v]'; %[x1 x2 x1' x2']
avector = [0 0];
Ihistory = zeros(1,4);
Ihistory(1,:) = [0 0 0 0];

for idx = 1:4
    % upate last state
    xlast = Xmatrix(end,:)';
    
    % update current
    coil.current = zeros(1,4);
    coilon = rem(idx,4); 
    if coilon == 0
        coilon = 4;
    end
    coil.current(coilon) = 512
    Ihistory(idx,:) = coil.current;
    tspan = [idx - 1, idx]/freq;
   
    [tvectortemp, Xmatrixtemp] = ode45(@(t,state) compute_magnet_deriv_4coil(t,state,magnet,coil,params), tspan, xlast);
    
    % update plot
    tvector = [tvector;tvectortemp(2:end)];
    Xmatrix = [Xmatrix;Xmatrixtemp(2:end,:)];    
end
%% 
colorvec = zeros(length(tvector),3);
colorvec(:,2) = linspace(0,1,length(tvector));

figure(1)
scatter(Xmatrix(:,1),Xmatrix(:,2),25,colorvec)
axis([-1 1 -1 1] * xycoil)
for i = 1:idx
    tidx = find(tvector<i/freq,1,'last')
    hold on 
    plot(Xmatrix(tidx,1),Xmatrix(tidx,2),'xb','MarkerSize',15)
end
%axis equal
figure(2)
plot(tvector,Xmatrix(:,1))  
%end




% Plot magnet posn and orientation
