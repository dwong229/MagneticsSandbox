% Plot 2D field for 4 coil configuration

%% Define parameters for 4 coils
R = 0.5;
mu0 = 4*pi * 10^-7;
I = 1;

ymin = -3;
ymax = 3;
xmin = -3;
xmax = 3;
numpts = 15;

xycoil = 8*R;

coords = [-xycoil 0;0 xycoil;xycoil 0;0 -xycoil];
% Euler angles - zxz : [local] = R*[world]
rotation = [pi/2 pi/2 -pi/2;...
            0 pi/2 -pi/2;...
            pi/2 -pi/2 pi/2;...
            0 -pi/2 pi/2]; 

%% Don't change below here
ylin = linspace(ymin,ymax,numpts); % default - 100 pts
xlin = linspace(xmin,xmax,numpts);

yres = diff(ylin(1:2));
xres = diff(xlin(1:2));

[x,y] = meshgrid(xlin,ylin);

% Initialize Bfield parameter holder of 4 coils
Btot = zeros(numpts,numpts,2); %:,:,[x,y]

for i = 1:4
    i
    % pack variables
    coil = struct('R',R,'current',I,'coords',coords(i,:),'rot',rotation(i,:));
    % compute Bfield
    Btemp = computeBfield(coil,x,y);
    Btot(:,:,1) = Btot(:,:,1) + Btemp(:,:,1);
    Btot(:,:,2) = Btot(:,:,2) + Btemp(:,:,2);
end

%%
Bx = Btot(:,:,1);
By = Btot(:,:,2);

figure('Position',[114 546 1120 420])
subplot(1,2,1)
Bsum = sqrt(By.^2 + Bx.^2);
title('Magnitude of B')
surfc(x,y,Bsum)
xlabel('x')
ylabel('y')
zlabel('Bsum')

subplot(1,2,2)
quiver(x,y,Bx,By,2);

% normalized quiver
normB = sqrt(Bx.^2 + By.^2);
quiver(x,y,Bx./normB,By./normB,.5);

n = 10;
ystart = [linspace(ymin,ymax,n),linspace(ymin,ymax,n),ymin*ones(1,n),ymax*ones(1,n)];
xstart = [xmin*ones(1,n),xmax*ones(1,n),linspace(xmin,xmax,n),linspace(xmin,xmax,n)];
hstream = streamline(x,y,Bx,By,xstart,ystart);
set(hstream,'Color','r')
hold on 
contour(x,y,Bsum,500)
%axis equal
axis([-xmax,xmax,ymin,ymax])
xlabel('x')
ylabel('y')
