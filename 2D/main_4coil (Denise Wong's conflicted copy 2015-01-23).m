% Plot 2D field for 4 coil configuration
close all
clear all

%% Define parameters for 4 coils
R = 0.5;

mu0 = 4*pi * 10^-7;
I = [1 1 1 1];

ymin = -3;
ymax = 3;
xmin = -3;
xmax = 3;
numpts = 30;
xycoil = 10*R;
% 4coil_plane_dwedit params
R = 0.001; 
ymin = -0.0024;
ymax = 0.0024;
xmin = -0.0024;
xmax = 0.0024;
numpts = round((ymax-ymin)/0.0001+1);
xycoil = 0.0025;


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
BxCoil = zeros(numpts,numpts,4);
ByCoil = zeros(numpts,numpts,4);
for i = 1:4
    %i
    % pack variables
    coil = struct('R',R,'current',I(i),'coords',coords(i,:),'rot',rotation(i,:));
    % compute Bfield
    Btemp = computeBfield(coil,x,y);
    BxCoil(:,:,i) = Btemp(:,:,1);
    ByCoil(:,:,i) = Btemp(:,:,2);
    Btot(:,:,1) = Btot(:,:,1) + Btemp(:,:,1);
    Btot(:,:,2) = Btot(:,:,2) + Btemp(:,:,2);
end
%keyboard
%% Plot Bfield
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
%quiver(x,y,Bx./normB,By./normB,.5);
quiver(x,y,Bx,By);


n = 20;
%ystart = [linspace(ymin,ymax,n),linspace(ymin,ymax,n),ymin*ones(1,n),ymax*ones(1,n)];
%xstart = [xmin*ones(1,n),xmax*ones(1,n),linspace(xmin,xmax,n),linspace(xmin,xmax,n)];
ystart = [linspace(-2,2,n),linspace(-2,2,n),-2*ones(1,n),2*ones(1,n)];
xstart = [-2*ones(1,n),2*ones(1,n),linspace(-2,2,n),linspace(-2,2,n)];
hstream = streamline(x,y,Bx,By,xstart,ystart);
set(hstream,'Color','c')
hold on 
contour(x,y,Bsum,n)
%axis equal
axis([-xmax,xmax,ymin,ymax])
axis square
xlabel('x')
ylabel('y')

%% Calculate and Plot force
% calculate direction of magnet
uvec = zeros(size(Bx,1),size(Bx,2),2); %[z,y,{k,j}]
uvec(:,:,1) = Bx./normB;
uvec(:,:,2) = By./normB;

%[ddxBxBx,~] = gradient(Bx.^2,xres,yres);
%[~,ddxByBy] = gradient(By.^2,xres,yres);

%Fx = ddxBxBx.*uvec(:,:,1)
%Fy = Btotdy.*uvec(:,:,2);

[Btotdx, ~] = gradient(Bx.^2,abs(x(1,1)-x(1,2)),abs(y(1,1)-y(2,1)));
[~, Btotdy] = gradient(By.^2,abs(x(1,1)-x(1,2)),abs(y(1,1)-y(2,1)));
Fx = Btotdx; %.*uvec(:,:,1);
Fy = Btotdy; %.*uvec(:,:,2);
Fnorm = sqrt(Fx.^2 + Fy.^2);

%% Force Plots
figure('Position',[167 572 1451 402])
subplot(1,3,1)
surfc(x,y,Fnorm)
xlabel('x')
ylabel('y')
zlabel('Magnitude of Force')
title('Force Magnitude Spatial Distribution')


subplot(1,3,2)
hf = quiver(x,y,Fx./Fnorm,Fy./Fnorm,.5);
set(hf,'Color','r')
hold on
hm = quiver(x,y,Bx./normB,By./normB,.5);
set(hm,'Color','b')
legend({'Force','Orientation'})
axis equal
axis([-xmax,xmax,ymin,ymax])
xlabel('x')
ylabel('y')
title('Force Direction and Magnet Orientation Plot')

% Force quiver plot
subplot(1,3,3)
hquiv = quiver(x,y,Fx,Fy,3);
set(hquiv,'Color','r')
hold on
[xstartmesh,ystartmesh] = meshgrid(linspace(-2,2,10),linspace(-2,2,10));
plot(xstartmesh,ystartmesh,'x','MarkerFaceColor',[0.5,0,0]);
hstream = streamline(x,y,Fx,Fy,xstartmesh,ystartmesh);
set(hstream,'Color',[0.5,0,0])
hold on 
contour(x,y,Bsum,10)
axis equal
axis([-xmax,xmax,ymin,ymax])
title('Force Vector Field')
xlabel('x')
ylabel('y')
%% Save figure
str = input('Save file? Y/N [Y]:','s');
if isempty(str)
    str = 'Y';
end
if str == 'Y'
    set(gcf,'PaperPositionMode','auto')
    filename = strcat('./figures/Current_',num2str(I(1)),'_',num2str(I(2)),'_',num2str(I(3)),'_',num2str(I(4)),'_Force');
    
    disp('Saving file:')
    disp(filename)
    print('-depsc','-tiff','-r0',filename)
else
    disp('Not Saving file.')
end
    



