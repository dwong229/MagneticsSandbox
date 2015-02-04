%% Visualize 2D magnetic flux
close all
clear all
clc

mu0 = 4*pi * 10^-7;
n = 1;
I = .5;
R = .001;

ymin = -.0025;
ymax = .0025;
zmin = -.0025;
zmax = .0025;
numpts = (ymax-ymin)/0.0005;

xycoil = 3*R;
%% Plot options:
plotB = true;
plotF = false;

%% Don't change below here
ylin = linspace(ymin,ymax,numpts); % default - 100 pts
zlin = linspace(zmin,zmax,numpts);

yres = diff(ylin(1:2));
zres = diff(zlin(1:2));

[z,y] = meshgrid(zlin,ylin);

kz2 = 4 * R * abs(y)./((R + abs(y)).^2 + z.^2);
Bz = mu0 * I./(2*pi) * 1./(sqrt((R + abs(y)).^2 + z.^2)) .* (ellipticF(pi/2,kz2) + (R^2 - y.^2 - z.^2)./((R-abs(y)).^2 + z.^2).* ellipticE(pi/2,kz2));

%contour3(z,y,Bz)
surfc(z,y,Bz)
%drawcoil(R,gcf)
xlabel('z')
ylabel('y')
zlabel('Bz')

%% dBz/Bz

%dB = diff(Bz,1,2);
%figure
%surfc(dB)
%figure

%% By
%R = 0.001;

%Byfun = @(x) sin(x)./(R^2 + y.^2 + z.^2 - 2.*y.*R.*sin(x)).^(3/2);
%By_int = mu0 * I * R .* z ./(4 * pi) .* integral(@(x) Byfun(x),0,2*pi,'ArrayValued',true);

% exact solution - not quite there
ky2 = -4 * R * y./((R-y).^2 + z.^2);
sqrtRYZ = sqrt(R^2 + y.^2 + z.^2);
By = mu0 * I * R .* z ./(4 * pi).* (1./(R.*y.*sqrtRYZ.*((R+y).^2 + z.^2)) .* (sqrtRYZ.^2 .* (ellipticE(pi/4,ky2) + ellipticE(3*pi/4,ky2)) - (sqrtRYZ.^2 + 2.*R.*y ).*(ellipticF(pi/4,ky2) + ellipticF(3*pi/4,ky2)))).* sqrtRYZ ./ sqrt((R-y).^2 + z.^2);
if sum(isnan(By(:))) > 0
    By(y==0) = 0;
end

figure
surfc(z,y,By)
title('By')
xlabel('z')
ylabel('y')
zlabel('By')

% subplot(1,2,1)
% surfc(z,y,By_int)
% title('By Integral')
% xlabel('z')
% ylabel('y')
% zlabel('By')
% subplot(1,2,2)
% surfc(z,y,By_ell)
% title('By Ellipse')
% xlabel('z')
% ylabel('y')
% zlabel('By')
%% R << (y^2 + z^2)^(1/2)
% R = 0.001;
% By_smallRadius = mu0*I/(4*pi)*3*pi*R^2.*y.*z.*(y.^2 + z.^2).^(-5/2);
% figure
% surfc(z,y,By_smallRadius)
% xlabel('z')
% ylabel('y')
% zlabel('By_smallradius')

%% Bmagnitude and Vector field plot
if plotB 
figure('Position',[114 546 1120 420])
subplot(1,2,1)
Btot = sqrt(By.^2 + Bz.^2);
title('Magnitude of B')
surfc(z,y,Btot)
xlabel('z')
ylabel('y')
zlabel('Btot')

subplot(1,2,2)
%quiver(z,y,Bz,By,2);

% normalized quiver
normB = sqrt(Bz.^2 + By.^2);
quiver(z,y,Bz./normB,By./normB,.5);

n = 10;
ystart = [linspace(ymin,ymax,n),linspace(ymin,ymax,n),ymin*ones(1,n),ymax*ones(1,n)];
zstart = [zmin*ones(1,n),zmax*ones(1,n),linspace(zmin,zmax,n),linspace(zmin,zmax,n)];
hstream = streamline(z,y,Bz,By,zstart,ystart);
set(hstream,'Color','r')
hold on 
contour(z,y,Btot,500)
axis equal
axis([0,zmax,ymin,ymax])
xlabel('z')
ylabel('y')
end

%% Force plot
if plotF
% assuming maximum force: 
uvec = zeros(size(Bz,1),size(Bz,2),2); %[z,y,{k,j}]
uvec(:,:,1) = Bz./normB;
uvec(:,:,2) = By./normB;

[Btotdz, Btotdy] = gradient(Btot,abs(z(1,1)-z(1,2)),abs(y(1,1)-y(2,1)));
Fz = Btotdz.*uvec(:,:,1)
Fy = Btotdy.*uvec(:,:,2);
Fnorm = sqrt(Fz.^2 + Fy.^2);

figure('Position',[114 546 1120 420])
subplot(1,2,1)
surfc(z,y,Fnorm)
xlabel('z')
ylabel('y')
zlabel('Magnitude of Force')

subplot(1,2,2)
hf = quiver(z,y,-Fz./Fnorm,-Fy./Fnorm,.5);
set(hf,'Color','r')
hold on
hm = quiver(z,y,Bz./normB,By./normB,.5);
set(hm,'Color','b')
legend({'Force','Orientation'})
axis equal
axis([0,zmax,ymin,ymax])
xlabel('z')
ylabel('y')
end



%% Gradient
% if false
% [Btotdz, Btotdy] = gradient(Btot,abs(z(1,1)-z(1,2)),abs(y(1,1)-y(2,1)));
% 
% figure
% subplot(2,1,1)
% surfc(Btotdy)
% subplot(2,1,2)
% surfc(Btotdz)
% figure
% quiver(z,y,Btotdz,Btotdy,5)
% ystart = [ymin:0.1:ymax,ymin:0.1:ymax];
% zstart = [zmin*ones(1,length(ystart)/2),zmax*ones(1,length(ystart)/2)];
% hstream = streamline(z,y,Btotdz,Btotdy,zstart,ystart)
% set(hstream,'Color','r')
% 
% gradB = sqrt(Btotdz.^2 + Btotdy.^2);
% figure
% contour(z,y,gradB,500)
% end




