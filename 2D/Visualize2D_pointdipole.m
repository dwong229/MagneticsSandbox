% Point Dipole Model

% R << (y^2 + z^2)^(1/2)

close all
clear all
clc

mu0 = 4*pi * 10^-7;
n = 1;
I = 1;
R = .001;
uy = 1;
uz = 1;

ymin = -2;
ymax = 2;
zmin = -2;
zmax = 2;

ylin = linspace(ymin,ymax,51); % default - 100 pts
zlin = linspace(zmin,zmax,51);

[z,y] = meshgrid(zlin,ylin);
%y(abs(y)<R*2) = 0;
%z(abs(z)<R*2) = 0;
%Fy = (3 * I *R^2 * (uy * z .*(-4 * y.^2 + z.^2) + uz * (y.^3 - 4 * y.* z.^2)))./(4 * (y.^2 + z.^2).^(7/2));
%Fz = (3 * I *R^2 * (uy * y.^3 + 3 * uz * y.^2 .* z - 4*uy*y.*z.^2 - 2*uz *z.^3))./(4 * (y.^2 + z.^2).^(7/2));
%Ftot = 3/4* sqrt((I^2 * R^4 * mu0^2*(-2 * uy * uz * z .*(y.^3 - 4*y.*z.^2) + uy^2*(y.^4 + 7 * y.^2 .* z.^2 + z.^4) + uz^2* (y.^4 + 4 * z.^4)))./(y.^2 + z.^2).^6);
%Ftot = sqrt(Fy.^2 + Fz.^2);

By = 3 * I * R^2 * mu0 * y .*z ./(4 * (y.^2 + z.^2).^(5/2));
Bz = I * R^2 * mu0 .* ( 2 - 3*y.^2 ./(y.^2 + z.^2))./(4 * (y.^2 + z.^2).^(3/2));
Btot = sqrt(By.^2 + Bz.^2);

[dBtotdz,dBtotdy] = gradient(Btot);
gradB = sqrt(dBtotdz.^2 + dBtotdy.^2);

%%
figure
subplot(1,2,1)
surfc(z,y,By)
subplot(1,2,2)
surfc(z,y,Bz)

figure

contour(z,y,Btot,5000)
axis equal
figure
contour(z,y,gradB,5000)

