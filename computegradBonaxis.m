function dB = computegradBonaxis(x,I,R,x0)

% compute gradient of magnetic field
% INPUTS:
% x: vector of positions along axis 
% I: current around coil
% R: radius of coil
% x0: center of coil in world coords along x
%
% OUTPUT:
% dB: vector of gradient of B same length as x

% constants
mu0 = 4*pi*10^-7; 

% position in coil coordinates
dx = x - x0; % compute distance from cetner of coil

dB = zeros(1,length(x));

dB = 3 * dx./((R^2 + dx.^2).^(5/2)) * (mu0/2*I*R^2);




