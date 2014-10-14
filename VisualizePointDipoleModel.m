close all
clear all
%% Visualize point dipole model %%
% KEY
% h : magnetic field [3x1] exists in R3
% p : p_c - p_a
% p_c : center of actuator magnet
% p_a : center of actuator
% m_a : dipole moment of permanent magnet

% B : on-axis magnetic field from a current loop

m_a = [0 0 1]'; % saturization magnetization per unit volume of material
v = 1; % magnetic volume of material
mu_0 = 4*pi*10^-7; % [Tm/A]
I = 10; % current [A] I = V/R
R = 0.01; % radius of coil, 1cm [m]

z = 0:.0001:0.5; % distance from center for coil [m]

%% On axis magnetic field of a current loop

for i = 1:length(z)
    B(i) = mu_0 * I *R^2/(2*(z(i)^2 + R^2)^(3/2));
    Fz(i) = -v*m_a(3)*mu_0 * I *R^2/2*3*z(i)/((z(i)^2 + R^2)^(5/2));
end

figure('Position',[67 230 766 716]);
subplot(2,1,1)
plot(z,B,'.-b','LineWidth',2)
xlabel('Distance from center for the coil [m]')
ylabel('Magnetic Field Strength [T]')
title('Magnetic field in direction of the axis along the axis of the coil')
subplot(2,1,2)
plot(z,Fz,'.r','LineWidth',2)
xlabel('Distance from center for the coil [m]')
ylabel('Force')
title('Magnetic field in direction of the axis along the axis of the coil')