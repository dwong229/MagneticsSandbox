function F = computeFelectromagnet(dB,M)

%% computes Force acting on a magnet with volume v and saturization magnetization M
% from the gradient of the magnetic field B using F = v(M dot grad(B))
% 
% INPUTS
% v: volume of magnet (m^3)
% M: saturization magnetization in direction of gradient
% dB: gradient of magnetic field (B) in direction of F

% OUTPUT:
% F: scalar of force in direction of dB and M

F = M*dB;

