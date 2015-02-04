function F = computeF_4coil(coils,magnet,params,x,v)

% INPUTS:
%  - struct for current configuration of coils
%  - magnets 

coords = coils.coords;
rotation = coils.rot;
I = coils.current;
R = coils.R;

Cd = params.Cd;
u0 = params.u0;

Area = magnet.Area;
u = magnet.u;
m = magnet.m;



%% Compute drag force:
rho = 500; % density of fluid
F_drag = (-1)*sign(v).* 1/2 .* rho .* v.^2 .* Cd .* Area;
%F_drag = 0;

%% define xmesh ymesh 1 dx about magnet posn
dx = .0001;
xlin = [x(1)-dx x(1) x(1)+dx];
ylin = [x(2)-dx x(2) x(2)+dx];
[xmesh,ymesh] = meshgrid(xlin,ylin);

Btot = zeros(size(xmesh,1),size(xmesh,1),2);
% Compute Bfield
for i = 1:4
    onecoil = struct('R',R,'current',I(i),'coords',coords(i,:),'rot',rotation(i,:));
    %compute Bfield
    Btemp = computeBfield(onecoil,xmesh,ymesh);
    %BxCoil(:,:,i) = Btemp(:,:,1);
    %ByCoil(:,:,i) = Btemp(:,:,2);
    Btot(:,:,1) = Btot(:,:,1) + Btemp(:,:,1);
    Btot(:,:,2) = Btot(:,:,2) + Btemp(:,:,2);
end
Bx = Btot(:,:,1);
By = Btot(:,:,2);   

% Compute gradient for F_coil
[dBxdx,~] = gradient(Bx,dx);
[~,dBydy] = gradient(By,dx);

F_coil = [dBxdx(2,2);dBydy(2,2)].*[Bx(2,2);By(2,2)] * u/sqrt(Bx(2,2)^2 + By(2,2)^2);

% compute inter magnet force:
F_mag = 0;

F = F_coil + F_mag + F_drag; % [n x 1]