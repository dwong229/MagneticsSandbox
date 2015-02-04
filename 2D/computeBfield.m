function [B] = computeBfield(coil,worldx,worldy)
% 
% INPUTS: 
% 
% given coil parameters as a struct, compute 

R = coil.R;
I = coil.current;
mu0 = 4*pi*10^-7;
worldcoil = [coil.coords';0];
th = coil.rot;
%Rmat = Rz Rx Rz;
Rz = [cos(th(1)) sin(th(1)) 0;-sin(th(1)) cos(th(1)) 0;0 0 1];
Rx = [1 0 0;0 cos(th(2)) sin(th(2));0 -sin(th(2)) cos(th(2)) ];
Rz2 = [cos(th(3)) sin(th(3)) 0;-sin(th(3)) cos(th(3)) 0;0 0 1];
Rmat = Rz2*Rx*Rz;

worldxy = [worldx(:)';worldy(:)';zeros(1,length(worldx(:)))];
localcoord = Rmat*bsxfun(@plus,worldxy,-worldcoil);

localz = localcoord(3,:);
localy = localcoord(2,:);
%ylin = localcoord(2,:);
%zlin = localcoord(3,:);

% compute in local frame
%[localz,localy] = meshgrid(zlin,ylin);

kz2 = 4 * R * abs(localy)./((R + abs(localy)).^2 + localz.^2);
Bz = mu0 * I./(2*pi) * 1./(sqrt((R + abs(localy)).^2 + localz.^2)) .* (ellipticF(pi/2,kz2) + (R^2 - localy.^2 - localz.^2)./((R-abs(localy)).^2 + localz.^2).* ellipticE(pi/2,kz2));

ky2 = -4 * R * localy./((R-localy).^2 + localz.^2);
sqrtRYZ = sqrt(R^2 + localy.^2 + localz.^2);
By = mu0 * I * R .* localz ./(4 * pi).* (1./(R.*localy.*sqrtRYZ.*((R+localy).^2 + localz.^2)) .* (sqrtRYZ.^2 .* (ellipticE(pi/4,ky2) + ellipticE(3*pi/4,ky2)) - (sqrtRYZ.^2 + 2.*R.*localy ).*(ellipticF(pi/4,ky2) + ellipticF(3*pi/4,ky2)))).* sqrtRYZ ./ sqrt((R-localy).^2 + localz.^2);
%if sum(isnan(By(:))) > 0
%    By(localy==0) = 0;
%end
%keyboard

%% map local frame into world coords
B = zeros(size(worldx,1),size(worldx,2),2);
%populate Bx
B(:,:,1) = reshape(Rmat(2,1)*By + Rmat(3,1)*Bz,size(worldx));
B(:,:,2) = reshape(Rmat(2,2)*By + Rmat(3,2)*Bz,size(worldx));

%localcoordgrid = [zeros(1,length(localz(:)));localy(:)';localz(:)'];
%worldcoordgrid = bsxfun(@plus,Rmat'*localcoordgrid,worldcoil);

%B(worldcoordgrid(1,:),worldcoordgrid(2,:),1) = By(localz(:),localy(:));
%B(worldcoordgrid(1,:),worldcoordgrid(2,:),2) = Bz(localz(:),localy(:));


end
