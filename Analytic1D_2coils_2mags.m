%% Analytically Verify 1D controllability
%  

close all
clear all
%% Visualize point dipole model %%
% KEY
% 

R = 0.01; % radius of coil, 1cm [m]
L = .1; % separation between 2 coils [m];
dx = 0.001; % resolution, 1mm

%% 
% idx = 1;
% for x1 = 0:dx:L-dx
%     for x2 = x1+dx:dx:L
%         x1vec(idx) = x1;
%         x2vec(idx) = x2;
%         z(idx) = (x1*x2-x1*L)/(R^4 + R^2*(x2-L)^2 + x1^2*(x2-L)^2 + x1^2*R^2)^(5/2) - (x1*x2 - x2*L)/(R^4 + R*x2^2 + x2^2*(x1-L)^2 + (x1-L)^2*R^2)^(5/2);
%         idx = idx + 1;
%     end
% end
% 
% plot3(x1vec,x2vec,z,'.b')

[x1,x2] = meshgrid([0:dx:L]);
z = (x1.*x2-x1*L)./((R^4 + R^2.*(x2-L).^2 + x1.^2.*(x2-L).^2 + x1.^2.*R.^2).^(5/2)) - (x1.*x2 - x2*L)./((R.^4 + R*2.*x2.^2 + x2.^2.*(x1-L).^2 + (x1-L).^2*R.^2).^(5/2));

% filter data
invalididx = find(x1(:)>=x2(:));
ij_invalididx = ind2sub(size(x1),invalididx);
x1(invalididx) = 0;
x2(invalididx) = 0;
z(invalididx) = NaN;

surf(x1,x2,z)
colorbar

hold on
% draw zero plane
plane(:,1) = [0,0,0];
plane(:,2) = [0,L,0];
plane(:,3) = [L,L,0];
plane(:,4) = [L,0,0];

fill3(plane(1,:),plane(2,:),plane(3,:),'b');
alpha(0.3);

%plot3(x1(z(:) == 0),x2(z(:) == 0),z(z(:) == 0),'ow')

title('det of matrix')
xlabel('x_1')
ylabel('x_2')
zlabel('det(A)')


figure
plot3(x1(z>0),x2(z>0),z(z>0),'.b')
hold on
plot3(x1(z<0),x2(z<0),z(z<0),'.r')
title('det of matrix')
xlabel('x_1')
ylabel('x_2')
zlabel('det(A)')
legend(['det>0';'det<0'])
%%%%
%% Given x1 and x2 solve for I1 and I2

inputx1 = 0.04;
inputx2 = 0.06;

[I1,I2] = solveforcurrent(inputx1,inputx2); 

fprintf('For a coil separation of L = %2.2f m with coils of radius R = %2.2f \n',L,R)
fprintf('To hold 2 magnets at x = %3.2f and %3.2f m\n',inputx1,inputx2)
fprintf('Current needed: I1 = %7.4f A, I2 = %7.4f A \n',I1,I2)

%% Plot for I1 and I2
idx = 1;
buffer = 3;
for x1 = buffer*dx:dx:L-buffer*dx
    for x2 = x1+10*dx:dx:L-buffer*dx
        x1vec(idx) = x1;
        x2vec(idx) = x2;
        [I1(idx),I2(idx)] = solveforcurrent(x1vec(idx),x2vec(idx)); 
        idx = idx + 1;
    end
end

figure
plot3(x1vec,x2vec,I1,'.b')
hold on
plot3(x1vec,x2vec,I2,'.r')
xlabel('x1')
ylabel('x2')
zlabel('Current')
legend(['I_1';'I_2'])