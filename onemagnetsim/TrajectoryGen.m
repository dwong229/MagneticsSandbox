% %% jtraj sandbox
% close all
% xstart = 0;
% xend = 10;
% vzero = 0;
% 
% [Q,QD,QDD] = jtraj(xstart,xend,100,vzero,vzero);
% 
% hold all
% plot(Q)
% plot(QD)
% plot(QDD)
% legend({'Q';'QD';'QDD'})

%% 

L = 10;
x1 = linspace(0,L,100);
x2 = x1;
[x,y] = meshgrid(0:L,0:L);
%for i = 1:length(x1)
    %for j = 1:length(x2)
        %A = [x1(i) x1(i)-L;x2 x2(i)-L];
        %lambda = eig(A);
        e1 = 1/2*(x + y - L - sqrt(L^2 + 2*L*x - 6*L*y+x.^2 + 2*x.*y+y.^2));
        e2 = 1/2*(x + y - L + sqrt(L^2 + 2*L*x - 6*L*y+x.^2 + 2*x.*y+y.^2));
        
        xlessthany = flipud(triu(ones(size(e1)),1));
%    end
%end


% subplot(2,1,1)
% xlabel('x1')
% ylabel('x2')
% plot3(x(xlessthany==1),y(xlessthany==1),real(e1(xlessthany==1)),'.')
% subplot(2,1,2)
% plot3(x(xlessthany==1),y(xlessthany==1),imag(e1(xlessthany==1)),'.')
% xlabel('x1')
% ylabel('x2')

subplot(2,2,1)
mesh(x,y,real(e1))
title('Re(e1)')
xlabel('x1')
ylabel('x2')
subplot(2,2,2)
mesh(x,y,imag(e1))
title('Im(e1)')
xlabel('x1')
ylabel('x2')

subplot(2,2,3)
mesh(x,y,real(e2))
title('Re(e2)')
xlabel('x1')
ylabel('x2')
subplot(2,2,4)
mesh(x,y,imag(e2))
title('Im(e2)')
xlabel('x1')
ylabel('x2')