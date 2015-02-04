cmap = ametrine(2);
c1 = cmap(1,:);
c2 = cmap(2,:);

figure
%hf = quiver(x,y,Fx./Fnorm,Fy./Fnorm,.5);
%set(hf,'Color','r')
% plot streamlines
[xstartmesh,ystartmesh] = meshgrid(3,linspace(-3,3,11));
hstream = streamline(x,y,Fx,Fy,xstartmesh,ystartmesh);
set(hstream,'Color',c1)
hold on
[xstartmesh,ystartmesh] = meshgrid(linspace(-3,3,11),[-3 3]);
hstream = streamline(x,y,Fx,Fy,xstartmesh,ystartmesh);
set(hstream,'Color',c1)



hold on
hm = quiver(x,y,Bx./normB,By./normB,.5);
set(hm,'Color',c2)
legend({'Force','Orientation'})
axis equal
axis([-xmax,xmax,ymin,ymax])
xlabel('x')
ylabel('y')
