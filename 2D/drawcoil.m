function [] = drawcoil(R,hfig)
    figure(hfig)
    hold on
    [x,y,z] = sphere(20);
    surf(zeros(size(y)),y,z,'FaceColor','none','MeshStyle','column','EdgeColor','k');
end