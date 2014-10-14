function [varargout] = drawmagnet(centerx,centery,w,h,mdirection,varargin)

%%%%% Draw a barmagnet in white (S) and red (N) where mdirection is the
%%%%% vector point from S to N

% coord of bottom left corner
cornerx = centerx-w/2;
cornery = centery-h/2;

if sign(mdirection) == 1
    facecolor = ['w','r'];
else
    facecolor = ['r','w'];
end

if nargin == 5
    %disp('Initialize magnet handle')
        
    hmagnet(1) = rectangle('Position',[cornerx,cornery,w/2,h],'LineWidth',2,'FaceColor',facecolor(1));
    hmagnet(2) = rectangle('Position',[centerx,cornery,w/2,h],'LineWidth',2,'FaceColor',facecolor(2));
    varargout = {hmagnet};
else
    %disp('Update magnet handle')
    hmagnet = varargin{1};
    posn1 = get(hmagnet(1),'Position');
    posn1(1:2) = [cornerx,cornery];
    set(hmagnet(1),'Position',posn1);
    
    posn2 = get(hmagnet(2),'Position');
    posn2(1:2) = [centerx,cornery];
    set(hmagnet(2),'Position',posn2);
    
    drawnow
    
    varargout = [];
    
end


