function [] = resizeFig(width, height)
% See http://tipstrickshowtos.blogspot.com/2010/08/how-to-get-rid-of-white-margin-in.html

% Make the figure boundaries tight
ti = get(gca, 'TightInset');
set(gca, 'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

% Adjust the paper size
set(gca, 'Units','inches');
pos = get(gca, 'Position');
pos = [pos(1) pos(2) width height];
set(gca, 'Position',pos);
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

end