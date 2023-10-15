function createfigure()

dX = [];
for p=123:-1:119
d = plot3D(p);
dX = [dX,d];
end

err = 0.543/4*ones(size(dX))/2;

% Create figure
close all;
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create errorbar
yyaxis left
errorbar((-1:-1:-5),dX,err,'MarkerSize',4,'Marker','o',...
    'LineStyle','none',...
    'LineWidth',1.5);

ylabel({['relative peak / [{\mu}m]']});

yyaxis right

err = 0.543/4*ones(size(diff(dX)))/2*sqrt(2)*2.5;
errorbar((-1.5:-1:-4.5),diff(dX)*2.5,err,'MarkerSize',4,'Marker','o',...
    'LineStyle','none',...
    'LineWidth',1.5);
ylabel({['peak speed / [{\mu}m/s]']});
% Create ylabel

% Create xlabel
xlabel({'relative picture number'});


% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.5 10.5]);
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'XGrid','on','YGrid','on', 'XTick', [-5:1:-1]);
% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,...
    %'Position',[0.591842796684779 0.806065795995443 0.290051027185368 0.0841836759717401]);

