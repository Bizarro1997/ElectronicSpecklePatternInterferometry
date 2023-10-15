function createfigure()

dPic = ([54,119,126,135,146,160,165,171,174,178]-15)*8.9/1280;
err = 7*ones(size(dPic))*8.9/1280;
dX = [];
for p=65:-2:47
d = correlation2D(67,p);
dX = [dX,d];
end
dIP=(1279-dX )*8.9/1280;

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create errorbar
errorbar(dPic,err,'DisplayName','reference scale','MarkerSize',4,'Marker','o',...
    'LineStyle','none',...
    'LineWidth',1.5);

% Create plot
plot(dIP,'DisplayName','maximum correlation','MarkerSize',5,'Marker','+',...
    'LineWidth',2,...
    'LineStyle','none');

% Create ylabel
ylabel({'horizontal displacement / [mm]'});

% Create xlabel
xlabel({'picture number'});

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.5 10.5]);
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'XGrid','on','YGrid','on');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.584699939541922 0.159807292594083 0.290051027185368 0.0841836759717402]);

