function scalebar(kwargs)
%scalebar('fig', 'ax', 'im', 'pixelSize', 'scaleBar', 'unit', 'color', 'location')

arguments
    kwargs.fig = 'none';
    kwargs.ax = 'none';
    kwargs.im = 'none';
    kwargs.pixelSize = 4.7e-6
    kwargs.scaleBar = 200 %in micron
    kwargs.unit = 'micron';
    kwargs.color = 'k';
    kwargs.location = 'bottom left';
end

if strcmp(kwargs.fig, 'none')
    if strcmp(kwargs.ax, 'none')
        fig = gcf();
    else
        fig = ancestor(kwargs.ax, {'figure'}, 'toplevel');
    end
else
    fig = kwargs.fig;
end

if strcmp(kwargs.ax, 'none')
    ax = axes('Parent', fig);
else
    ax = kwargs.ax;
end

hold(ax, 'on')

switch kwargs.unit
    case 'micron'
        kwargs.pixelSize = kwargs.pixelSize * 1e6;
    case 'mm'
        kwargs.pixelSize = kwargs.pixelSize * 1e3;
    case 'm'
        kwargs.pixelSize = kwargs.pixelSize;
end
        
if ~strcmp(kwargs.im, 'none')
    axisObjects = kwargs.im;
else
    axisObjects = ax.Children;
    
data = axisObjects.CData;

right = size(data,2)*kwargs.pixelSize;
top = size(data,1)*kwargs.pixelSize;

pixelPerUnit = 1/kwargs.pixelSize;

% nX = floor(right/500);
% nY = floor(top/500);

% set(axisObjects, 'XData', 1:size(data, 2) * kwargs.pixelSize)
% set(axisObjects, 'YData', 1:size(data, 1) * kwargs.pixelSize)

switch kwargs.location
    case 'bottom left'
        xStart = 0.1*right*pixelPerUnit;
        yStart = 0.15*top*pixelPerUnit;
        textShift = 0.05*top*pixelPerUnit;
    case 'bottom right'
        xStart = 0.85*right;
        yStart = 0.15*top;
        textShift = 0.05*top;
    case 'top left'
        xStart = 0.1*right;
        yStart = 0.9*top;
        textShift = -0.04*top;        
    case 'top right'
        xStart = 0.85*right;
        yStart = 0.9*top;
        textShift = -0.04*top;
end

xEnd = xStart + kwargs.scaleBar * pixelPerUnit;
yEnd = yStart + 0.07*top* pixelPerUnit;
yCenter = mean([yStart, yEnd]);

% left
plot(ax, [xStart, xStart],[yStart, yEnd], '-', 'LineWidth',2, 'color', kwargs.color)
% right
plot(ax, [xEnd, xEnd],[yStart, yEnd], '-', 'LineWidth',2, 'color', kwargs.color)
% line
plot(ax, [xStart, xEnd],[yCenter, yCenter], '-', 'LineWidth',2, 'color', kwargs.color)
text(ax, xStart+(kwargs.scaleBar * pixelPerUnit/2), yStart-textShift, ...
    sprintf('%.0f %s', kwargs.scaleBar, strrep(kwargs.unit, 'micron', '\mum')), ...
    'color', kwargs.color, 'HorizontalAlignment', 'center')

hold(ax, 'off');
end    