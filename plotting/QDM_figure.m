function  map_figure = QDM_figure(data, kwargs)
%CREATEFIGURE(cdata1)
% Create figure

arguments
    data
    kwargs.figure = 'none';
    kwargs.ax = 'none';
    kwargs.nROI = 'none';
    kwargs.filter_hot_pixels = 0
end

if kwargs.figure == 'none'
    map_figure = figure('Name','QDM map','units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
else
    map_figure = kwargs.figure;
end

if kwargs.ax == 'none'
    ax = axes('Parent',map_figure);
else
    ax = kwargs.ax;
end

if kwargs.filter_hot_pixels
    data = filter_hot_pixels(data, 'cutOff', kwargs.filter_hot_pixels, 'winSize', nan);
end

%%

med = nanmedian(data,'all'); st = nanstd(data,[],'all'); mx = max(abs(data), [], 'all');

% Create axes
axis off
hold(ax,'on');

% prefilter data for hot pixels
data = filter_hot_pixels(data, 'cutOff', 12);

% Create image
imagesc(data,'Parent',ax,'CDataMapping','scaled');

% Create title
title({'QDM DATA'});

box(ax,'on');
axis(ax,'tight');


% Set the remaining axes properties
set(ax,'CLim',[-1 1] * (med + 4*st));
axis equal, axis tight, axis xy

if iscell(kwargs.nROI)
    add_ROI(kwargs.nROI, 'ax', ax)
end

% Create colorbar
colorbar(ax);

