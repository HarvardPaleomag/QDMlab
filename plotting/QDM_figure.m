function  map_figure = QDM_figure(data, kwargs)
%CREATEFIGURE(cdata1)
% Create figure

arguments
    data
    kwargs.figure = 'none';
    kwargs.ax = 'none';
end

if kwargs.figure == 'none'
    map_figure = figure('Name','QDM map','units','normalized','outerposition',[0 0 1 1]);
else
    map_figure = kwargs.figure;
end


if kwargs.ax == 'none'
    ax = axes('Parent',kwargs.figure);
else
    ax = kwargs.ax;
end

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

mx = max(abs(data), [], 'all');

% Set the remaining axes properties
set(ax,'CLim',[-1 1] * mx/10);
axis equal, axis tight, axis xy

% Create colorbar
colorbar(ax);

