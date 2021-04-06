function  map_figure = QDM_figure(data, kwargs)
%CREATEFIGURE(cdata1)
% Create figure

arguments
    data
    kwargs.fig = 'none';
    kwargs.ax = 'none';
    kwargs.nROI = 'none';
    kwargs.fitSuccess = 'none';
    kwargs.filter_hot_pixels = 0;
    kwargs.title = 'QDM DATA';
end

if kwargs.fig == 'none'
    if kwargs.ax == 'none'
        map_figure = figure('Name','QDM map','units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
    else
        map_figure = ancestor(kwargs.ax,{'figure'},'toplevel');
    end
else
    map_figure = kwargs.fig;
end

if ~all(all(data > 5))
    data = filter_hot_pixels(data);
end

if kwargs.ax == 'none'
    ax = axes('Parent',map_figure);
else
    ax = kwargs.ax;
end

if kwargs.filter_hot_pixels
    data = filter_hot_pixels(data, 'cutOff', kwargs.filter_hot_pixels, 'winSize', nan);
end

if ~strcmp(kwargs.fitSuccess, 'none')
    data(~kwargs.fitSuccess) = nan;
end
%%
% Create axes
axis off
hold(ax,'on');

% Create image
pcolor(data,'Parent',ax);
shading flat;
set(ax, 'ydir', 'reverse');
% imagesc(data,'Parent',ax,'CDataMapping','scaled');

% Create title
title(kwargs.title);

box(ax,'on');
axis(ax,'tight');

% Set the remaining axes properties
med = abs(median(data,'all','omitnan')); st = std(data,[],'all','omitnan'); 
mx = max(abs(data), [], 'all'); mn = min(abs(data), [], 'all');

if ~all(data>0)
    fprintf('<>    setting Clim: +-%.3f, according to: median (%.3f) + 4*std (%.3f)\n', med + 4*st, med, st);
    set(ax,'CLim',[-1 1] * (med + 4*st));
else
    delta = mx-mn;
    set(ax,'CLim',[med-delta/10 med+delta/10]);
    fprintf('<>    setting Clim: (%.3f, %.3f) according to: median (%.3f) +- (max(%.3f)-min(%.3f))/10\n',med-delta/2, med+delta/2, med, mn, mx);
end

axis equal, axis tight, axis xy

if iscell(kwargs.nROI)
    add_ROI(kwargs.nROI, 'ax', ax)
end

% Create colorbar
colorbar(ax);

