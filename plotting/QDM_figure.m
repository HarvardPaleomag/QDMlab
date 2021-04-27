function map_figure = QDM_figure(data, kwargs)
% Create figure
% Parameter
% ---------
%     fig: figure ['none'];
%         figure object will be used if passed to function. Otherwise one
%         will be created
%     ax: axis ['none']
%         axis object will be used if passed to function, otherwise created
%         at runtime.
%     nROI: ['none']
%         adds ROI to the plot if passed
%     pixelAlerts: array ['none']
%         If pixelAlerts array is passed pixels that, during fitting created
%         an alert (see. fit_resonance) will be replaced by nan
%     filter_hot_pixels: [0]
%         If value (n) >0 pixels will be filtered according to the data
%         with n standard deviations and replaced by nan.
%     title: ['QDM DATA']
%         Title of the axis
%     cbTitle: ['       B_z (T)'];
%         Title of the color bar.
%     axis: ['on']
%         if 'on' x/y labels and box around plot are created,
%         if 'off' x/y labels and box around plot will NOT created
%     return: ['fig']
%         if 'fig' function returns the figure object
%         if 'ax' function returns the axis object. Useful for adding data
%         to a plot

arguments
    data
    kwargs.fig = 'none';
    kwargs.ax = 'none';
    kwargs.nROI = 'none';
    kwargs.pixelAlerts = 'none';
    kwargs.filter_hot_pixels = 0;
    kwargs.title = 'QDM DATA';
    kwargs.cbTitle = 'B_z (T)';
    kwargs.axis = 'on';
    kwargs.return = 'fig';
end

if kwargs.fig == 'none'
    if kwargs.ax == 'none'
        map_figure = figure('Name', 'QDM map', 'units', 'normalized', 'outerposition', [0.2, 0.2, 0.4, 0.6]);
    else
        map_figure = ancestor(kwargs.ax, {'figure'}, 'toplevel');
    end
else
    map_figure = kwargs.fig;
end

if ~all(all(data > 5))
    data = filter_hot_pixels(data);
end

if kwargs.ax == 'none'
    ax = axes('Parent', map_figure);
else
    ax = kwargs.ax;
end

if kwargs.filter_hot_pixels
    data = filter_hot_pixels(data, 'cutOff', kwargs.filter_hot_pixels, 'winSize', nan);
end

if ~strcmp(kwargs.pixelAlerts, 'none')
    data(kwargs.pixelAlerts) = nan;
end

%%
% Create axes
axis off
hold(ax, 'on');

% Create image
% pcolor(data, 'Parent', ax);
imAlpha=ones(size(data));
imAlpha(isnan(data))=0;
imagesc(data,'Parent',ax,'CDataMapping','scaled','AlphaData',imAlpha);

colormap(jet);
shading flat;
set(ax, 'ydir', 'reverse');

% Create title
title(kwargs.title, 'Fontsize', 12);

% box(ax, 'on');
axis(ax, 'tight');

% Set the remaining axes properties
med = abs(median(data, 'all', 'omitnan'));
st = std(data, [], 'all', 'omitnan');
mx = max(abs(data), [], 'all');
mn = min(abs(data), [], 'all');

if ~all(data > 0)
    fprintf('<>   setting Clim: +-%.3f, according to: median (%.3f) + 4*std (%.3f)\n', med+4*st, med, st);
    set(ax, 'CLim', [-1, 1]*(med + 4 * st));
else
    delta = mx - mn;
    set(ax, 'CLim', [med - delta / 10, med + delta / 10]);
    fprintf('<>   setting Clim: (%.3f, %.3f) according to: median (%.3f) +- (max(%.3f)-min(%.3f))/10\n', med-delta/2, med+delta/2, med, mn, mx);
end

axis equal, axis tight, axis xy

if strcmp(kwargs.axis, 'off')
    axis off
end

if iscell(kwargs.nROI)
    add_ROI(kwargs.nROI, 'ax', ax)
end

% Create colorbar
cb = colorbar(ax);
title(cb, kwargs.cbTitle, 'Fontsize', 12);

if strcmp(kwargs.return, 'ax')
    map_figure = ax;
end
