function [fig, ax, im] = QDM_figure(data, kwargs, filter) 
%[fig, ax, im] = QDM_figure(data; 'fig', 'ax', 'led', 'nROI', 'pixelAlerts', 'title', 'cbTitle', 'axis', 'std', 'scaleBar', 'filterStruct', 'preThreshold')
% Creates a QDM figure
%
% Parameter
% ---------
%     data: double [false]
%         Provide either LED or B111/Bz data. If nothing is provided, you
%         get to pick the dataFile to be loaded.
%     fig: figure ['none'];
%         figure object will be used if passed to function. Otherwise one
%         will be created
%     ax: axis ['none']
%         axis object will be used if passed to function, otherwise created
%         at runtime.
%     led: bool [false]
%         to plot LED data
%     nROI: ['none']
%         adds ROI to the plot if passed
%     pixelAlerts: array ['none']
%         If pixelAlerts array is passed pixels that, during fitting created
%         an alert (see. fit_resonance) will be replaced by nan
%     std: int [6]
%         defines how many standard deviations are used to calculate the
%         clims
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
%     filterStruct: struct [empty]
%         a structure to be passed to filter_hot_pixels
%     preThreshold: int [5]
%         Thresholding the map by 'preThreshold' (Gauss), applied before
%         filtering (i.e. see 'filterStruct')

arguments
    data = false;
    kwargs.fig = 'none';
    kwargs.ax = 'none';
    kwargs.led = false;
    kwargs.nROI = 'none';
    kwargs.pixelAlerts = 'none';
    kwargs.title = 'QDM DATA';
    kwargs.cbTitle = 'B_z (T)';
    kwargs.axis = 'on';
    kwargs.std {mustBeInteger} = 10;
    kwargs.scaleBar = false
    
    filter.filterStruct struct = struct();
    filter.preThreshold = 5
end

%% check for data
if isequal(data, false)
    dataFile = automatic_input_ui__('none', 'single', true, 'type', 'file');
    expData = load(dataFile);
    [bool, dataName,ledName] = is_B111(expData);
    
    if kwargs.led
        data = expData.(ledName);
    else
        data = expData.(dataName);
    end
end

%% find figure
if kwargs.fig == 'none'
    if kwargs.ax == 'none'
        fig = figure('Name', 'QDM map', 'units', 'normalized', 'outerposition', [0.2, 0.2, 0.4, 0.6]);
    else
        fig = ancestor(kwargs.ax, {'figure'}, 'toplevel');
    end
else
    fig = kwargs.fig;
end

if filter.preThreshold && ~kwargs.led
    data = filter_hot_pixels(data, 'threshold', filter.preThreshold);
end

if kwargs.ax == 'none'
    ax = gca();
else
    ax = kwargs.ax;
end

if ~all( structfun(@isempty, filter.filterStruct))
    filterProps = namedargs2cell(filter.filterStruct);
    data = filter_hot_pixels(data, filterProps{:});
end

if ~strcmp(kwargs.pixelAlerts, 'none')
    data(kwargs.pixelAlerts) = nan;
end

%%
% Create axes
axis(ax, kwargs.axis);
hold(ax, 'on');

% Create image
% pcolor(data, 'Parent', ax);
imAlpha=ones(size(data));
imAlpha(isnan(data)) = 0;
imagesc(data,'Parent',ax,'CDataMapping','scaled','AlphaData',imAlpha);

xc = 1:size(data, 2);
yc = 1:size(data, 1);

colormap(ax, turbo(512));
shading flat;
set(ax, 'ydir', 'reverse');

% Create title
title(kwargs.title, 'Fontsize', 12);

% box(ax, 'on');
axis(ax, 'tight');
axis equal, axis tight, axis xy

if kwargs.led
    colormap(ax, bone(512));
    return
end

% Set the remaining axes properties
if isnumeric(kwargs.std)
    med = abs(median(data, 'all', 'omitnan'));
    st = std(data, [], 'all', 'omitnan');
    mx = max(data, [], 'all', 'omitnan');
    mn = min(data, [], 'all', 'omitnan');
    
    if (med + kwargs.std * st) > max(abs([mx,mn]))
        msg = sprintf('Clim values exceeds min/max');
        logMsg('debug',msg,1,0);
    elseif ~all(data > 0, 'all')
        msg = sprintf('setting Clim: +-%.3e, according to: median (%.3e) + %i*std (%.3e)', med+kwargs.std*st, med,kwargs.std, st);
        logMsg('debug',msg,1,0);
        set(ax, 'CLim', [-1, 1]*(med + kwargs.std * st));
    else
        msg = sprintf('setting Clim: %.3e:%.3e, according to: median (%.3e) +- %i*std (%.3e)', ...
                     med-kwargs.std*st, med+kwargs.std*st, med,kwargs.std, st);
        logMsg('info',msg,1,0);
        set(ax, 'CLim', [med - kwargs.std * st, med + kwargs.std * st]);
    end
end

if strcmp(kwargs.axis, 'off')
    axis off
end

if iscell(kwargs.nROI)
    add_ROI(kwargs.nROI, 'ax', ax)
end

% Create colorbar
cb = colorbar(ax);
title(cb, kwargs.cbTitle, 'Fontsize', 12);

if ~isequal(kwargs.scaleBar, false)
    msg = sprintf('adding scalebar to the figure. NOTE this is set to a default pixelSize of 4.7e-6 and should be called separately if the size is wrong.');
    logMsg('info',msg,1,0);
    msg = sprintf('e.g. >> [f,a,i] = QDM_figure(Bz); scalebar(''ax'', a, ''scaleBar'', 250, ''location'', ''bottom left'')');
    logMsg('info',msg,1,1);
    scalebar('ax', ax, 'scaleBar', kwargs.scaleBar)
end
