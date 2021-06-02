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
    kwargs.pixelAlerts = 'none';

    kwargs.title = 'QDM DATA';
    kwargs.cbTitle = 'B_z';
    kwargs.unit = 'G'; % assuming the input data is in G
    
    kwargs.nROI = 'none';

    kwargs.axis = 'on';
    kwargs.xc = false;
    kwargs.yc = false;
    kwargs.alpha = false;
    
    kwargs.std {mustBeInteger} = 10;
    
    kwargs.scaleBar = false
    kwargs.pixelSize = 4.7e-6
    
    filter.filterProps struct = struct();
    filter.preThreshold = 5
    filter.mustBe = false;
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

if filter.preThreshold & ~kwargs.led
    data = filter_hot_pixels(data, 'threshold', filter.preThreshold);
end

if kwargs.ax == 'none'
    ax = gca();
else
    ax = kwargs.ax;
end

if ~all( structfun(@isempty, filter.filterProps))
    filterProps = namedargs2cell(filter.filterProps);
    data = filter_hot_pixels(data, filterProps{:});
end

if ~strcmp(kwargs.pixelAlerts, 'none')
    data(kwargs.pixelAlerts) = nan;
end

%% remove +/- alues if they should not be there
if ~isequal(filter.mustBe, false)
    switch filter.mustBe
        case 'neg'
            data(data>0) = nan;
        case 'pos'
            data(data<0) = nan;
    end
end

%% Create image
imAlpha=ones(size(data));
imAlpha(isnan(data)) = 0;

if isequal(kwargs.xc, false)
    xc = 1:size(data, 2);
else
    xc= kwargs.xc;
end

if isequal(kwargs.yc, false)
    yc = 1:size(data, 1);
else
    yc= kwargs.yc;
end

data = convert_to(data, kwargs.unit);
im = imagesc(xc,yc, data,'Parent',ax,'CDataMapping','scaled','AlphaData',imAlpha);

colormap(ax, parula(512));

% Create title
title(kwargs.title, 'Fontsize', 12);

%% 
% Create axes
axis(ax, kwargs.axis);
box(ax, 'on');
axis(ax, 'tight');
axis(ax, 'equal')
axis(ax, 'xy')
hold(ax, 'on');

%% led
if kwargs.led
    colormap(ax, gray(512));
end

if iscell(kwargs.nROI)
    add_ROI(kwargs.nROI, 'ax', ax)
    if kwargs.led
        return
    end
end

% Set the remaining axes properties
if isnumeric(kwargs.std)
    med = median(abs(data), 'all', 'omitnan');
    st = std(data, [], 'all', 'omitnan');
    mx = max(data, [], 'all', 'omitnan');
    mn = min(data, [], 'all', 'omitnan');
    
    try
        if (med + kwargs.std * st) > max(abs([mx,mn]))
            msg = sprintf('Clim values exceeds min/max');
            logMsg('debug',msg,1,0);
%         elseif all(data(~isnan(data)) > 0, 'all') | all(data(~isnan(data)) < 0, 'all')
%             set(ax, 'CLim', [mn, mx]);
        elseif ~all(data(~isnan(data)) > 0, 'all')
            msg = sprintf('setting Clim: +-%.3e, according to: median (%.3e) + %i*std (%.3e)', med+kwargs.std*st, med,kwargs.std, st);
            logMsg('debug',msg,1,0);
            set(ax, 'CLim', [-1, 1]*(med + kwargs.std * st));
        else
            msg = sprintf('setting Clim: %.3e:%.3e, according to: median (%.3e) +- %i*std (%.3e)', ...
                         med-kwargs.std*st, med+kwargs.std*st, med,kwargs.std, st);
            logMsg('info',msg,1,0);
            set(ax, 'CLim', [med - kwargs.std * st, med + kwargs.std * st]);
        end
    catch
        return
    end
end

%% set alpha
if ~isequal(kwargs.alpha, false)
    im.AlphaData = abs(data) ./ max(abs(data),[],'all');
    alpha 'none'
%     alpha(im, kwargs.alpha * abs(data) ./ max(abs(data),[],'all'));
end

%% Create colorbar
if ~ isequal(kwargs.cbTitle, false)
    cb = colorbar(ax);
    if isequal(kwargs.led, false)
        title(cb, sprintf('%s (%s)', kwargs.cbTitle, strrep(kwargs.unit, 'micro', '\mu')), 'Fontsize', 12);
    else
        title(cb, '', 'Fontsize', 12);
    end
end

%% scalebar
if ~isequal(kwargs.scaleBar, false)
    msg = sprintf('adding scalebar to the figure. NOTE this is set to a default pixelSize of 4.7e-6 and should be called separately if the size is wrong.');
    logMsg('info',msg,1,0);
    msg = sprintf('e.g. >> [f,a,i] = QDM_figure(Bz); scalebar(''ax'', a, ''scaleBar'', 250, ''location'', ''bottom left'')');
    logMsg('info',msg,1,1);
    scalebar('ax', ax, 'scaleBar', kwargs.scaleBar, 'pixelSize', kwargs.pixelSize)
end

end

function data = convert_to(data, unit)
    switch unit
        case 'T'
            conv = 0.0001;
        case 'microT'
            conv = 0.1;
        case 'nT'
            conv = 100;
        case 'G'
            conv = 1;
    end
    msg = sprintf('converting 1 G -> %i%s: ', conv, unit);
    logMsg('debug',msg,1,0);
    data = data * conv;
end
