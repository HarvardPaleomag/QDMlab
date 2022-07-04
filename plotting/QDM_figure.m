function [fig, ax, im] = QDM_figure(data, kwargs, filter) 
%[fig, ax, im] = QDM_figure('data', 'fig', 'ax', 'led', 'pixelAlerts', 'title', 'cbTitle', 'unit', 'nROI', 'axis', 'xc', 'yc', 'alpha', 'colormap', 'cLim', 'symmetricCLim', 'method', 'std', 'scaleBar', 'pixelSize', 'nOutlier', 'filterProps', 'preThreshold', 'mustBe')
% Creates a QDM figure
%
% Parameters
% ----------
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
%     nOutlier: int [1]
%         Used to filter the nOutlier lowest and highest values before the
%         colorscale is calculated.

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
    kwargs.colormap = 'parula';
    kwargs.cLim = false;
    kwargs.symmetricCLim = true
    
    kwargs.method = 'std'
    %kwargs.std {mustBeInteger} = 8;
    kwargs.std = 8;
    
    kwargs.scaleBar = false
    kwargs.pixelSize = 4.7e-6
    kwargs.nOutlier = 0
    
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
        kwargs.cbTitle = false
    else
        data = expData.(dataName);
    end
end

%% find figure
if kwargs.fig == 'none'
    if kwargs.ax == 'none'
        fig = figure('Name', 'QDM map', 'units', 'normalized', 'position', [0.1, 0.1, 0.4, 0.6]);
    else
        fig = ancestor(kwargs.ax, {'figure'}, 'toplevel');
    end
else
    fig = kwargs.fig;
end

movegui(fig,'center')

if ~isequal(filter.preThreshold,'none') & ~kwargs.led & ~all(logical(~rem(data(:),1)))
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

if ~isequal(kwargs.unit, 'G') & ~strcmp(kwargs.unit, 'none')
    data = convert_to(data, kwargs.unit);
end

im = imagesc(xc,yc, data,'Parent',ax,'CDataMapping','scaled','AlphaData',imAlpha);

switch kwargs.colormap
    case 'parula'
        colormap(ax, parula(512));
    case 'turbo'
        colormap(ax, turbo(512));
    case 'jet'
        colormap(ax, jet(512));
    case 'rwb'
        colormap(ax, BWR2(512));
end

% Create title
title(ax, kwargs.title, 'Fontsize', 12);

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
if ~strcmp(kwargs.cLim, 'none') & ~isnumeric(kwargs.cLim)
    kwargs.cLim = get_colorscale(data, kwargs.method, 'symmetric', kwargs.symmetricCLim,...
        'std', kwargs.std, 'nOutlier', kwargs.nOutlier, 'mustBe', filter.mustBe);
end

if ~strcmp(kwargs.cLim, 'none')
    set(ax, 'CLim', kwargs.cLim);
end

%% set alpha
if ~isequal(kwargs.alpha, false)
    set(im, 'AlphaData', ~isnan(data));
%     alpha 'none'
%     alpha(im, ~isnan(data)); %kwargs.alpha * abs(data) ./ max(abs(data),[],'all'));
end


%% Create colorbar
if ~isequal(kwargs.cbTitle, false)
    cb = colorbar(ax);
    if isequal(kwargs.led, false)
        unit = strrep(kwargs.unit, 'micro', '\mu');
        unit = strrep(unit, 'mu', '\mu');
        if strcmp(kwargs.unit, 'none')
            title(cb, sprintf('%s', kwargs.cbTitle), 'Fontsize', 12);
        else
            title(cb, sprintf('%s (%s)', kwargs.cbTitle, unit), 'Fontsize', 12);
        end
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
hold(ax, 'off');

end

function bool = isLED(data)
%[bool] = isLED(data)
    data = data(~isnan(data));
    bool = all(logical(~rem(data(:),1)));
end

function h = BWR2(m)
%[h] = BWR2(m)
% red-white-blue color map
%   BWR2(M) returns an M-by-3 matrix containing the colormap.
%   BWR2, by itself, is the same length as the current figure's
%   colormap. If no figure exists,  it creates one as long as MATLAB defautly does.
%-----------------------
% Syntax
%------------------------
%   %Example #1
%   ------------------------------
%   figure;
%   imagesc(peaks(150));
%   colormap(BWR2(20)), colorbar
%   ------------------------------
%   Example #2
%   ------------------------------
%   figure
% bwr1(1)=subplot(1,2,1)
% imagesc(peaks(150), [0 10])
% bwr1(2)=subplot(1,2,2)
% imagesc(peaks(150), [0 10])
% try,colormap(bwr1(1),BWR), title(bwr1(1),'BWR');end;colorbar(bwr1(1));
% colormap(bwr1(2),BWR2);colorbar(bwr1(2)); title(bwr1(2),'BWR2')
%   ------------------------------
%   See also HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG,
%   COLORMAP, RGBPLOT.
%   ------------------------------
%------------------------
%Communication to:  
%------------------------
% Maosheng He
% https://orcid.org/0000-0001-6112-2499 
% Leibniz Institute of Atmospheric Physics  
% Schlossstrasse 6, 18225 Kuehlungsborn
%  hmq512@gmail.com
%  11-Feb., 2020
if nargin < 1, m = size(get(gcf,'colormap'),1); end
x0=[0,1,2,3,4]/4;
ColorCore= [0 0 0.5;0 0.5 1;1 1 1;1 0 0;0.5 0 0];
if mod(m,2)==1,
    h=interp1(x0,ColorCore,[0:1/(m-2):1]);
else
    h1=interp1(x0,ColorCore,[0:1/(m-2):0.5]);
    h2=interp1(x0,ColorCore,0.5+[0:1/(m-2):0.5]);
    h = [h1;h2];
end
end