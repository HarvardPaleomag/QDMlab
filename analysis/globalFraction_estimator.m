function globalFraction = globalFraction_estimator(expData, kwargs)
%[globalFraction] = globalFraction_estimator('expData', 'header', 'binSize', 'nRun', 'nRes')
%
% Parameters
% ----------
%     expData: path, struct
%         The path to 'run0000n.mat' or the structure retured from load(run0000n.mat)
%     binSize: int
%         binning size (can be 1)
%     nRes: int
%         number of resonance. Low frequencies = 1, High frequencies = 2
    
%% load/prepare data

arguments
    expData = 'none'
    kwargs.header = 'none'
    kwargs.binSize = 'none'
    kwargs.nRun = 'none'
    kwargs.nRes = 'none'
end


%% load data
expData = automatic_input_ui__(expData, 'type', 'file', 'single', true, 'title','select run_000*.mat file');

% if no expData structure was passed, the file needs to be loaded
if ~isstruct(expData)
    msg = sprintf('loading: %s', expData);
    logMsg('debug',msg,1,0);
    expData = load(expData);
end

defaults = struct('binSize', 4, 'nRun',1, 'nRes', 1);
kwargs = ask_arguments(kwargs, defaults);

[binDataNorm, freq] = prepare_raw_data(expData, kwargs.binSize, kwargs.nRes, kwargs.header);
nCol = size(binDataNorm, 2);

%% remove non fitting pixels
% transform data back into picel by pixel col by col data
nFreq = size(binDataNorm, 3);
imgStack = QDMreshape_reverse(binDataNorm, nFreq);
stdD = std(imgStack, 'omitnan');
pixErr = stdD < 1e-3; % find idx of pixels with std close to zero -> all freqs have same value

%% determine indices of binned image
[~, idxMin] = min(imgStack, [], 'omitnan');
idxMin(pixErr) = nan;

[~, iMin] = min(idxMin,[], 'omitnan'); % left most minimum
[~, iMax] = max(idxMin,[], 'omitnan'); % right most minimum

%%
globalMean = squeeze(mean(binDataNorm, [1,2]));

%% figure
% initialize UI figure
f = figure('Units', 'normalized', ...
           'Position',[0.1 0.2 0.8 0.4],'NumberTitle', 'off', 'Name', 'global correction: 0.00');

% define axes objects
ax1 = axes('Parent',f,'position',[0.1 0.3  0.25 0.54]);
ax2 = axes('Parent',f,'position',[0.4 0.3  0.25 0.54]);
ax3 = axes('Parent',f,'position',[0.7 0.3  0.25 0.54]);

%% plot pixels
%% left most pixel
[y1,x1] = index2xy(iMin, size(binDataNorm), 'type', 'raw');
plot(ax1, freq, squeeze(binDataNorm(y1,x1,:)), 'k','lineWidth',1);
hold(ax1, 'on');
plot(ax1, freq, globalMean, 'b:','lineWidth',1)
p1 = plot(ax1, freq, squeeze(binDataNorm(y1,x1,:)),'lineWidth',1);
title(ax1, 'min(min) pixel')

%% center pixel random
randInt = randi(size(idxMin,2));
[y2,x2] = index2xy(randInt, size(binDataNorm), 'type', 'raw');
p2_data = plot(ax2, freq, squeeze(binDataNorm(y2,x2,:)), 'k','lineWidth',1);
hold(ax2, 'on');
globalPlot = plot(ax2,freq, globalMean, 'b:','lineWidth',1);
p2 = plot(ax2, freq, squeeze(binDataNorm(y2,x2,:)),'lineWidth',1);
title(ax2, sprintf('random pixel [%i] (%i,%i)',randInt,x2,y2));
xlim(ax2, [min(freq), max(freq)]);
legend([p2_data, globalPlot, p2], 'data', 'global', 'corrected','Location','southeast')
legend('boxoff')

%% right most pixel
[y3,x3] = index2xy(iMax, size(binDataNorm), 'type', 'raw');
plot(ax3, freq, squeeze(binDataNorm(y3,x3,:)), 'k','lineWidth',1);
hold(ax3, 'on');
plot(ax3, freq, globalMean, 'b:','lineWidth',1);
p3 = plot(ax3, freq, squeeze(binDataNorm(y3,x3,:)),'lineWidth',1);
title(ax3, 'max(min) pixel');

%% plot setup
for ax = [ax1 ax2 ax3]
    xlim(ax, [min(freq), max(freq)]);
    xlabel(ax, 'f (GHz)')
end

%% sliders
% select pixel option
nPixels = size(binDataNorm,1) * size(binDataNorm,2);
idx = uicontrol(f, ...
    'Units', 'normalized', 'Position',[0.15 0.91 0.8 0.05],...
    'Style','slider','Min',1,'Max',nPixels,'SliderStep',[1/nPixels nCol/nPixels]);
idx.Value = randInt;
label1 = uicontrol('style','text', ...
                   'HorizontalAlignment', 'right', ...,
                   'FontSize',12,...
                   'Units', 'normalized', 'Position',[0.04 0.855 0.1 0.1]);
label1.String = 'index';

% select globalFraction
sld = uicontrol(f, ...
            'Units', 'normalized', 'Position',[0.15 0.11 0.8 0.05],...
            'Style','slider','Min',0,'Max',1,'SliderStep',[0.01 0.10]);
sld.Value = 0;
label1 = uicontrol('style','text', ...
                   'HorizontalAlignment', 'right', ...
                   'FontSize',12,...
                   'Units', 'normalized', 'Position',[0.04 0.055 0.1 0.1]);
label1.String = 'globalFraction [0.00]';

% listeners for callback
addlistener(sld, 'Value', 'PostSet',@(src,event) updatePlot(src, event, binDataNorm, globalMean, x1,x2,x3,y1,y2,y3, p1,p2,p3,f));
addlistener(idx, 'Value', 'PostSet',@(src,event) updateRand(src, event, binDataNorm, globalMean, p2, p2_data, sld));

% functions
function updateRand(src, event, binDataNorm, globalMean, p2,p2_data, sld)
%updateRand(src, event, binDataNorm, globalMean, p2, p2_data, sld)

    idx.Value = round(idx.Value);
    [y2,x2] = index2xy(idx.Value, size(binDataNorm), 'type', 'raw');

    drand = correct_global(binDataNorm(y2,x2,:), sld.Value, 'mean', globalMean);
    set(p2_data, 'ydata', squeeze(binDataNorm(y2,x2,:)))
    set(p2, 'ydata', squeeze(drand));
    title(ax2, sprintf('random pixel [%i] (%i,%i)',randInt,x2,y2));

end

% Create ValueChangedFcn callback
function updatePlot(src, event, binDataNorm, globalMean, x1,x2,x3, y1,y2,y3, p1,p2,p3, f)
%updatePlot(src, event, binDataNorm, globalMean, x1, x2, x3, y1, y2, y3, p1, p2, p3, f)
    glob = sld.Value;
    dmin = correct_global(binDataNorm(y1,x1,:), glob, 'mean', globalMean);
    drand = correct_global(binDataNorm(y2,x2,:), glob, 'mean', globalMean);
    dmax = correct_global(binDataNorm(y3,x3,:), glob, 'mean', globalMean);

    set(p1,'ydata',squeeze(dmin));
    set(p2, 'ydata', squeeze(drand));
    set(p3, 'ydata', squeeze(dmax));
    set(f, 'Name', sprintf('global correction: %.2f', glob));
    set(label1, 'String', sprintf('globalFraction [%.2f]', glob));

    drawnow;
end
end