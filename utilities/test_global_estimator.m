expData = load('/Users/mike/Desktop/NRM/run_00000.mat');
%%
[binDataNorm, freq] = prepare_raw_data(expData, 1, 1);
%%
stdD = nanstd(expData.imgStack1);
pixErr = find(stdD < 1e-5); % find idx of pixels with std close to zero -> all freqs have same value
[~, idxMin] = min(expData.imgStack1);
idxMin(pixErr) = nan;
[~, iMin] = nanmin(idxMin);
[~, iMax] = nanmax(idxMin);
%%
close all
f = figure('units','normalized','outerposition',[0.2 0.4 0.7 0.3],'NumberTitle', 'off', 'Name', 'This is the figure title');

globalMean = squeeze(mean(binDataNorm,[1 2]));
for glob = 0.1:0.1:1
    [x,y] = index2xy(iMin, 1920, 'type', 'binDataNorm');
    dmin = correct_global(binDataNorm(y,x,:), glob, 'mean', globalMean);
    subplot(1,3,1)
    hold on
    plot(squeeze(dmin), 'DisplayName',num2str(glob))
    title('min(min) pixel')
    
    [x,y] = index2xy(iMax, 1920, 'type', 'binDataNorm');
    dmax = correct_global(binDataNorm(y,x,:), glob, 'mean', globalMean);
    subplot(1,3,3)
    hold on
    plot(squeeze(dmax), 'DisplayName',num2str(glob))
    title('max(min) pixel')

    [x,y] = index2xy(randi(size(idxMin,2)), 1920, 'type', 'binDataNorm');
    drand = correct_global(binDataNorm(y,x,:), glob, 'mean', globalMean);
    subplot(1,3,2)
    hold on
    plot(squeeze(drand), 'DisplayName',num2str(glob))
    title('random pixel')
end
%%
close all
f = uifigure('Position',[200 200 800 275],'NumberTitle', 'off', 'Name', 'global correction: 0.00');
% pnl = uipanel(f);

ax1 = axes('Parent',f,'position',[0.1 0.39  0.25 0.54]);
ax2 = axes('Parent',f,'position',[0.4 0.39  0.25 0.54]);
ax3 = axes('Parent',f,'position',[0.7 0.39  0.25 0.54]);

[x1,y1] = index2xy(iMin, 1920, 'type', 'binDataNorm');
plot(ax1, squeeze(binDataNorm(y1,x1,:)), 'k','lineWidth',1)
hold(ax1, 'on')
plot(ax1, globalMean, 'b:')
p1 = plot(ax1, squeeze(binDataNorm(y1,x1,:)),'lineWidth',1)
title(ax1, 'min(min) pixel')

[x2,y2] = index2xy(randi(size(idxMin,2)), 1920, 'type', 'binDataNorm');
plot(ax2, squeeze(binDataNorm(y2,x2,:)), 'k','lineWidth',1)
hold(ax2, 'on')
plot(ax2, globalMean, 'b:')
p2 = plot(ax2, squeeze(binDataNorm(y2,x2,:)),'lineWidth',1);
title(ax2, 'random pixel')


[x3,y3] = index2xy(iMax, 1920, 'type', 'binDataNorm');
plot(ax3, squeeze(binDataNorm(y3,x3,:)), 'k','lineWidth',1)
hold(ax3, 'on')
plot(ax3, globalMean, 'b:')
p3 = plot(ax3, squeeze(binDataNorm(y3,x3,:)),'lineWidth',1);
title(ax3, 'max(min) pixel')

sld = uislider(f, 'Position',[50 50 700 1], 'ValueChangedFcn',@(sld,event) updateGauge(sld, binDataNorm, globalMean, x1,x2,x3,y1,y2,y3, p1,p2,p3, f));
sld.Limits = [0 1];
sld.Value = 0;

% Create ValueChangedFcn callback
function updateGauge(sld, binDataNorm, globalMean, x1,x2,x3,y1,y2,y3, p1,p2,p3, f)
	glob = sld.Value;
    dmin = correct_global(binDataNorm(y1,x1,:), glob, 'mean', globalMean);
    drand = correct_global(binDataNorm(y2,x2,:), glob, 'mean', globalMean);
    dmax = correct_global(binDataNorm(y3,x3,:), glob, 'mean', globalMean);

    set(p1,'ydata',squeeze(dmin));
    set(p2, 'ydata', squeeze(drand));
    set(p3, 'ydata', squeeze(dmax));
    set(f, 'Name', sprintf('global correction: %.3f', glob));
    drawnow;
end

%%
function myslider
x = 1:10;
hplot = plot(x,0*x);
h = uicontrol('style','slider','units','pixel','position',[20 20 300 20]);
addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject, event,x,hplot));
end
function makeplot(hObject,event,x,hplot)
n = get(hObject,'Value');
set(hplot,'ydata',x.^n);
drawnow;

end