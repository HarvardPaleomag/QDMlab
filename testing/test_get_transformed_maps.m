load('/Users/mike/Dropbox/science/harvard/QDM2_data/electra/mike/MIL3/MIL3folders_mac.mat')
%%
transFormFile = 'none';
fixedIdx = 1;
removeHotPixels = 0;
includeHotPixel = 0;
chi = 0;
winSize = 4;
reverse = 0;
upCont = {0,100};
checkPlot = 1;
%% tranformation / filtering
[transformedData, nFiles] = get_transformed_maps(fov1(1:2), ...
                  'fileName', 'B111dataToPlot', 'transFormFile', transFormFile,...
                  'fixedIdx', fixedIdx, 'removeHotPixels', removeHotPixels,...
                  'includeHotPixel', includeHotPixel,  'chi', chi, ...
                  'upCont', upCont,...
                  'winSize', winSize, 'reverse', reverse, ...
                  'checkPlot', checkPlot);
%%
keys = transformedData.keys();
a = keys{1}; b = keys{2};
%%
close all
f = figure('units','normalized','outerposition',[0.2 0.4 0.5 0.5],'NumberTitle', 'off', 'Name', 'This is the figure title');
ax1 = subplot(2,3,1);
ref = transformedData(a).refData;
ref = filter_hot_pixels(ref);
pcolor(ref);
axis xy; axis equal; axis tight; shading flat;
title('reference Data')

ax2 = subplot(2,3,2);
target = transformedData(a).targetData;
target = filter_hot_pixels(target);
pcolor(target);
axis xy; axis equal; axis tight;shading flat;
title('target Data')

ax3 = subplot(2,3,3);
trans = transformedData(a).transData;
trans = filter_hot_pixels(trans);
pcolor(trans);
axis xy; axis equal; axis tight;shading flat;

ax4 = subplot(2,3,4);
pcolor(transformedData(a).refLed);
axis xy; axis equal; axis tight; shading flat;
colormap(ax4, bone)

title('reference LED')

ax5 = subplot(2,3,5);
pcolor(transformedData(a).transLed);
axis xy; axis equal; axis tight;shading flat;
colormap(ax5, bone)
title('target LED')

ax6 = subplot(2,3,6);
pcolor(transformedData(a).transLed);

axis xy; axis equal; axis tight;shading flat;
title('transformed LED')
colormap(ax6, bone)

linkaxes([ax1 ax2 ax3]);
linkaxes([ax4 ax5 ax6]);