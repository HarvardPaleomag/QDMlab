function ODMR_to_B111_plot(fits, savePath)
%ODMR_to_B111_plot(fits, savePath)

%% PLOTS
rng = .03;

r1 = nanmean(fits.B111para(~fits.pixelAlerts)) - rng;
r2 = nanmean(fits.B111para(~fits.pixelAlerts)) + rng;

[yc, xc] = size(fits.negDiff);
xc = 1:xc; yc = 1:yc;

%f1=figure; imagesc( (negDiff) ); axis equal tight; caxis([-r2 -r1]); colorbar; colormap(gca, turbo(512)); title('Negative current B_{111} (gauss)'); set(gca,'YDir','normal');
[f1, ~, ~] = QDM_figure(fits.negDiff, 'preThreshold', 14, 'title', 'Negative current B_{111} (gauss)', 'cbTitle', 'B_{111} (G)', 'xc', xc,'yc', yc, 'mustBe', 'neg');

%f2=figure; imagesc( (posDiff) ); axis equal tight; caxis([r1 r2]); colorbar; colormap(gca, turbo(512)); title('Positive current B_{111} (gauss)'); set(gca,'YDir','normal');
[f2, ~, ~] = QDM_figure(fits.posDiff, 'preThreshold', 14, 'title', 'Positive current B_{111} (gauss)', 'cbTitle', 'B_{111} (G)', 'xc', xc,'yc', yc, 'mustBe', 'pos');

%f3=figure; imagesc( B111ferro ); axis equal tight; caxis(-.1 + [-rng rng]); colorbar; colormap(gca, turbo(512)); title('Positive + negative ferro B_{111} (gauss)'); set(gca,'YDir','normal');
[f3, ~, ~] = QDM_figure(fits.B111ferro, 'preThreshold', 14, 'title', 'Positive + negative ferro B_{111} (gauss)', 'cbTitle', 'B_{111} (G)');

%f4=figure; imagesc(ledImg); axis equal tight; colorbar; colormap(gca, gray(512)); caxis auto; title('LED image'); set(gca,'YDir','normal');
[f4, ~, ~] = QDM_figure(fits.ledImg, 'preThreshold', 14, 'led', true, 'title', 'LED image', 'xc', xc,'yc', yc);

%f5=figure; imagesc( B111para ); axis equal tight; caxis([r1 r2]); colorbar; colormap(gca, turbo(512)); title('Positive + negative para B_{111} (gauss)'); set(gca,'YDir','normal');
[f5, ~, ~] = QDM_figure(fits.B111para, 'preThreshold', 14, 'title', 'Positive - negative ferro B_{111} (gauss)', 'cbTitle', 'B_{111} (G)', 'xc', xc,'yc', yc, 'mustBe', 'pos');

%% combined figure
f6 = figure('Name', 'data', 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
set(gca, 'YDir', 'normal');
s1 = subplot(2, 2, 1); %imagesc( (negDiff) ,'hittest', 'off'); axis equal tight; caxis auto; colorbar; colormap(s1,turbo(512)); title('Negative current B_{111} (gauss)'); set(gca,'YDir','normal');
QDM_figure(fits.negDiff, 'ax', s1, 'preThreshold', 14, 'title', 'Negative current B_{111} (gauss)', 'cbTitle', 'B_{111} (G)', 'xc', xc,'yc', yc, 'mustBe', 'neg');
s2 = subplot(2, 2, 2); %imagesc( (posDiff) ,'hittest', 'off'); axis equal tight; caxis auto; colorbar; colormap(s2,turbo(512)); title('Positive current B_{111} (gauss)'); set(gca,'YDir','normal');
QDM_figure(fits.posDiff, 'ax', s2, 'preThreshold', 14, 'title', 'Positive current B_{111} (gauss)', 'cbTitle', 'B_{111} (G)', 'xc', xc,'yc', yc, 'mustBe', 'pos');
s3 = subplot(2, 2, 3); %imagesc( B111ferro ,'hittest', 'off'); axis equal tight;  caxis(mean2(B111ferro) + [-rng rng]); colorbar; colormap(s3,turbo(512)); title('Positive + negative ferro B_{111} (gauss)'); set(gca,'YDir','normal');
QDM_figure(fits.B111ferro, 'ax', s3, 'preThreshold', 14, 'title', 'Positive + negative ferro B_{111} (gauss)', 'cbTitle', 'B_{111} (G)', 'xc', xc,'yc', yc);
s4 = subplot(2, 2, 4); %imagesc( (ledImg) ,'hittest', 'off'); axis equal tight; colorbar; colormap(s4,gray(512)); caxis auto; title('LED image'); set(gca,'YDir','normal');
QDM_figure(fits.ledImg, 'ax', s4, 'preThreshold', 14, 'led', true, 'title', 'LED image', 'xc', xc,'yc', yc);
sgtitle(', B111 points up and out of page');
ax = [s1, s2, s3, s4];
linkaxes(ax);

map = [1, 1, 1; 1, 0, 0];
f7 = figure;
imagesc(xc,yc, fits.pixelAlerts);
axis equal tight;
colormap(gca, map);
title('pixel alerts');
set(gca, 'YDir', 'normal');

if ~ isequal(savePath, false)
    saveas(f1, fullfile(savePath, 'negCurrent.png'), 'png');
    saveas(f2, fullfile(savePath, 'posCurrent.png'), 'png');
    saveas(f3, fullfile(savePath, 'ferromagImg.png'), 'png');
    saveas(f4, fullfile(savePath, 'ledImg.png'), 'png');
    saveas(f5, fullfile(savePath, 'paramagImg.png'), 'png');
    saveas(f6, fullfile(savePath, 'allPlots.png'), 'png');
    saveas(f7, fullfile(savePath, 'pixelAlerts.png'), 'png');
end
end