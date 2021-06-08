%expData = load('/Volumes/backups/antares/mike/MIL/MIL3/FOV1/1000G/final_fits_(2x2).mat');
% expData = load('/Volumes/backups/antares/mike/MIL/MIL3/FOV1/NRM/14G/4x4Binned_GF2/final_fits_(4x4).mat');
% expData = load('Z:\antares\mike\MIL\MIL3\FOV1\NRM\4x4Binned\final_fits_(4x4).mat');
% expData = load('/Volumes/backups/antares/tiled_blanks/mariner/FOV_2_2/final_fits_(4x4).mat');
% expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/blanks/2021_05_11/final_fits_(4x4).mat');
expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/speleothems/10-20-10/final_fits_(4x4).mat');
b = expData.B111ferro;
%%  
close all 
clc
data = B111(expData, 'fieldsPlot', 0, 'centerShiftPlot', 0, 'resonanceFieldPlot', 0, 'centerShiftDistribution',1);
d1_ = filter_hot_pixels(data.ferro, 'threshold',false,'cutoff', 15);
d3_ = QuadBGsub(data.centerShift);
% [d1 d2 d3 d4] = B111(expData, 'centerShiftPlot', 1);
% [d1 d2 d3 d4] = B111(expData, 'resonanceFieldPlot', 1);
%%
close all
ferro = QuadBGsub(data.ferro);
ferro = filter_hot_pixels(ferro,'cutoff',2, 'win', nan);
CS = QuadBGsub(data.centerShift);
CS = filter_hot_pixels(CS,'cutoff',2, 'win', nan);
d = [reshape(ferro,[],1), reshape(CS,[],1)];

d = rmmissing(d);
%%
close all
yline(0)
hold on
plot(reshape(d3,[],1), reshape(d1,[],1), '.')
ylim([-0.01,0.01])
xlim([-0.1,0.2])
ylabel('B111ferro')
xlabel('centerShift')
%%
close all
f = figure('units','normalized','outerposition',[0.2 0.4 0.6 0.8],'NumberTitle', 'off', 'Name', 'title');
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,2);
ax3 = subplot(2,2,3);
ax4 = subplot(2,2,4);
GFs = linspace(0,0.01,500);
RMS = rand(size(GFs));

for i = 1: size(GFs,2)
    GF = GFs(i);
    d = d1_ - GF*d3_;
    
    RMS(i) = rms(d,'all');
    
    if i~= 1 && RMS(i) > 1.01*RMS(i-1)
        break
    end
end
[m,idx] = min(RMS);
QDM_figure(d1_, 'ax',ax1, 'title', 'data');
QDM_figure(d1_ - m*d3_, 'ax',ax2, 'title', 'corrected');
QDM_figure(m*d3_, 'ax',ax3, 'title', '\Delta');
plot(ax4, GFs(1:i), RMS(1:i));
set(ax4,'Position', [ax4.Position(1:2), ax1.Position(3:4)]);
xlabel('GF');
ylabel('RMS');
hold(ax4, 'on');
plot(GFs(idx), RMS(idx), 'Xr','LineWidth',2);
title(ax4, sprintf('optimal: %.4f, %s RMS: %.2f %%', GFs(idx), '\Delta', 100*(RMS(idx)/RMS(1)-1)))
linkaxes([ax1 ax2 ax3])

%%
close all
f = figure('units','normalized','outerposition',[0.2 0.4 0.6 0.8],'NumberTitle', 'off', 'Name', 'title');
ax1 = subplot(2,2,1)
QDM_figure(d1_, 'ax',ax1, 'title', 'data', 'clim',[-6e-3,6e-3]);
ax2 = subplot(2,2,2);
[f, ax2, im2] = QDM_figure(d1_ - 0*d3_, 'ax',ax2, 'title', 'corrected', 'clim',[-6e-3,6e-3]);    
ax3 = subplot(2,2,3);
[f, ax3, im3] = QDM_figure(0*d3_, 'ax',ax3, 'title', '\Delta', 'clim',[-1e-3,1e-3]);
ax4 = subplot(2,2,4);

GFs = linspace(0,0.03,10);
RMS = rand(size(GFs));

for i = 1: size(GFs,2)

    GF = GFs(i);
    d = d1_ - GF*d3_;
    
    set(im2, 'CData', d1_ - GF*d3_);
    set(im3, 'CData', GF*d3_);
    
    RMS(i) = rms(d,'all');
    plot(ax4, GFs(1:i), RMS(1:i),'.k','LineWidth',1.5);
    set(ax4,'Position', [ax4.Position(1:2), ax1.Position(3:4)]);
    xlabel(ax4, 'GF');
    ylabel(ax4, 'RMS');
    xlim(ax4, [0,max(GFs)])
    drawnow
end
p = polyfit(GFs,RMS,2);
r = roots(polyder(p));
[m,idx] = min(RMS);
set(im2, 'CData', d1_ - GFs(idx)*d3_);
set(im3, 'CData', GFs(idx)*d3_);
xline(ax4, GFs(idx));
hold(ax4, 'on');
plot(linspace(0,max(GFs),100), polyval(p, linspace(0,max(GFs),100)))
plot(r, rms(d1_ - r*d3_,'all'), 'Xr','LineWidth',2);
plot(GFs(idx), RMS(idx), 'or','LineWidth',2);
title(ax4, sprintf('optimal: %.4f, %s RMS: %.2f %%', GFs(idx), '\Delta', 100*(RMS(idx)/RMS(1)-1)))
linkaxes([ax1 ax2 ax3]);

%%
close all
imagesc(filter_hot_pixels( d1-0*d3, 'threshold', false, 'cutoff', 5))
figure
imagesc(filter_hot_pixels( d1-0.01*d3, 'threshold', false, 'cutoff', 5))

%%
imagesc(filter_hot_pixels( d1/d3, 'threshold', false, 'cutoff', 5))
%%
close all
clc
f = figure('units','normalized','outerposition',[0.15 0.4 0.6 0.6],'NumberTitle', 'off', 'Name', '');

ax1 = subplot(1,1,1);
[f,ax,im] = QDM_figure(expData.ledImg, 'ax', ax1, 'led', true, 'xc', 1:960, 'yc', 1:600, 'title', 'LED');
hold on

QDM_figure(b, 'alpha', 1e-13, 'ax', ax1, 'xc', 1:960, 'yc', 1:600, 'title', 'QDM')