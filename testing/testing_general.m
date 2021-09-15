% testing fit functions
dataFolder = 'D:\data\mike\NRM\'
expData = load('D:\data\mike\NRM\run_00000.mat');
%% fit resonance
% NEW gaussian pre fir
fResNew = fit_resonance(expData, 8, 1);
% OLD local fit
fResLocal = fit_resonance(expData, 8, 1, 'type', 1);
% OLD local fit with gaussian
fResLocalGauss = fit_resonance(expData, 8, 1, 'type', 1, 'gaussianFit', 1);
% OLD global
fResGlobal = fit_resonance(expData, 8, 1, 'type', 0);
%%
figure
ax1 = subplot(2,2,1);
imagesc(fResGlobal.resonance);
axis xy tight equal;
colorbar()

ax2 = subplot(2,2,2);
imagesc(fResLocal.resonance);
axis xy tight equal;
colorbar();

ax3 = subplot(2,2,3);
imagesc(fResLocalGauss.resonance);
axis xy tight equal;
colorbar();

ax4 = subplot(2,2,4);
imagesc(fResNew.resonance);
axis xy tight equal;
colorbar();
linkaxes([ax1 ax2 ax3 ax4]);

for ax = [ax1 ax2 ax3 ax4]
    set(ax, 'CLim', [2.83, 2.86]);
end

%%
gpuFitNew = GPU_fit(dataFolder, 8);
%% Load data
milnrm = '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/4x4Binned/B111dataToPlot.mat';
milnrm_ = '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/4x4Binned';
milnrmBz = '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/4x4Binned/Bz_uc0.mat';
blank = '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/blanks/Mar6_2020_2/4x4Binned/B111dataToPlot.mat';
raw_data = '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/run_00000.mat';
% load data
MILNRM = load(milnrm);
MILNRMBZ = load(milnrmBz);
BLANK = load(blank);
expData = load(raw_data);
%%
globalFraction_estimator('expData', expData)

%% return_bin_data
close all
clc
cIdx = 359; rIdx = 129;
[binDataNorm, freq] = prepare_raw_data(expData, 4, 1);
binned = squeeze(binDataNorm(rIdx,cIdx,:));
idx = return_bin_data(cIdx, rIdx);
[x,y] = index2xy(idx, 300);
unbinned = expData.imgStack1(:, idx);
unbinned = unbinned ./ max(unbinned);
plot(unbinned)
hold on
plot(binned, 'k-', 'lineWidth', 2)
hold off
QDM_figure(MILNRM.B111ferro)
hold on
plot(x,y,'Xr')
%% crop_map
cropped_map = crop_map('filePath', '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/NRM/4x4Binned/Bz_uc0.mat','save',false);
%% subtract blank
close all
d = subtract_blank('nFiles', '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/4x4Binned/B111dataToPlot.mat', ...
    'blankFile', '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/blanks/Mar6_2020_2/4x4Binned/B111dataToPlot.mat');
dkeys = d.keys();
k = dkeys{1};
BLANKSUBTRACTED = d(k);
f = figure('units', 'normalized', 'outerposition', [0.2, 0.4, 0.5, 0.5], 'NumberTitle', 'off', 'Name', 'This is the figure title');
ax1 = subplot(1, 3, 1);
QDM_figure(MILNRM.B111ferro, 'ax', ax1)
axis xy; axis equal; axis tight;
ax2 = subplot(1, 3, 2);
QDM_figure(BLANK.B111ferro, 'ax', ax2)
ax3 = subplot(1, 3, 3);
QDM_figure(BLANKSUBTRACTED.B111ferro, 'ax', ax3)
linkaxes([ax1, ax2, ax3]);

%% subtract_constant
close all
d = subtract_constant('filePath', milnrm, 'editFigure', false);
f = figure('units', 'normalized', 'outerposition', [0.2, 0.4, 0.5, 0.5], 'NumberTitle', 'off', 'Name', 'This is the figure title');
ax1 = subplot(1, 2, 1);
QDM_figure(MILNRM.B111ferro, 'ax', ax1);
axis xy; axis equal; axis tight;
ax2 = subplot(1, 2, 2);
QDM_figure(d.B111ferro, 'ax', ax2);
linkaxes([ax1, ax2]);

%% subtract_source
clc
close all
dipSub = subtract_source('filePath', milnrmBz);
f = figure('units','normalized','outerposition',[0.2 0.4 0.5 0.5],'NumberTitle', 'off', 'Name', 'This is the figure title');
ax1 = subplot(1,2,1);
QDM_figure(dipSub.Bz_original, 'ax', ax1);
axis xy; axis equal; axis tight;

ax2 = subplot(1,2,2);
QDM_figure(dipSub.Bz, 'ax', ax2);
linkaxes([ax1 ax2]);

%% dipole fit
close all; clc;
dfit = fit_source('filePath', milnrmBz, 'fitOrder', 1, 'cropFactor', 10, 'save',0, 'xy', [334,254]);
pause(2)
close all
dfitc = fit_source('filePath', milnrmBz, 'fitOrder', 1, 'cropFactor', 10, 'save',0, 'constrained', true, 'xy', [334,254]);
pause(2)
close all
dfitcMinMax = fit_source('filePath', milnrmBz, 'fitOrder', 1, 'cropFactor', 10, 'save',0, 'constrained', true, 'maxheight', 10e-6, 'minheight', 3e-6, 'xy', [334,254]);
pause(2)
close all
f = figure('units','normalized','outerposition',[0.2 0.4 0.4 0.5],'NumberTitle', 'off', 'Name', 'Dipole fitting comparison');
ax1 = subplot(2,2,1);
QDM_figure(dfit.data, 'ax', ax1);
axis xy; axis equal; axis tight;

ax2 = subplot(2,2,2);
QDM_figure(dfit.model, 'ax', ax2, 'title', 'unconstrained');
ax3 = subplot(2,2,3);
QDM_figure(dfitc.model, 'ax', ax3, 'title', 'constrained');
ax4 = subplot(2,2,4);
QDM_figure(dfitcMinMax.model, 'ax', ax4, 'title', 'constrained h = [3e-6,10e-6]');
linkaxes([ax1 ax2 ax3 ax4]);
%% fit_source_series
clc
close all
dipfitSer = fit_source_series(strrep(milnrmBz, '/Bz_uc0.mat', ''));
dipfitSerConst = fit_source_series(strrep(milnrmBz, '/Bz_uc0.mat', ''), 'constrained', true, 'minheight', 5e-6, 'maxheight', 5e-6, 'checkPlot', true);
f = figure('units','normalized','outerposition',[0.2 0.4 0.3 0.5],'NumberTitle', 'off', 'Name', 'fit_source_series constrains');
ax1 = subplot(2,2,1);
QDM_figure(dipfitSer.data{1}, 'ax', ax1);
axis xy; axis equal; axis tight;
ax2 = subplot(2,2,2);
QDM_figure(dipfitSer.model{1}, 'ax', ax2);
ax3 = subplot(2,2,3);
QDM_figure(dipfitSerConst.data{1}, 'ax', ax3);
axis xy; axis equal; axis tight;
ax4 = subplot(2,2,4);
QDM_figure(dipfitSerConst.model{1}, 'ax', ax4);

linkaxes([ax1 ax2])
linkaxes([ax3 ax4]);

%%
clear all
folders = { ...
    '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/NRM/4x4Binned', ...
    '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/150G/4x4Binned', ...
    '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/280G_2/4x4Binned', ...
    '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/400G/4x4Binned', ...
    '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/500G/4x4Binned', ...
    '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/600G/4x4Binned', ...
    '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/1000G/4x4Binned', ...
    };
dseries = fit_source_series(folders);
%%
clc
close all
for n = 1:size(dseries.nFiles,1)
    f = figure('units','normalized','outerposition',[0.05 0.4 0.9 0.4],'NumberTitle', 'off', 'Name', 'Dipole fit Series');
    axes = [];
    for i = 1:size(dseries.nFiles,2)
        ax = subplot(2,7,i);
        QDM_figure(dseries.data{n,i}, 'ax', ax);
        ax = subplot(2,7,7+i);
        QDM_figure(dseries.model{n,i}, 'ax', ax);
        axes = [axes ax];
    end
    linkaxes(axes);
end

%% demag_behavior
% only argument
res = demag_behavior({milnrm_, milnrm_});
demag_behavior_plot(res);
%% filter
res = demag_behavior({milnrm_, milnrm_}, 'filterProps', struct('threshold', 1, 'cutOff', 5));
demag_behavior_plot(res);
%% nROI
nROI = pick_box(MILNRM.B111ferro);
res = demag_behavior({milnrm_, milnrm_}, 'filterProps', struct('threshold', 1, 'cutOff', 5), 'nROI', nROI);
demag_behavior_plot(res);
%% error
res = demag_behavior({milnrm_, milnrm_}, 'filterProps', struct('threshold', 1, 'cutOff', 5), 'nROI', nROI, 'bootStrapN', 500, 'pixelShift', 1);
demag_behavior_plot(res);
