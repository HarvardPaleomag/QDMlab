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

%% crop_map
crop_map('filePath', '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/NRM/4x4Binned/Bz_uc0.mat','save',false)

milnrm = '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/4x4Binned/B111dataToPlot.mat';
milnrmBz = '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/4x4Binned/Bz_uc0.mat';
blank = '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/blanks/Mar6_2020_2/4x4Binned/B111dataToPlot.mat';
% load data
MILNRM = load(milnrm);
MILNRMBZ = load(milnrmBz);
BLANK = load(blank);

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
dfit = dipole_fit('filePath', milnrmBz, 'fitOrder', 1, 'cropFactor', 10, 'save',0, 'xy', [334,254]);
pause(2)
close all
dfitc = dipole_fit('filePath', milnrmBz, 'fitOrder', 1, 'cropFactor', 10, 'save',0, 'constrained', true, 'xy', [334,254]);
pause(2)
close all
dfitcMinMax = dipole_fit('filePath', milnrmBz, 'fitOrder', 1, 'cropFactor', 10, 'save',0, 'constrained', true, 'maxheight', 10e-6, 'minheight', 3e-6, 'xy', [334,254]);
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

%%
clc
% dipfitSer = dipole_fit_series(strrep(milnrmBz, '/Bz_uc0.mat', ''));
dipfitSerConst = dipole_fit_series(strrep(milnrmBz, '/Bz_uc0.mat', ''), 'constrained', true, 'minheight', 5e-6, 'maxheight', 5e-6);
%%
f = figure('units','normalized','outerposition',[0.2 0.4 0.3 0.5],'NumberTitle', 'off', 'Name', 'dipole_fit_series constrains');
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