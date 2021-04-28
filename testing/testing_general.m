% testing fit functions
dataFolder = 'D:\data\mike\NRM\'
expData = load('D:\data\mike\NRM\run_00000.mat');
%% fit resonance
% NEW gaussian pre fir
fResNew = fit_resonance(expData, 8, 1)
% OLD local fit
fResLocal = fit_resonance(expData, 8, 1, 'type', 1)
% OLD local fit with gaussian
fResLocalGauss = fit_resonance(expData, 8, 1, 'type', 1, 'gaussianFit', 1)
% OLD global
fResGlobal = fit_resonance(expData, 8, 1, 'type', 0)
%%
figure
ax1 = subplot(2,2,1)
imagesc(fResGlobal.resonance)
axis xy tight equal
colorbar()

ax2 = subplot(2,2,2)
imagesc(fResLocal.resonance)
axis xy tight equal
colorbar()

ax3 = subplot(2,2,3)
imagesc(fResLocalGauss.resonance)
axis xy tight equal
colorbar()

ax4 = subplot(2,2,4)
imagesc(fResNew.resonance)
axis xy tight equal
colorbar()
linkaxes([ax1 ax2 ax3 ax4])

for ax = [ax1 ax2 ax3 ax4]
    set(ax, 'CLim', [2.83, 2.86]);
end

%%
gpuFitNew = GPU_fit(dataFolder, 8)

%% crop_map
crop_map('filePath', '/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/NRM/4x4Binned/Bz_uc0.mat','save',false)