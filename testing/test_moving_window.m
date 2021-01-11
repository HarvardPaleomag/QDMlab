% expData = load('D:\data\mike\NRM\run_00000.mat');
expData = load('/Users/mike/Desktop/NRM/run_00000.mat');
%%
binSize = 4;
nRes = 1;
%%
data = zeros(expData.imgNumRows, expData.imgNumCols, expData.numFreqs);

% reshape and transpose each image
for y = 1:expData.numFreqs
    data(:,:,y) = transpose(reshape(expData.imgStack1(y, :), [expData.imgNumCols, expData.imgNumRows]));
end
%% only moving Win
bin = moving_bin(data, 4);

%%
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes, 'winType', 'boxcar');
%% OLD style
[binDataNorm, freq] = prepare_raw_data(expData, 2, nRes, 'winType', 'none');
%%
binDataNormC = correct_global(binDataNorm, 0.1);
%%
close all
r = randi(1200,1); c = randi(1900,1);
hold on
plot(squeeze(binDataNorm(r,c,:)), '.k')
plot(squeeze(binDataNormC(r,c,:)), 'r--')
title(sprintf('X:%i Y:%i',c,r))
%%
m1 = nanmean(binDataNorm,1);
m2 = nanmean(squeeze(m1),1);
meanData = squeeze(m2);
%%
dataStack = expData.(sprintf('imgStack%i',1));
sweeplength = size(dataStack,1);
sizeX = size(binDataNorm,2); % binned image x-dimensions
sizeY = size(binDataNorm,1); % binned image y-dimensions
imgpts = sizeX*sizeY; % number of (x,y) pixels
gpudata = reshape(binDataNorm,[imgpts,sweeplength]); % make it into 2d matrix
%%
initialGuess = global_guess(binDataNorm, freq); % initial guess for GPUfit
%%
disp('LOCAL GAUSSIAN')
[fit, guess, preParams] = fit_resonance(expData, binSize, nRes, 'type',2);
%%
disp('MOVING window')
[fitNew, guessNew, preParamsNew] = fit_resonance(expData, binSize, nRes, 'type',2, 'winType','boxcar');
%%
figure
ax1 = subplot(1,2,1)
imagesc()
axis xy; axis equal; axis tight
ax1 = subplot(1,2,2)
imagesc()
axis xy; axis equal; axis tight
linkaxes([ax1 ax2])
