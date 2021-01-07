expData = load('D:\data\mike\NRM\run_00000.mat');
%%
binSize = 4;
nRes = 1;
%%
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes, 'winType', 'sine');
%%
binDataNormC = correct_global(binDataNorm, 0.5);
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
[fit, guess, preParams] = fit_resonance(expData, binSize, n, 'type',2);
%%
disp('MOVING window')
[fitNew, guessNew, preParamsNew] = fit_resonance(expData, binSize, n, 'type',2, 'winType','boxcar');