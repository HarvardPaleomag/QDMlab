expData = load('D:\data\mike\NRM\run_00000.mat');
%%
binSize = 4; nRes = 1;

%% create data like fit_resonance
%% data preparation
% this step could easily be skipped, the only thing one needs to figure out
% is how to get the
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes);

[sizeY, sizeX, sweepLength] = size(binDataNorm); % binned image x-dimensions
imgPts = sizeX * sizeY; % number of (x,y) pixels

%%
for GF = linspace(0.1, 0.6, 5)
    binDataNorm = correct_global(binDataNorm, GF);
end
%% first determine global guess
meanData = squeeze(mean(binDataNorm, [1, 2]));

if kwargs.type ~= 2
    initialGuess = global_guess(binDataNorm, freq); % initial guess for GPUfit
end

%% prepare GPUfit data
gpudata = reshape(binDataNorm, [imgPts, sweepLength]); % make it into 2d matrix
gpudata = transpose(gpudata); %transpose to make it 51 x pixels
gpudata = single(gpudata);

%%
specData = squeeze(mean(binDataNorm,[1 2]));%global spectrum 51x1 array
cdnorm = correct_global(binDataNorm(1,1,:), 0.1, 'mean', specData);
%%
cgpu = correct_global_vec(gpudata(:,1), 0.1, 'mean', specData);
