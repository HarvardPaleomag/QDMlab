clear all
% expData = load('D:\data\mike\NRM\run_00000.mat');
expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/run_00000.mat');
%%
binSize = 16; nRes = 1;

%% create data like fit_resonance
%% data preparation
% this step could easily be skipped, the only thing one needs to figure out
% is how to get the
[gpuData, freq] = prepare_raw_data(expData, binSize, nRes, 'gpudata',true);
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes, 'gpudata',false);

[sizeY, sizeX, sweepLength] = size(binDataNorm); % binned image x-dimensions
imgPts = sizeX * sizeY; % number of (x,y) pixels

%%

GFs = linspace(0.1, 0.1, 1);
data = [];

binDataNorm = correct_global(binDataNorm, 0.1);
data = correct_global_vec(gpuData, 0.1);
all(transpose(reshape(binDataNorm, [], sweepLength)) == data,'all')

% reshapedData = reshape(data, 51,[],1);
meanData = squeeze(mean(binDataNorm, [1, 2]));
mdata = squeeze(mean(reshapedData, 2));
close all
plot(meanData);
hold on
plot(mdata);
%%
if kwargs.type ~= 2
    initialGuess = global_guess(binDataNorm, freq); % initial guess for GPUfit
end

%% prepare GPUfit data
gpudata = reshape(binDataNorm, [imgPts, sweepLength]); % make it into 2d matrix
gpudata = transpose(gpudata); %transpose to make it 51 x pixels
% gpudata = single(gpudata);

%% testing gc vectorized
specData = squeeze(mean(binDataNorm,[1 2]));%global spectrum 51x1 array
cdnorm = correct_global(binDataNorm, 0.1, 'mean', specData);
cdnorm = transpose(reshape(cdnorm, [], sweepLength));
cgpu = correct_global_vec(gpudata, 0.1, 'mean', specData);
all(cdnorm == cgpu,'all')

close all
plot(squeeze(mean(cdnorm, [1,2])));
hold on
plot(squeeze(mean(cgpu, 2)));