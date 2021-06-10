clear all
close all
clc
% expData = load('D:\data\mike\NRM\run_00000.mat');
expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/run_00000.mat');
%%
binSize = 16; nRes = 1;
[sizeY, sizeX, nFreq] = size(binDataNorm); % binned image x-dimensions
nPixel = sizeX * sizeY; % number of (x,y) pixels

%% create data like fit_resonance
%% make data
[gpuData, freq] = prepare_raw_data(expData, binSize, nRes, 'gpudata',true);
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes, 'gpudata',false);

%% test GPUdata
clc
gpudata = reshape(binDataNorm, [nPixel, nFreq]); % make it into 2d matrix
gpudata = transpose(gpudata); %transpose to make it 51 x pixels
fprintf('<>   GPU data is correct: %i\n', all(gpudata == gpuData,'all'))
%% testing gc vectorized
clc
specData = squeeze(mean(binDataNorm,[1 2]));%global spectrum 51x1 array
cdnorm = correct_global(binDataNorm, 0.1, 'mean', specData);
cdnorm = transpose(reshape(cdnorm, [], nFreq));
cgpu = correct_global_vec(gpuData, 0.1, 'mean', specData);
fprintf('<>   GF correction iscorrect: %i\n', all(cdnorm - cgpu < 0.0001,'all'))
%% GF series test
n = 6;
GFs = linspace(0.1, 0.6, n);
data = [];

for GF = GFs
    data = cat(2, data, correct_global_vec(gpuData, GF));
end

d = reshape(data, nFreq, nPixel, []); % reshape into (f, pixel, gf)
fprintf('<>   GF correction iscorrect: %i\n', all(cdnorm - d(:,:,1) < 0.0001, 'all'))

meanData = squeeze(mean(cdnorm, 2));
mdata = squeeze(mean(d, 2));
close all
plot(meanData,'-');
hold on
plot(mdata(:,1),'--');

%%
preGuess = get_initial_guess(data, freq, 'N14');

%%
initialGuess = [];
for i = 1:6
    initialGuess = cat(2, initialGuess, global_guess(d(:,:,i), freq)); % initial guess for GPUfit
end

%%
fit_resonance(expData, 32, 1, 'type',1)

%% fitting helper functions
function initialGuess = get_initial_guess(gpudata, freq, diamond)
    %[initialGuess] = get_initial_guess(gpudata, freq, diamond)
    
%     if ndims(gpudata) == 3
%         [sizeY, sizeX, nFreq] = size(binDataNorm); % binned image x-dimensions
%         gpudata = reshape(gpudata, nFreq, []);
%     end
    
    initialGuess = zeros(4, size(gpudata, 2), 'single');
    %     n = 1; % cut off outer points
    %     gpudata = gpudata(n:end-n,:);

    % amplitude
    mx = nanmax(gpudata);
    mn = nanmin(gpudata);
    initialGuess(1, :) = -abs(((mx - mn)./mx));

    % center frequency
    [~, idx] = sort(gpudata);
    l = 7; % lowest n values
    mxidx = max(idx(1:l, :));
    mnidx = min(idx(1:l, :));

    if strcmp(diamond, 'N15')
        cIdx = int16((mxidx+mnidx)/2);
    else
        cIdx = int16(mean(cat(1, mxidx, mnidx)));
    end

    center = freq(cIdx);
    initialGuess(2, :) = center;

    % width
    initialGuess(3, :) = 0.0005;

    % offset
    initialGuess(4, :) = mx;
end