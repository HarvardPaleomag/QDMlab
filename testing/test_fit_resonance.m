clear all
close all
clc
% expData = load('D:\data\mike\VISC\ALH84001\S10\FOV1_VISC_20-20\run_00000.mat');

expData = load('D:\data\mike\NRM\run_00000.mat');
% expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/run_00000.mat');
%%
binSize = 64; nRes = 1;

%% create data like fit_resonance
%% make data
[gpuData, freq] = prepare_raw_data(expData, binSize, nRes, 'gpudata',true);
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes, 'gpudata',false);
[sizeY, sizeX, nFreq] = size(binDataNorm); % binned image x-dimensions
nPixel = sizeX * sizeY; % number of (x,y) pixels
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
%%
close all
clc
c = 4; r = 7;
idx = xy2index(r,c,size(binDataNorm),'type','gpu')
a = squeeze(binDataNorm(r,c,:));
b = gpuData(:,idx);
all(a==b,'all')
plot(freq,a)
hold on
plot(freq, b)
%% GF series test
nGF = 2;
GFs = linspace(0.1, 1, nGF);
data = [];

for GF = GFs
    data = cat(2, data, correct_global_vec(gpuData, GF));
end

d = reshape(data, nFreq, nPixel, []); % reshape into (f, pixel, gf)
[normData,mn,mx] = normalize_gpu_data(data, nFreq, nPixel,nGF);
d_ = reshape(normData, nFreq, nPixel, []); % reshape into (f, pixel, gf)

fprintf('<>   GF correction iscorrect: %i\n', all(cdnorm - d(:,:,1) < 0.0001, 'all'))

meanData = squeeze(mean(cdnorm, 2));
mdata = squeeze(mean(d, 2));
mdata_ = squeeze(mean(d_, 2));
% close all
% plot(meanData,'-');
% hold on
% plot(mdata(:,1),'--');
plot(mdata_(:,1),':');
%% normalize GF data
clc
[normData,mn,mx] = normalize_gpu_data(data, nFreq, nPixel,nGF);
d = reshape(data, nFreq, nPixel, []); % reshape into (f, pixel, gf)

%% gaussian fit
clc
close all
hold on
plot(freq, mdata_(:,1),'-');
plot(freq, model_GPU([2.8437,0.0005,1,1,1,1], freq))
%%
preGuess = get_initial_guess(data, freq, 'N14');

%%
initialGuess = [];
for i = 1:6
    initialGuess = cat(2, initialGuess, global_guess(d(:,:,i), freq)); % initial guess for GPUfit
end

%%
close all 
clc
nGF = 3;
GFs = linspace(0.2,0.5,nGF);
fit = fit_resonance(expData, 16, 1, 'type',2, 'globalFraction', GFs,'checkPlot',1);
%%
close all
ax = [];
for i = 1: nGF
    a = subplot(ceil(nGF/5),5,i);
    QDM_figure(100*(fit.resonance(:,:,i)./fit.resonance(:,:,1)-1),...
        'ax',a, 'st',2, ...
        'clim',[-1,1]*1e-2,'title',sprintf('res: GF(%.2f)', GFs(i)));
    ax = [ax a];
end
linkaxes(ax)
%%
figure
ax = [];
for i = 1: nGF
    a = subplot(ceil(nGF/5),5,i);
    QDM_figure(fit.chiSquares(:,:,i)-fit.chiSquares(:,:,1), 'ax',a, 'st',2, ...
        'clim',[-1e-4,1e-4],'title',sprintf('Chi: GF(%.2f)', GFs(i)));
    ax = [ax a];
end
% colormap(turbo(11));
linkaxes(ax)

%%
close all
clc
comp = fit.resonance(:,:,1);
whichGF = ones(size(fit.chiSquares,[1,2])) * GFs(1);
gamma = 0.0028;
zfs = 2.870;

 for i = 2: nGF
    d = fit.resonance(:,:,i);
    idx = fit.chiSquares(:,:,i) < fit.chiSquares(:,:,i-1);
    whichGF(idx) = GFs(i);
    comp(idx) = d(idx);
 end
comp = zfs - comp;
comp = comp/gamma;

f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'title');
ax1 = subplot(1,3,1);
QDM_figure((zfs-fit.resonance(:,:,1))/gamma-9, 'ax', ax1,'title', sprintf('GF:%.4f',GFs(1)),'preThreshold', 1);
ax2 = subplot(1,3,2);
QDM_figure(comp-9, 'ax', ax2, 'title', 'composite','preThreshold', 1);
ax3 = subplot(1,3,3);
imagesc(whichGF);
colormap(ax3, turbo(nGF))
colorbar(ax3)
linkaxes([ax1 ax2 ax3]);
axis equal xy tight
%%
close all
 for i = 2: nGF
    idx = fit.chiSquares(:,:,i) < fit.chiSquares(:,:,i-1);
    comp(idx) = GFs(i);
 end

imagesc(comp);
colormap(jet(nGF))
colorbar
axis equal xy tight
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


function [normData, mn, mx] = normalize_gpu_data(data, nFreq, nPixel, nGF)
    % double normalizes the data (0-1)
    % Parameters
    % ----------
    mn = min(data);
    normData = data - mn;
    mx = max(normData);
    normData = normData ./ mx;
end