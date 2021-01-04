function [initialGuess, badPixels] = fit_resonance(expData, binSize, freq, n, kwargs)

arguments
    expData struct
    binSize double
    freq double
    n (1,1) int16
    % keyword arguments
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 0
    kwargs.globalFraction (1,1) {mustBeNumeric} = 0.5
    kwargs.forceGuess (1,1) {mustBeMember(kwargs.forceGuess, [1, 0])} = false
    kwargs.checkPlot (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = false
    kwargs.gaussianFit (1,1) {mustBeMember(kwargs.gaussianFit, [1, 0])} = false
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.nucSpinPol (1,1) {mustBeMember(kwargs.nucSpinPol, [1, 0])} = false
end

disp('<> --------------------------------------------------------------------')
tStart = tic;

dataStack = expData.(sprintf('imgStack%i',n));

%% output setup
badPixels = struct('x',{},'y',{},'nRes', {}, 'fitFlg', {}, 'fName',{});

%% fittiing related
tolerance = 1e-30;
max_n_iterations = 100;

%% data preparation
% X/Y of unbinned data
% Note: X = COL; Y = ROW -> (y,x) for (row, col) matlab convention
spanXTrans = 1:expData.imgNumCols;
spanYTrans = 1:expData.imgNumRows;

% check for 101 frequencies. File includes imgStack3
if isfield(expData, 'imgStack3')      
    % combine 1&2 or 3&4
    dataStacka = expData.(sprintf('imgStack%i',n)); 
    dataStackb = expData.(sprintf('imgStack%i',n+1));
    dataStack = [dataStacka; dataStackb];
end

data = zeros(expData.imgNumRows, expData.imgNumCols, expData.numFreqs);

% crop
data = data(spanYTrans,spanXTrans,:);

% reshape and transpose each image
for y = 1:expData.numFreqs
    data(:,:,y) = transpose(reshape(dataStack(y, :), [expData.imgNumCols, expData.imgNumRows]));
end

% binning
fprintf('<>   %i: binning data >> binSize = %i\n', n, binSize);

sizeXY = size(BinImage(data(:,:,1),binSize));
binData = zeros(sizeXY(1),sizeXY(2),length(freq));

for y = 1:length(freq)
    binData(:,:,y) = BinImage(data(:,:,y),binSize);
end

sizeX = size(binData,2); % binned image x-dimensions
sizeY = size(binData,1); % binned image y-dimensions

% Correct for severely non-unity baseline by dividing pixelwise by
% average of all frequency points

binDataNorm = zeros(size(binData));
NormalizationFactor = mean(binData,3);    % compute average

for y = 1:length(freq)
    binDataNorm(:,:,y) = binData(:,:,y) ./ NormalizationFactor;
end

%% 1. GAUSSIAN BLUR
% default = 0
% gaussian filter on the guesses, can lead to some problems if high gradient
if kwargs.gaussianFilter ~= 0
    fprintf('<>   %i: smoothing data using gaussian blur: %.1f\n', n, gaussianFilter)
    gFilter = fspecial('gaussian',[20,20], gaussianFilter);
    binDataNorm = imfilter(binDataNorm, gFilter, 'symmetric', 'conv');
end

fprintf('<>   %i: starting parameter estimation\n', n);

%% global spectra subtraction
binDataNorm = correct_global(binDataNorm, kwargs.globalFraction);

%% first determine global guess
meanData = squeeze(mean(mean(binDataNorm,1),2));

initialGuess = global_guess(binDataNorm, freq); % initial guess for GPUfit

%% prepare GPUfit data
sweeplength = size(dataStack,1);
imgpts = sizeX*sizeY; % number of (x,y) pixels
gpudata = reshape(binDataNorm,[imgpts,sweeplength]); % make it into 2d matrix
gpudata = transpose(gpudata); %transpose to make it 51 x pixels
gpudata = single(gpudata);
xValues = single(freq');

%% GUESS INITIAL FIT PARAMETERS
%% gloabl guess
if kwargs.type == 0 % reshape into [6 x numpoints]
    initialGuess = reshape(initialGuess,[imgpts,6]);
    initialGuess = transpose(initialGuess);
end
%% local guess -> guess parameter for each pixel
if kwargs.type == 1 %% old local/gaussian guess 
    fprintf('<>   %i: local guess estimation\n', n);

    sizeX = size(binData,2); % binned image x-dimensions
    sizeY = size(binData,1); % binned image y-dimensions

    %% generating the initialguess for fitting
    % iterate over all pixels with index (x,y)
    for x = 1:sizeX
        for y = 1:sizeY 
            pixelData = squeeze(binDataNorm(y,x,:)); 

            %% Nuclear Spin polarization
            % default = 0
            if kwargs.nucSpinPol
                initialGuess(y,x,:) = guessNucSpinPol(pixelData, freq);    

            %% no Nuclear Spin Polarization    
            else
                % try getting peak positions from smoothed data
                % LEFT
                [pkVal, pkLoc, fitFlg] = guess_peaks(pixelData, meanData, freq, ...
                                                    'smoothDegree', kwargs.smoothDegree, ...
                                                    'forceGuess', kwargs.forceGuess,...
                                                    'gaussianFit', kwargs.gaussianFit, ...
                                                    'pixel', [y x n]);
                % check if find peaks returned 3 peaks
                % add them to badpixels
                if fitFlg == 1 | fitFlg == 2 % 1 == gauss 2= global (i.e. local failed)
                    badPixels(end+1).x = x; % new entry
                    badPixels(size(badPixels, 2)).y = y;
                    badPixels(size(badPixels, 2)).nRes = n;
                    badPixels(size(badPixels, 2)).fitFlg = fitFlg;                    
                end

                % replace guess with new guess if fitFlg is not global
                % (i.e. fitFlg = 2)
                if fitFlg ~= 2
                    % if it returns 3 -> replace the global guess with
                    % the local
                    resonance  = (pkLoc(1)+pkLoc(2)+pkLoc(3))/3;% in GHz
                    width = 0.0005;
                    contrast = (mean(pixelData(1:10)) + pkVal-1)';
                    baseline = mean(pixelData(1:10))-1;
                    initialGuess(y,x,:) = [resonance width contrast baseline];
                end
            end
        end
    end
    initialGuess = reshape(initialGuess,[imgpts,6]);
    initialGuess = transpose(initialGuess);
%     initialGuess(1,:) = initialGuess(1,:)+res1;
end
%% GPU pre fits
if kwargs.type == 2
    model_id = ModelID.GAUSS_1D;
    
    % initial parameters
    initialPreGuess = get_initial_guess(gpudata, freq);
    
    % single gaus fit for initial parameters
    [initialGuess, states, chiSquares, n_iterations, time] = gpufit(gpudata, [], ...
                       model_id, initialPreGuess, tolerance, 10000, ...
                       [], EstimatorID.MLE, xValues);
    initialGuess = parameters_to_guess(initialGuess);
    badPixels = struct();
    badPixels.chi = chiSquares;
    badPixels.state = states;

end

if size(badPixels,2) > 0
    s = size(binDataNorm);
    fprintf('<>      INFO: %i / %i pixels had to be substituted\n', size(badPixels,2), s(1)*s(2));
end

fprintf('<>      INFO: initial parameter estimation complete in: %.1f s\n', toc(tStart)');

%     % final GPU fits
%     model_id = ModelID.ESR3RT;
%     max_n_iterations = 100;
%     fprintf('<>      INFO: initial parameter estimation complete in: %.1f s\n', toc(start)');
%     
%     %% FINAL GPU FIT
% 
%     disp('<>   Starting final GPU fits');
%     % run Gpufit - Res 1
%     [parameters, states, chiSquares, n_iterations, time] = gpufit(gpudata, [], ...
%         model_id, initialGuess1, tolerance, max_n_iterations, [], EstimatorID.LSE, xValues);
end

%% fitting helper functions
function initialGuess = get_initial_guess(gpudata, freq)
    initialGuess = zeros(4, size(gpudata,2), 'single');
    n = 5;
    freq = freq(n:end-n);
    for i = 1:size(gpudata,2)
        data = gpudata(n:end-n,i);
        data = smooth(data, 20);
        mx = max(data);
        mn = min(data);
        initialGuess(1,i) = -(mx-mn)/mx; % amplitude
        initialGuess(2,i) = freq(find(data==mn,1)); %center
        initialGuess(3,i) = 0.002; % width
        initialGuess(4,i) = mean(data(1:10)); % offset -- mean of highest 10 maybe???
    end
end

function guess = parameters_to_guess(parameters)
    guess = zeros(6, size(parameters,2)); 
    guess(1,:) = parameters(2,:); % location
    guess(2,:) = 0.0004; % width
    guess(3,:) = -parameters(1,:); % amplitude (contrast)
    guess(4,:) = -parameters(1,:); % amplitude (contrast)
    guess(5,:) = -parameters(1,:); % amplitude (contrast)
    guess(6,:) = parameters(4,:)-1; % baseline
    guess = single(guess);
end

