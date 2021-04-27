function [fit, initialGuess] = fit_resonance(expData, binSize, nRes, kwargs)
% fits a single resonance frequency (i.e. low/high frequency range) of
% either positive or negative field.
%
% hint
% ----
%  This is what the function does:
%
%  1. :code:`prepare_raw_data`
%  2. gaussian filter (:code:`imfilter`) if `gaussianFilter` == 1
%  3. Global illumination correction with :code:`correct_global`
%  4. global_guess if `type` ~= 2
%  5. reshapes data for gpufit -> all pixels in a row
%  6. if *type* == 1: guess peaks; if *type* == 2: **(a)** get_initial_guess -> creates pre guess for single gaussian GPU fit; **(b)** gpu_fit (GAUSS_1D); **(c)** parameters_to_guess calculates the initial guess from the fitted parameters of 6.b.
%  7. :code:`gpu_fit` calculates lorentzian fits
%  8. :code:`reshape_fits` creates the (y,x) sized array out of the fitted parameters
%  9. *checkPlot* of fits -> needs to be closed to proceed
%
%
% Parameters
% ----------
%     expData: struct
%         Data of :code:`load(run0000n.mat)`
%     binSize: int
%         binning size (can be 1)
%     nRes: int
%         number of resonance. Low frequencies = 1, High frequencies = 2
%     type: int (2)
%         type of initial guess:
%         1: global
%         2: local
%         3: gaussian
%     globalFraction: double (0.5)
%         Ammount of global signal to be corrected for (see. correct_global)
%     forceGuess: int (0)
%         Used for forcing a guess (NOT IMPLEMENTED)
%     checkPlot: int (0)
%         Creates an interactive plot to check the fits
%     gaussianFit: int (0)
%         In case the :code:`type = local` and the MATLAB function find_peaks does
%         not find 3 peaks.
%         **if 0**: the global guess will be used for that pixel
%         **if 1**: a gaussian fit is used to find peak positions
%     gaussianFilter: int (0)
%         Determines if a gaussian filter is applied before fitting
%     smoothDegree: int (2)
%         The ammount of smoothing if gaussianFilter == 1
%     diamond: str (N14)
%         The type of diamond. Choses the type of fitting.
%
%  state definitions: CONVERGED = 0, MAX_ITERATION = 1,
%                     SINGULAR_HESSIAN = 2, NEG_CURVATURE_MLE = 3,
%                     GPU_NOT_READY = 4,
%                     X2 > 1e-4 = 5, fRes > +- max(frequency) = 6

arguments
    expData struct
    binSize double
    nRes(1, 1) int16
    % keyword arguments
    kwargs.type(1, 1) {mustBeMember(kwargs.type, [0, 1, 2])} = 2
    kwargs.globalFraction(1, 1) {mustBeNumeric} = 0.5
    kwargs.forceGuess(1, 1) {mustBeMember(kwargs.forceGuess, [1, 0])} = 0
    kwargs.checkPlot(1, 1) {mustBeBoolean(kwargs.checkPlot)} = 0
    kwargs.gaussianFit(1, 1) {mustBeBoolean(kwargs.gaussianFit)} = 0
    kwargs.gaussianFilter(1, 1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree(1, 1) {mustBeNumeric, mustBePositive} = 2
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14'])} = 'N14';
end

disp('<> --------------------------------------------------------------------')
tStart = tic;

%% check type/diamond combination
if kwargs.type ~= 2 && strcmp(kwargs.diamond, 'N15')
    disp('<>   ERROR: Determining the initial parameters for a fit with this method is not supported for N15 diamonds, yet')
end

dataStack = expData.(sprintf('imgStack%i', nRes));

%% output setup
pixelAlerts = struct('x', {}, 'y', {}, 'nRes', {}, 'fitFlg', {}, 'fName', {});

%% data preparation
% this step could easily be skipped, the only thing one needs to figure out
% is how to get the
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes);

sizeX = size(binDataNorm, 2); % binned image x-dimensions
sizeY = size(binDataNorm, 1); % binned image y-dimensions

%% 1. GAUSSIAN BLUR
% default = 0
% gaussian filter on the guesses, can lead to some problems if high gradient
if kwargs.gaussianFilter ~= 0
    fprintf('<>   %i: smoothing data using gaussian blur: %.1f\n', nRes, gaussianFilter)
    gFilter = fspecial('gaussian', [20, 20], gaussianFilter);
    binDataNorm = imfilter(binDataNorm, gFilter, 'symmetric', 'conv');
end

fprintf('<>   %i: starting parameter estimation (%s)\n', nRes, kwargs.diamond);

%% global spectra subtraction
binDataNorm = correct_global(binDataNorm, kwargs.globalFraction);

%% first determine global guess
meanData = squeeze(mean(binDataNorm, [1, 2]));

if kwargs.type ~= 2
    initialGuess = global_guess(binDataNorm, freq); % initial guess for GPUfit
end

%% prepare GPUfit data
sweepLength = size(freq, 2);
imgPts = sizeX * sizeY; % number of (x,y) pixels
gpudata = reshape(binDataNorm, [imgPts, sweepLength]); % make it into 2d matrix
gpudata = transpose(gpudata); %transpose to make it 51 x pixels
gpudata = single(gpudata);
xValues = single(freq');

%% GUESS INITIAL FIT PARAMETERS

%% gloabl guess
if kwargs.type == 0 % reshape into [6 x numpoints]
    initialGuess = reshape(initialGuess, [imgPts, 6]);
    initialGuess = transpose(initialGuess);
end

%% local guess -> guess parameter for each pixel
if kwargs.type == 1 %% old local/gaussian guess
    fprintf('<>   %i: local guess estimation\n', nRes);

    sizeX = size(binData, 2); % binned image x-dimensions
    sizeY = size(binData, 1); % binned image y-dimensions

    %% generating the initialguess for fitting
    % iterate over all pixels with index (x,y)
    for x = 1:sizeX
        for y = 1:sizeY
            pixelData = squeeze(binDataNorm(y, x, :));

            % try getting peak positions from smoothed data
            % LEFT
            [pkVal, pkLoc, fitFlg] = guess_peaks(pixelData, meanData, freq, ...
                'smoothDegree', kwargs.smoothDegree, ...
                'forceGuess', kwargs.forceGuess, ...
                'gaussianFit', kwargs.gaussianFit, ...
                'pixel', [y, x, nRes]);
            % check if find peaks returned 3 peaks
            % add them to pixelAlerts
            if fitFlg == 1 || fitFlg == 2 % 1 == gauss 2= global (i.e. local failed)
                pixelAlerts(end+1).x = x; % new entry
                pixelAlerts(size(pixelAlerts, 2)).y = y;
                pixelAlerts(size(pixelAlerts, 2)).nRes = nRes;
                pixelAlerts(size(pixelAlerts, 2)).fitFlg = fitFlg;
            end

            % replace guess with new guess if fitFlg is not global
            % (i.e. fitFlg = 2)
            if fitFlg ~= 2
                % if it returns 3 -> replace the global guess with
                % the local
                resonance = (pkLoc(1) + pkLoc(2) + pkLoc(3)) / 3; % in GHz
                width = 0.0005;
                contrast = (mean(pixelData(1:10)) + pkVal - 1)';
                baseline = mean(pixelData(1:10)) - 1;
                initialGuess(y, x, :) = [resonance, width, contrast, baseline];
            end
        end
    end
    initialGuess = reshape(initialGuess, [imgPts, 6]);
    initialGuess = transpose(initialGuess);
end

%% GPU pre fits

%% fittiing related
tolerance = 1e-10;

if kwargs.type == 2
    % initial parameters
    initialPreGuess = get_initial_guess(gpudata, freq, kwargs.diamond);

    %initiate badPixels structure
    pixelAlerts = struct();
    if strcmp(kwargs.diamond, 'N14')
        % single gaus fit for initial parameters
        model_id = ModelID.GAUSS_1D;
        [initialGuess, states, chiSquares, n_iterations, time] = gpufit(gpudata, [], ...
            model_id, initialPreGuess, tolerance, 1000, ...
            [], EstimatorID.MLE, xValues);
        initialGuess = parameters_to_guess(initialGuess, kwargs.diamond);
        pixelAlerts.chi = chiSquares;
        pixelAlerts.state = states;
    elseif strcmp(kwargs.diamond, 'N15')
        initialGuess = parameters_to_guess(initialPreGuess, kwargs.diamond);
        pixelAlerts.state = zeros(size(initialPreGuess));
    end
    pixelAlerts.initialPreGuess = initialPreGuess;
end

nBadPixels = numel(nonzeros(pixelAlerts.state));

if nBadPixels > 1
    if kwargs.type == 2 && strcmp(kwargs.diamond, 'N14')
        fprintf('<>      WARNING: %i / %i pixels failed the pre-guess. See fits.states\n', nBadPixels, imgPts);
    else
        fprintf('<>      WARNING: %i / %i pixels had to be substituted\n', nBadPixels, imgPts);
    end
end

fprintf('<>      INFO: initial parameter estimation complete in: %.1f s\n', toc(tStart)');

%% FINAL GPU FIT
if strcmp(kwargs.diamond, 'N14')
    model_id = ModelID.ESR3RT;
elseif strcmp(kwargs.diamond, 'N15')
    model_id = ModelID.ESR15N;
end

max_n_iterations = 1000;

fprintf('<>   %i: starting GPU fit, model: %s\n', nRes);
% run Gpufit - Res 1
[parameters, states, chiSquares, n_iterations, time] = gpufit(gpudata, [], ...
    model_id, initialGuess, tolerance, max_n_iterations, [], EstimatorID.LSE, xValues);
% make parameters double again
parameters = double(parameters);

% failed fits for pixel with extrem chiSquared or values outside of the
% measurement freq. range
states(chiSquares > 1e-4) = 5;
states(parameters(1, :) > max(freq)) = 6;
states(parameters(1, :) < min(freq)) = 6;

fit = make_fit_struct(initialPreGuess, initialGuess, parameters, states, ...
    chiSquares, n_iterations, nRes, sizeX, sizeY, kwargs.diamond);
fit.freq = freq;
fit.binSize = binSize;

fprintf('<>      INFO: final GPU fitting complete in: %.1f s\n', toc(tStart)');

if numel(nonzeros(states)) > 0
    badPre = numel(nonzeros(states));
    fprintf('<>      WARNING: %i / %i pixels failed the final fit. See fit.states!\n', badPre, imgPts);
end

if kwargs.checkPlot
    fprintf('<>>>>>> INFO: close figure to continue\n');
    fig = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, kwargs.diamond);
    waitfor(fig)
end

end

%% fitting helper functions
function initialGuess = get_initial_guess(gpudata, freq, diamond)
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
initialGuess(3, :) = 0.0015;
% offset
initialGuess(4, :) = mx;
end

function guess = parameters_to_guess(parameters, diamond)
if strcmp(diamond, 'N14')
    guess = zeros(6, size(parameters, 2));
    guess(1, :) = parameters(2, :); % location
    guess(2, :) = 0.0004; % width
    guess(3, :) = -parameters(1, :); % amplitude (contrast)
    guess(4, :) = -parameters(1, :); % amplitude (contrast)
    guess(5, :) = -parameters(1, :); % amplitude (contrast)
    guess(6, :) = parameters(4, :) - 1; % baseline
    guess = single(guess);
elseif strcmp(diamond, 'N15')
    guess = zeros(5, size(parameters, 2));
    guess(1, :) = parameters(2, :); % location
    guess(2, :) = 0.0004; % width
    guess(3, :) = -parameters(1, :); % amplitude (contrast)
    guess(4, :) = -parameters(1, :); % amplitude (contrast)
    guess(5, :) = parameters(4, :) - 1; % baseline
    guess = single(guess);
end
end

function fit = make_fit_struct(preGuess, initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY, diamond)

% initialize struct
fit = struct();

fprintf('<>   %i: INFO: reshaping data into (%4i, %4i)\n', nRes, sizeY, sizeX);

%make parameters matrix into 3d matrix with x pixels, y pixels, and parameters
%     fit.preGuess = reshape(preGuess, [], sizeY, sizeX); % initial guess
fit.initialGuess = reshape(initialGuess, [], sizeY, sizeX); % initial guess
fit.parameters = reshape(parameters, [], sizeY, sizeX); % fitted parameters

% matricies with 2 dimensions for x and y pixels:

fit.resonance = squeeze(fit.parameters(1, :, :));
fit.width = squeeze(fit.parameters(2, :, :));
fit.contrastA = squeeze(fit.parameters(3, :, :));
fit.contrastB = squeeze(fit.parameters(4, :, :));

% check for diamond type: N14 has 6 parameters, N15 has only 5
nParams = size(initialGuess, 1);

if nParams == 6 % for N14 diamonds
    fit.contrastC = squeeze(fit.parameters(5, :, :));
    fit.baseline = squeeze(fit.parameters(6, :, :)+1);
else % for N15 diamonds
    fit.baseline = squeeze(fit.parameters(5, :, :)+1);
end

fit.states = reshape(states, [sizeY, sizeX]);
fit.chiSquares = reshape(chiSquares, [sizeY, sizeX]);
fit.n_iterations = reshape(n_iterations, [sizeY, sizeX]);
fit.nRes = nRes;

fit.p = parameters;
fit.g = initialGuess;
fit.pg = parameters_to_guess(preGuess, diamond);
end
