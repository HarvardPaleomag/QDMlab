function fit = fit_resonance(expData, binSize, nRes, kwargs)
%[fit] = fit_resonance(expData, binSize, nRes; 'type', 'globalFraction', 'forceGuess', 'checkPlot', 'gaussianFit', 'gaussianFilter', 'smoothDegree', 'diamond', 'slopeCorrection')
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
%  6. if *type* == 1: guess peaks; if *type* == 2: **(a)** get_preGuess -> creates pre guess for single gaussian GPU fit; **(b)** gpu_fit (GAUSS_1D); **(c)** parameters_to_guess calculates the initial guess from the fitted parameters of 6.b.
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
    nRes (1, 1) int16
    kwargs.type (1, 1) {mustBeMember(kwargs.type, [0, 1, 2])} = 2
    kwargs.globalFraction {mustBeNumeric} = 0.5
    kwargs.forceGuess (1, 1) {mustBeMember(kwargs.forceGuess, [1, 0])} = false;
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.gaussianFit (1, 1) {mustBeBoolean(kwargs.gaussianFit)} = false;
    kwargs.gaussianFilter (1, 1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0;
    kwargs.smoothDegree (1, 1) {mustBeNumeric, mustBePositive} = 2
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14'])} = 'N14';
    kwargs.slopeCorrection = false;
end
show_references()

msg = sprintf('--------------------------------------------------------------------');
logMsg('info',msg,1,0);
tStart = tic;

%% fitting related
tolerance = 1e-13;
initialPreGuess = 'none';
fit = struct();

%% check type/diamond combination
if kwargs.type ~= 2 && strcmp(kwargs.diamond, 'N15')
    msg = sprintf('Determining the initial parameters for a fit with this method is not supported for N15 diamonds, yet');
    logMsg('error',msg,1,0);
end

%% data preparation
% this step could easily be skipped, the only thing one needs to figure out
% is how to get the
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes);
[sizeY, sizeX, nFreq] = size(binDataNorm); % binned image x-dimensions
xValues = single(freq');

%% 1. GAUSSIAN BLUR
% default = 0
% gaussian filter on the guesses, can lead to some problems if high gradient
if kwargs.gaussianFilter ~= 0
    msg = sprintf('%i: smoothing data using gaussian blur: %.1f', nRes, gaussianFilter');
    logMsg('info',msg,1,0);
    gFilter = fspecial('gaussian', [20, 20], gaussianFilter);
    binDataNorm = imfilter(binDataNorm, gFilter, 'symmetric', 'conv');
end

% make GPU data
nPixel = sizeX * sizeY; % number of (x,y) pixels
gpudata = reshape(binDataNorm, [nPixel, nFreq]); % make it into 2d matrix
gpudata = transpose(gpudata); %transpose to make it 51 x pixels

msg = sprintf('%i: starting parameter estimation (%s)', nRes, kwargs.diamond');
logMsg('info',msg,1,0);

%% globalFraction subtraction
nGF = size(kwargs.globalFraction,2);
data = [];

for GF = kwargs.globalFraction
    data = cat(2, data, correct_global_vec(gpudata, GF));
end

%% prepare GPUfit data
if kwargs.slopeCorrection
    data = slope_correction(data, freq, kwargs.slopeCorrection);
end

[data, mn, mx] = normalize_gpu_data(data);

initialGuess = global_guess(data, freq); % initial guess for GPUfit

%% GUESS INITIAL FIT PARAMETERS
fit.initialGuess.states = zeros(size(data, 2),1);

%% gloabl guess
switch kwargs.type 
    case 0 % Global
        initialGuess = global_guess(data, freq);
    case 1 % old local
        msg = sprintf('%i: local guess estimation', nRes');
        logMsg('info',msg,1,0);

        %% generating the initialguess for fitting
        % iterate over all pixels
        initialGuess = global_guess(data, freq);
        for i = 1:size(data,2)
            pixelData = squeeze(data(:,i));
            
            % number of pixel in that image
            % Note: for more than one GF factor there will ne n times as
            % many pixels in total
            idx = mod(i, nPixel);
            if idx == 1
                msg = sprintf('%i: Global fraction (%.1f)', nRes', kwargs.globalFraction(1+floor(i/nPixel)));
                logMsg('info',msg,1,0);
            end
            [x,y] = index2xy(idx, [sizeY, sizeX, nFreq]);
            
            [pkVal, pkLoc, fitFlg] = guess_peaks(pixelData, freq, ...
                                    'gaussianFit', kwargs.gaussianFit, ...
                                    'smoothDegree', kwargs.smoothDegree, ...
                                    'forceGuess', kwargs.forceGuess, ...
                                    'pixel', [x,y, nRes]);

            % check if find peaks returned 3 peaks
            % add them to fit.initialGuess.states
            if fitFlg ~= 0 % 1 == gauss 2= global (i.e. local failed)
                fit.initialGuess.states(i) = fitFlg;
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
                initialGuess(:,i) = [resonance, width, contrast, baseline];
            end
        end
        
    case 2 % new local
        %% initial preGuess
        initialPreGuess = get_preGuess(data, freq, kwargs.diamond);
        fit.initialPreGuess = initialPreGuess;

        if strcmp(kwargs.diamond, 'N14')
            % single gaus fit for initial parameters
            model_id = ModelID.GAUSS_1D;

            [initialGuess, states, chiSquares, n_iterations, time] = gpufit(single(data), [], ...
                model_id, initialPreGuess, tolerance, 1000, ...
                [], EstimatorID.LSE, xValues);
            
            initialGuess = parameters_to_guess(initialGuess, kwargs.diamond);
            fit.initialGuess.chi = chiSquares;
            fit.initialGuess.states = states;
            
        elseif strcmp(kwargs.diamond, 'N15')
            msg = sprintf('determining initial guess only from (N14) preInitialGuess'); 
            logMsg('debug',msg,1,0);
            initialGuess = parameters_to_guess(initialPreGuess, kwargs.diamond);
            fit.initialGuess.states = zeros(size(initialGuess));
        end
end

%% local guess -> guess parameter for each pixel

%% GPU pre fits

nBadPixels = numel(nonzeros(fit.initialGuess.states));

if nBadPixels > 1
    if kwargs.type == 2 && strcmp(kwargs.diamond, 'N14')
        msg = sprintf('%i: %i / %i (%i globalFractions) pixels failed the pre-guess. See fits.initialGuess.states', nRes, nBadPixels, numel(fit.initialGuess.states), nGF);
        logMsg('warn',msg,1,0);
    else
        msg = sprintf('%i: %i / %i (%i globalFractions) pixels had to be substituted', nRes, nBadPixels, numel(fit.initialGuess.states), nGF);
        logMsg('warn',msg,1,0);
    end
end
msg = sprintf('%i: initial parameter estimation complete in: %.1f s', nRes, toc(tStart));
logMsg('info',msg,1,0);

%% FINAL GPU FIT
if strcmp(kwargs.diamond, 'N14')
    model_id = ModelID.ESR14N;
elseif strcmp(kwargs.diamond, 'N15')
    model_id = ModelID.ESR15N;
end

max_n_iterations = 1000;

msg = sprintf('%i: starting GPU fit, model: %s', nRes);
logMsg('info',msg,1,0);

% run Gpufit - Res 1
[parameters, states, chiSquares, n_iterations, time] = gpufit(single(data), [], ...
    model_id, single(initialGuess), tolerance, max_n_iterations, [], EstimatorID.LSE, xValues);

% failed fits for pixel with extrem chiSquared or values outside of the
% measurement freq. range
% states(chiSquares > 1e-4) = 5;
% states(parameters(1, :) > max(freq)) = 6;
% states(parameters(1, :) < min(freq)) = 6;

fit = make_fit_struct(fit, initialPreGuess, initialGuess, parameters, states, ...
    chiSquares, n_iterations, nRes, sizeX, sizeY, nPixel, nGF, kwargs.diamond,mn,mx);
fit.freq = freq;
fit.binSize = binSize;

fit.GFs = kwargs.globalFraction;

msg = sprintf('%i: final GPU fitting complete in: %.1f s', nRes, toc(tStart));
logMsg('info',msg,1,0);

if numel(nonzeros(states)) > 0
    badPre = numel(nonzeros(states));
    msg = sprintf('%i: %i / %i (%i globalFractions) pixels failed the final fit. See fit.states!', nRes, badPre, numel(states), nGF);
    logMsg('warn',msg,1,0);
end

if kwargs.checkPlot
    msg = sprintf('close figure to continue');
    logMsg('ATTENTION',msg,1,0);
    fig = gpu_fit_checkPlot(fit, data, freq, binSize, kwargs.diamond, kwargs.globalFraction);
    waitfor(fig)
end

end

%%
function data = slope_correction(data, freq, nPoints)
%[data] = slope_correction(data, freq, nPoints)
% calculates slope between 1st - last pixel and removes this from data
    msg = sprintf('correcting slope of for the initial guess calculation');
    logMsg('debug',msg,1,0);
    
    d1 = mean(data(1:nPoints+1,:));
    dend = mean(data(end+1-nPoints:end,:));
    
    delta = (dend-d1);
    slope = delta/numel(freq);
    
    correction = zeros(size(data));

    for i = 1:numel(freq)
        correction(i,:) = slope*(i-1);
    end

    data = data - correction;
end

function [normData, mn, mx] = normalize_gpu_data(data)
    % double normalizes the data (0-1)
    % Parameters
    % ----------
    mn = min(data);
    normData = data - mn;
    mx = max(normData);
    normData = normData ./ mx;
end

%% fitting helper functions
function initialGuess = get_preGuess(gpudata, freq, diamond)
    %[initialGuess] = get_preGuess(gpudata, freq, diamond)
    initialGuess = zeros(4, size(gpudata, 2), 'single');
    %     n = 1; % cut off outer points
    %     gpudata = gpudata(n:end-n,:);

    % amplitude
    mx = max(gpudata,[], 'omitnan');
    mn = min(gpudata,[], 'omitnan');
    initialGuess(1, :) = -double(abs(((mx - mn)./mx)));

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
    initialGuess(4, :) = double(1);
end

function guess = parameters_to_guess(parameters, diamond)
    %[guess] = parameters_to_guess(parameters, diamond)
    if strcmp(diamond, 'N14')
        guess = zeros(6, size(parameters, 2));
        guess(1, :) = parameters(2, :); % location
        guess(2, :) = 0.0003; % width
        guess(3, :) = abs(parameters(1, :)); % amplitude (contrast)
        guess(4, :) = abs(parameters(1, :)); % amplitude (contrast)
        guess(5, :) = abs(parameters(1, :)); % amplitude (contrast)
        guess(6, :) = parameters(4, :)-1 ; % baseline
        guess = single(guess);

    elseif strcmp(diamond, 'N15')
        guess = zeros(5, size(parameters, 2));
        guess(1, :) = parameters(2, :); % location
        guess(2, :) = 0.0005; % width
        guess(3, :) = -parameters(1, :); % amplitude (contrast)
        guess(4, :) = -parameters(1, :); % amplitude (contrast)
        guess(5, :) = parameters(4, :) - 1; % baseline
        guess = single(guess);
    end
end

function fit = make_fit_struct(fit, preGuess, initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY, nPixel, nGF, diamond, mn, mx)
    %[fit] = make_fit_struct(fit, preGuess, initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY, diamond)
    msg = sprintf('%i: reshaping data into (%4i, %1i)', nRes, nPixel, nGF);
    logMsg('debug',msg,1,0);
    
    shape = [sizeY, sizeX, nGF];
    %make parameters matrix into 3d matrix with x pixels, y pixels, and parameters
    %     fit.preGuess = reshape(preGuess, [], sizeY, sizeX); % initial guess
    fit.initialGuess.parameters = double(reshape(initialGuess, [], sizeY, sizeX, nGF)); % initial guess
    fit.initialGuess.p = double(initialGuess);

    fit.parameters = double(reshape(parameters, [], sizeY, sizeX, nGF)); % fitted parameters
    fit.p = double(parameters);

    % matricies with 2 dimensions for x and y pixels:

    fit.resonance = reshape(fit.parameters(1, :), shape);
    fit.width = reshape(fit.parameters(2, :), shape);
    fit.contrastA = reshape(fit.parameters(3, :), shape);
    fit.contrastB = reshape(fit.parameters(3, :), shape);

    % check for diamond type: N14 has 6 parameters, N15 has only 5
    nParams = size(initialGuess, 1);

    if nParams == 6 % for N14 diamonds
        fit.contrastC = reshape(fit.parameters(5, :), shape);
        fit.baseline = 1 + reshape(fit.parameters(6, :), shape);
    else % for N15 diamonds
        fit.baseline = 1 + reshape(fit.parameters(5, :), shape);
    end

    fit.states = reshape(states, shape);
    fit.chiSquares = double(reshape(chiSquares, shape));
    fit.n_iterations = reshape(n_iterations,shape);
    fit.nRes = nRes;

    if ~strcmp(preGuess, 'none')
        fit.pg = double(parameters_to_guess(preGuess, diamond));
    else
        fit.pg = preGuess;
    end

end
