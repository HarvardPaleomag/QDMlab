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
    nRes (1, 1) int16
    kwargs.type (1, 1) {mustBeMember(kwargs.type, [0, 1, 2])} = 2
    kwargs.globalFraction (1, 1) {mustBeNumeric} = 0.5
    kwargs.forceGuess (1, 1) {mustBeMember(kwargs.forceGuess, [1, 0])} = false;
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.gaussianFit (1, 1) {mustBeBoolean(kwargs.gaussianFit)} = false;
    kwargs.gaussianFilter (1, 1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0;
    kwargs.smoothDegree (1, 1) {mustBeNumeric, mustBePositive} = 2
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14', 'DAC'])} = 'N14';
    kwargs.slopeCorrection = false;
    kwargs.crop = 'none'

end
show_references()

msg = sprintf('--------------------------------------------------------------------');
logMsg('info',msg,1,0);
tStart = tic;

% preallocate fit
fit = struct();

%% check type/diamond combination
if kwargs.type ~= 2 && strcmp(kwargs.diamond, 'N15')
    msg = sprintf('Determining the initial parameters for a fit with this method is not supported for N15 diamonds, yet');
    logMsg('error',msg,1,0);
end

%% data preparation
% this step could easily be skipped, the only thing one needs to figure out
% is how to get the
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes, 'crop', kwargs.crop);


sizeX = size(binDataNorm, 2); % binned image x-dimensions
sizeY = size(binDataNorm, 1); % binned image y-dimensions

%% 1. GAUSSIAN BLUR
% default = 0
% gaussian filter on the guesses, can lead to some problems if high gradient
if kwargs.gaussianFilter ~= 0
    msg = sprintf('%i: smoothing data using gaussian blur: %.1f', nRes, gaussianFilter');
    logMsg('info',msg,1,0);
    gFilter = fspecial('gaussian', [20, 20], gaussianFilter);
    binDataNorm = imfilter(binDataNorm, gFilter, 'symmetric', 'conv');
end

msg = sprintf('%i: starting parameter estimation (%s)', nRes, kwargs.diamond');
logMsg('info',msg,1,0);

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

if kwargs.slopeCorrection
    gpudata = slope_correction(gpudata, freq, kwargs.slopeCorrection);
end

xValues = single(freq');

%% GUESS INITIAL FIT PARAMETERS

%% gloabl guess
if kwargs.type == 0 % reshape into [6 x numpoints]
    initialGuess = reshape(initialGuess, [imgPts, 6]);
    initialGuess = transpose(initialGuess);
end

%% local guess -> guess parameter for each pixel

%% output setup
fit.initialGuess.states = zeros(size(binDataNorm, [1 2]));

if kwargs.type == 1 %% old local/gaussian guess
    msg = sprintf('%i: local guess estimation', nRes');
    logMsg('info',msg,1,0);
    
    sizeX = size(binDataNorm, 2); % binned image x-dimensions
    sizeY = size(binDataNorm, 1); % binned image y-dimensions

    %% generating the initialguess for fitting
    % iterate over all pixels with index (x,y)
    for x = 1:sizeX
        for y = 1:sizeY
            pixelData = squeeze(binDataNorm(y, x, :));

            [pkVal, pkLoc, fitFlg] = guess_peaks(pixelData, meanData, freq, ...
                'smoothDegree', kwargs.smoothDegree, ...
                'forceGuess', kwargs.forceGuess, ...
                'gaussianFit', kwargs.gaussianFit, ...
                'pixel', [y, x, nRes]);
            
            % check if find peaks returned 3 peaks
            % add them to fit.initialGuess.states
            if fitFlg ~= 0 % 1 == gauss 2= global (i.e. local failed)
                fit.initialGuess.states(y,x) = fitFlg;
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

%% fitting related
tolerance = 1e-13;
initialPreGuess = 'none';

if kwargs.type == 2
    %% initial preGuess
    initialPreGuess = get_initial_guess(gpudata, freq, kwargs.diamond);
    fit.initialPreGuess = initialPreGuess;
    
    if find(strcmp(kwargs.diamond, {'N14', 'DAC'}))
        % single gaus fit for initial parameters
        model_id = ModelID.GAUSS_1D;
        [initialGuess, states, chiSquares, n_iterations, time] = gpufit(gpudata, [], ...
            model_id, initialPreGuess, tolerance, 1000, ...
            [], EstimatorID.MLE, xValues);
        initialGuess = parameters_to_guess(initialGuess, kwargs.diamond);
        fit.initialGuess.chi = chiSquares;
        fit.initialGuess.states = states;
    elseif find(strcmp(kwargs.diamond, {'N15'}))
        msg = sprintf('determining initial guess only from (N14) preInitialGuess'); 
        logMsg('debug',msg,1,0);
        % Note: DAC initial guess is different than N15
        % 
        initialGuess = parameters_to_guess(initialPreGuess, kwargs.diamond);
        fit.initialGuess.states = zeros(size(initialGuess));
    end
end

nBadPixels = numel(nonzeros(fit.initialGuess.states));

if nBadPixels > 1
    if kwargs.type == 2 && strcmp(kwargs.diamond, 'N14')
        msg = sprintf('%i: %i / %i pixels failed the pre-guess. See fits.states', nRes, nBadPixels, imgPts);
        logMsg('warn',msg,1,0);
    else
        msg = sprintf('%i: %i / %i pixels had to be substituted', nRes, nBadPixels, imgPts);
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
elseif strcmp(kwargs.diamond, 'DAC')
    model_id = ModelID.GAUSS_1D;
end

max_n_iterations = 1000;

msg = sprintf('%i: starting GPU fit, model: %s', nRes);
logMsg('info',msg,1,0);

% run Gpufit
[parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
    model_id, initialGuess, tolerance, max_n_iterations, [], EstimatorID.MLE, xValues);

% failed fits for pixel with extrem chiSquared or values outside of the
% measurement freq. range
states(chiSquares > 1e-4) = 5;
% states(parameters(1, :) > max(freq)) = 6;
% states(parameters(1, :) < min(freq)) = 6;

fit = make_fit_struct(fit, initialPreGuess, initialGuess, parameters, states, ...
    chiSquares, n_iterations, nRes, sizeX, sizeY, kwargs.diamond);
fit.freq = freq;
fit.binSize = binSize;

msg = sprintf('%i: final GPU fitting complete in: %.1f s', nRes, toc(tStart));
logMsg('info',msg,1,0);

if numel(nonzeros(states)) > 0
    badPre = numel(nonzeros(states));
    msg = sprintf('%i: %i / %i pixels failed the final fit. See fit.states!', nRes, badPre, imgPts);
    logMsg('warn',msg,1,0);
end

if kwargs.checkPlot
    msg = sprintf('close figure to continue');
    logMsg('ATTENTION',msg,1,0);
    fig = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, kwargs.diamond);
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

%% fitting helper functions
function initialGuess = get_initial_guess(gpudata, freq, diamond)
%[initialGuess] = get_initial_guess(gpudata, freq, diamond)
initialGuess = zeros(4, size(gpudata, 2), 'single');


% amplitude
mx = nanmax(gpudata);
mn = nanmin(gpudata);
initialGuess(1, :) = -abs(((mx - mn)./mx));

% center frequency
[~, idx] = sort(gpudata);
l = 7; % lowest n values
mxidx = max(idx(1:l, :));
mnidx = min(idx(1:l, :));

if strcmp(diamond, 'N15') |  strcmp(diamond, 'DAC')
    cIdx = int16((mxidx+mnidx)/2);
else
    cIdx = int16(mean(cat(1, mxidx, mnidx)));
end

center = freq(cIdx);
initialGuess(2, :) = center;

% width
if strcmp(diamond, 'DAC')
    initialGuess(3, :) = 0.005;
else
    initialGuess(3, :) = 0.0005;
end
    

% offset
initialGuess(4, :) = mx;
end

function guess = parameters_to_guess(parameters, diamond)
%[guess] = parameters_to_guess(parameters, diamond)
if strcmp(diamond, 'N14')
    guess = zeros(6, size(parameters, 2));
    guess(1, :) = parameters(2, :); % location
    guess(2, :) = 0.0002; % width
    guess(3, :) = abs(parameters(1, :)); % amplitude (contrast)
    guess(4, :) = abs(parameters(1, :)); % amplitude (contrast)
    guess(5, :) = abs(parameters(1, :)); % amplitude (contrast)
    guess(6, :) = parameters(4, :) - 1; % baseline
    guess = single(guess);
    
elseif find(strcmp(diamond, {'N15'}))
    guess = zeros(5, size(parameters, 2));
    guess(1, :) = parameters(2, :); % location
    guess(2, :) = 0.0002; % width
    guess(3, :) = -parameters(1, :); % amplitude (contrast)
    guess(4, :) = -parameters(1, :); % amplitude (contrast)
    guess(5, :) = parameters(4, :) - 1; % baseline
    guess = single(guess);
    
elseif strcmp(diamond, 'DAC')
    guess = parameters;
end
end

function fit = make_fit_struct(fit, preGuess, initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY, diamond)
%[fit] = make_fit_struct(fit, preGuess, initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY, diamond)
msg = sprintf('%i: reshaping data into (%4i, %4i)', nRes, sizeY, sizeX);
logMsg('debug',msg,1,0);

%make parameters matrix into 3d matrix with x pixels, y pixels, and parameters
%     fit.preGuess = reshape(preGuess, [], sizeY, sizeX); % initial guess
fit.initialGuess.parameters = double(reshape(initialGuess, [], sizeY, sizeX)); % initial guess
fit.initialGuess.p = double(initialGuess);

fit.parameters = double(reshape(parameters, [], sizeY, sizeX)); % fitted parameters
fit.p = double(parameters);

if ~ strcmp(diamond, 'DAC')
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
else
    fit.resonance = squeeze(fit.parameters(2, :, :));
    fit.width = squeeze(fit.parameters(3, :, :));
    fit.contrastA = -squeeze(fit.parameters(1, :, :));
    fit.baseline = squeeze(fit.parameters(4, :, :));
end


fit.states = reshape(states, [sizeY, sizeX]);
fit.chiSquares = double(reshape(chiSquares, [sizeY, sizeX]));
fit.n_iterations = reshape(n_iterations, [sizeY, sizeX]);
fit.nRes = nRes;

if ~strcmp(preGuess, 'none')
    fit.pg = double(parameters_to_guess(preGuess, diamond));
else
    fit.pg = preGuess;
end

end
