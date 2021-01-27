function [fit, initialGuess, badPixels] = fit_resonance(expData, binSize, nRes, kwargs)
% fits a single resonance frequency (i.e. low/high frequency range) of
% either positive or negative field
%
%   1. prepare_raw_data
%   2. gaussian filter (imfilter) if `gaussienFilter` == 1
%   3. correction of global data
%   4. global_guess if type ~= 2
%   5. reshapes data for gpufit -> all pixels in a row
%   6. type = 1: guess peaks
%      type = 2: a. get_initial_guess -> creates pre guess for single gaussian GPU fit
%                b. gpu_fit (GAUSS_1D)
%                c. parameters_to_guess calculates the initial guess from the fitted parameters of 6.b.
%   7. gpu_fit calculates lorentzian fits %todo add option for N15
%   8. reshape_fits creates the (y,x) sized array out of the fitted parameters
%  (9.) checkPlot of fits -> needs to be closed to proceed
%
% Parameters
% ----------
%     required
%     ========
%     expData: struct
%         Data of load(run0000n.mat)
%     binSize: int
%         binning size (can be 1)
%     nRes: int
%         number of resonance. Low frequencies = 1, High frequencies = 2
%     keyword
%     =======
%     type: int (2)
%         type of initial guess:
%         1: global
%         2: local
%         3: gaussian
%     globalFraction: double (0.5)
%         Ammount of global signal to be corrected for (see. correct_global)
%     forceGuess: int (0)
%         Used for forcing a guess (NOT IMP{LEMENTED) %todo
%     checkPlot: int (0)
%         Creates an interactive plot to check the fits
%     gaussianFit: int (0)
%         In case the type = local and the MATLAB function find_peaks does
%         not find 3 peaks:
%             if 0: the global guess will be used for that pixel
%             if 1: a gaussian fit is used to find peak positions
%     gaussianFilter: int (0)
%         Determines if a gaussian filter is applied before fitting
%     smoothDegree: int (2)
%         The ammount of smoothing if gaussianFilter == 1
%     nucSpinPol: int (0)
%         Not quite sure this is a remnant of the prev. code
                       
arguments
    expData struct
    binSize double
    nRes (1,1) int16
    % keyword arguments
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 2
    kwargs.globalFraction (1,1) {mustBeNumeric} = 0.5
    kwargs.forceGuess (1,1) {mustBeMember(kwargs.forceGuess, [1, 0])} = 0
    kwargs.checkPlot (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0
    kwargs.gaussianFit (1,1) {mustBeMember(kwargs.gaussianFit, [1, 0])} = 0
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.nucSpinPol (1,1) {mustBeMember(kwargs.nucSpinPol, [1, 0])} = 0
end

disp('<> --------------------------------------------------------------------')
tStart = tic;

fit = struct();
dataStack = expData.(sprintf('imgStack%i',nRes));

%% output setup
badPixels = struct('x',{},'y',{},'nRes', {}, 'fitFlg', {}, 'fName',{});

%% fittiing related
tolerance = 1e-20;
max_n_iterations = 100;

%% data preparation
% this step could easily be skipped, the only thing one needs to figure out
% is how to get the 
[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes);

sizeX = size(binDataNorm,2); % binned image x-dimensions
sizeY = size(binDataNorm,1); % binned image y-dimensions

%% 1. GAUSSIAN BLUR
% default = 0
% gaussian filter on the guesses, can lead to some problems if high gradient
if kwargs.gaussianFilter ~= 0
    fprintf('<>   %i: smoothing data using gaussian blur: %.1f\n', nRes, gaussianFilter)
    gFilter = fspecial('gaussian',[20,20], gaussianFilter);
    binDataNorm = imfilter(binDataNorm, gFilter, 'symmetric', 'conv');
end

fprintf('<>   %i: starting parameter estimation\n', nRes);

%% global spectra subtraction
binDataNorm = correct_global(binDataNorm, kwargs.globalFraction);


%% first determine global guess
meanData = squeeze(mean(binDataNorm, [1,2]));

if kwargs.type ~= 2
    initialGuess = global_guess(binDataNorm, freq); % initial guess for GPUfit
end


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
    fprintf('<>   %i: local guess estimation\n', nRes);

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
                                                    'pixel', [y x nRes]);
                % check if find peaks returned 3 peaks
                % add them to badpixels
                if fitFlg == 1 | fitFlg == 2 % 1 == gauss 2= global (i.e. local failed)
                    badPixels(end+1).x = x; % new entry
                    badPixels(size(badPixels, 2)).y = y;
                    badPixels(size(badPixels, 2)).nRes = nRes;
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
    
    % initial parameters
    initialPreGuess = get_initial_guess(gpudata, freq);  
    
    % single gaus fit for initial parameters
    model_id = ModelID.GAUSS_1D;
    [initialGuess, states, chiSquares, n_iterations, time] = gpufit(gpudata, [], ...
                       model_id, initialPreGuess, tolerance, 100, ...
                       [], EstimatorID.MLE, xValues);
    initialGuess = parameters_to_guess(initialGuess);

    badPixels = struct();
    badPixels.initialPreGuess = initialPreGuess;
    badPixels.chi = chiSquares;
    badPixels.state = states;

end

if size(badPixels,2) > 0
    s = size(binDataNorm);
    if kwargs.type == 2
        badPre = numel(nonzeros(badPixels.state));
        fprintf('<>      WARNING: %i / %i pixels failed the pre-guess. See badPixels.states\n', badPre, s(1)*s(2));
    else
        fprintf('<>      WARNING: %i / %i pixels had to be substituted\n', size(badPixels,2), s(1)*s(2));
    end
end

fprintf('<>      INFO: initial parameter estimation complete in: %.1f s\n', toc(tStart)');

% final GPU fits
model_id = ModelID.ESR3RT;
max_n_iterations = 100;

%% FINAL GPU FIT

fprintf('<>   %i: starting GPU fit\n', nRes);
% run Gpufit - Res 1
[parameters, states, chiSquares, n_iterations, time] = gpufit(gpudata, [], ...
    model_id, initialGuess, tolerance, max_n_iterations, [], EstimatorID.LSE, xValues);

fit = reshape_fits(initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY);
fit.freq = freq;
fit.binSize = binSize;

fprintf('<>      INFO: final GPU fitting complete in: %.1f s\n', toc(tStart)');

if numel(nonzeros(states)) > 0
    badPre =  numel(nonzeros(states));
    fprintf('<>      WARNING: %i / %i pixels failed the final fit. See fit.states!\n', badPre, s(1)*s(2));
end

if kwargs.checkPlot
    fprintf('<>>>>>> INFO: close figure to continue\n');
    fig = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize);
    waitfor(fig)
end

end

%% fitting helper functions
function initialGuess =  get_initial_guess(gpudata, freq)
    initialGuess = zeros(4, size(gpudata,2), 'single');
    n = 5; % cut off outer points
    gpudata = gpudata(n:end-n,:);
    
    % amplitude
    mx = nanmax(gpudata);
    mn = nanmin(gpudata);
    initialGuess(1,:) = -2*((mx-mn)./mx);
    
    % center frequency
    l = 10; % lowest n values
    [~, idx] = sort(gpudata);
    idx = int16(median(idx(1:l,:)));
    center = zeros(1, numel(idx));
    
    parfor i = 1:numel(idx)
        center(i) = freq(idx(i));
    end
    
    initialGuess(2,:) = center;
    % width
    initialGuess(3,:) = 0.003;
    % offset 
    initialGuess(4,:) = 1.002;
end
function initialGuess = get_initial_guess_OLD(gpudata, freq)
    initialGuess = zeros(4, size(gpudata,2), 'single');
    n = 5;
    freq = freq(n:end-n);
    for i = 1:size(gpudata,2)
        data = gpudata(n:end-n,i);
        mx = nanmax(data);
        mn = nanmin(data);
        
        [~, idx] = sort(data);
        idx = sort(idx(1:10));
        idx = int16(median(idx));
        
        initialGuess(1,i) = -2*(mx-mn)/mx; % amplitude
        initialGuess(2,i) = freq(idx); % center
        initialGuess(3,i) = 0.003; % width
        initialGuess(4,i) = 1.002;%mean(data); % offset 
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

function fit = reshape_fits(initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY)
    
    % initialize struct
    fit = struct(); 
    
    fprintf('<>   %i: INFO: reshaping data into (%4i, %4i)\n', nRes, sizeY, sizeX);


    %make parameters matrix into 3d matrix with x pixels, y pixels, and parameters
    output = zeros(6,sizeY,sizeX);
    output(1,:,:) = reshape(parameters(1,:),[sizeY,sizeX]);
    output(2,:,:) = reshape(parameters(2,:),[sizeY,sizeX]);
    output(3,:,:) = reshape(parameters(3,:),[sizeY,sizeX]);
    output(4,:,:) = reshape(parameters(4,:),[sizeY,sizeX]);
    output(5,:,:) = reshape(parameters(5,:),[sizeY,sizeX]);
    output(6,:,:) = reshape(parameters(6,:),[sizeY,sizeX]);
    fit.parameters = output; 

    ig = zeros(6,sizeY,sizeX);
    ig(1,:,:) = reshape(initialGuess(1,:),[sizeY,sizeX]);
    ig(2,:,:) = reshape(initialGuess(2,:),[sizeY,sizeX]);
    ig(3,:,:) = reshape(initialGuess(3,:),[sizeY,sizeX]);
    ig(4,:,:) = reshape(initialGuess(4,:),[sizeY,sizeX]);
    ig(5,:,:) = reshape(initialGuess(5,:),[sizeY,sizeX]);
    ig(6,:,:) = reshape(initialGuess(6,:),[sizeY,sizeX]);
    fit.initialGuess = ig;

    % matricies with 2 dimensions for x and y pixels:

    fit.resonance = squeeze(output(1,:,:));
    fit.width = squeeze(output(2,:,:));    
    fit.contrastA = squeeze(output(3,:,:));
    fit.contrastB = squeeze(output(4,:,:));
    fit.contrastC = squeeze(output(5,:,:));
    
    fit.baseline = squeeze(output(6,:,:)+1);
    fit.states = reshape(states,[sizeY,sizeX]);
    fit.chiSquares = reshape(chiSquares,[sizeY,sizeX]);
    fit.n_iterations = reshape(n_iterations,[sizeY,sizeX]);
    fit.nRes = nRes;
    
    fit.p = parameters;
    fit.g = initialGuess;
end