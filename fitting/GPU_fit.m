function bad_pixels = GPU_fit(dataFolder, binSize, kwargs)
%{
Parameters
----------
    required
    ========
    dataFolder: str
        location of the data folder
    binSize: int
        number of pixels to be binned into one.
        Uses 'BinImage' function

    optional
    ========
    fielpolarity: int 
        default: 0 
        0: both polarities
    type: int
        default: 0
        0: ONLY global guess parameters
        1: local guess parameters
           use '
        3: manual guess parameters %redundant???
            todo:(set values in guess1 and guess2 matricies)
    forceGuess: bool (false)
    gaussianFit: bool (false)
        Only used if the findpeaks function fails to find 3 peaks.
        if true: uses a gaussian fit to estimate the center peak of the
                 triplet. 
        if false: uses the locations from the global estimation
    gaussianFilter: numeric (0)
        if 0: no filter is applied
        if != 0: applies gaussian filter with a standard deviation of
                'gaussianFilter'. Previous versions of the code used 0.5.
    checkPlot: bool (false)
    plotGuessSpectra: int (1)
        display the spectrum used to establish the initial guesses for diagnostics
    smoothDegree: int(2)
    globalFraction: double (0.5)
        amount of global illumination signal to subtract from data. 
        Calls 'correct_global' function
    save: bool (true)
        if true the results are saved to 'dataFolder'

    nucSpinPol: bool (false)
        this is used for nuclear spin polarization -> NMR. Calls the
        function guessNucSpinPol. Uses code in original state (< Nov 2020).

Notes
-----
dataStack1:     2 dimentional wwith xy pixels, frequencies
data1:          3 dimentional with y pixels, x pixels, frequencies
binData1:       3 dimentional with binned y pixels, binned x pixels, frequencies
binDataNorm1:   3 dimentional with normalized binned x pixels, normalized binned, y pixels, frequencies
   if gaussianFilter: 3 dimentional with gaussian blurred normalized binned y pixels, gaussian blurred normalized binned, x pixels, frequencies
    
%}
% arguments
arguments
    dataFolder char
    binSize double
    kwargs.fieldPolarity (1,1) {mustBeMember(kwargs.fieldPolarity,[0,1,2])} = 0
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 0
    kwargs.globalFraction (1,1) {mustBeNumeric} = 0.5
    kwargs.forceGuess (1,1) {mustBeMember(kwargs.forceGuess, [1, 0])} = 0
    kwargs.checkPlot (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0
    kwargs.plotGuessSpectra (1,1) {mustBeMember(kwargs.plotGuessSpectra, [1, 0])} = 0
    kwargs.gaussianFit (1,1) {mustBeMember(kwargs.gaussianFit, [1, 0])} = 0
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.nucSpinPol (1,1) {mustBeMember(kwargs.nucSpinPol, [1, 0])} = 0
    kwargs.save (1,1) {mustBeMember(kwargs.save, [1, 0])} = 1
end

fieldPolarity = kwargs.fieldPolarity;
type = kwargs.type;

%% global variables
gamma = 0.0028;  % NV gyromagnetic ratio, in GHz / gauss
zfs = 2.870;     % NV zero-field splitting, in GHz
Ahyp = 0.002158; % longitudinal hyperfine for 14N
globalFraction = kwargs.globalFraction;

% smoothing and filter related
gaussianFilter = kwargs.gaussianFilter;
smoothDegree = kwargs.smoothDegree;  % degree of smoothing for 'sgolay' smoothing
nucSpinPolFlg = kwargs.nucSpinPol;  % nuclear spin polarization -> NMR

quadBGsubFlg = 1 ;
LEDimgFlg = 0;
diagplots = 0;

Mag = 15;

%% determine how many files need to be loaded
polarityMap = containers.Map([0,1,2,4],{[1 2];[1 1];[2 2];[1 4]});
startEnd = polarityMap(fieldPolarity);
startN = startEnd(1);
endN = startEnd(2);
%%

disp(['<>   WORKING DIR: << ' dataFolder ' >>']);
headerFiles = dir(fullfile(dataFolder,'*_header.txt'));
dataFiles = dir(fullfile(dataFolder,'run_0000*.mat'));
% filter files that are not run0000.mat
idx = ~cellfun('isempty', regexpi({dataFiles.name}, 'run_[0-9]{5}\.mat$','match'));

dataFiles = dataFiles(idx);

%pixels where guess_peak had to use global guess
badPixels = [];

%% GUESS PARAMETER ESTIMATION
for fileNum=startN:1:endN
    start = tic; % for timing 
    
    %%% select header and data file
    dataFile = dataFiles(fileNum).name;
    badPixels1 = {}; % for testing
    badPixels2 = {}; % for testing

    %%%
    LEDimgFile = 'laser.csv';%'LED_beforerun.csv';
    
    disp(['<>   loading data file:  ' fullfile(dataFolder, dataFile)]);
    expData = load(fullfile(dataFolder, dataFile));
    disp(['<>     loading complete']);
    
    SpanXTrans = 1:expData.imgNumCols;
    SpanYTrans = 1:expData.imgNumRows;
    
    if LEDimgFlg==1
        transImg = load(fullfile(dataFolder,LEDimgFile));
        transImg = transImg(SpanYTrans,SpanXTrans);
        transImgBin = BinImage(transImg, binSize);
        transImgBinUint = im2uint8forExportDG(transImgBin,min(min(transImgBin)), max(max(transImgBin)) );
    end

    % define variables from expData
    Freqs1 = expData.freqList(1:expData.numFreqs) / 1E9;   %everything is in GHz
    Freqs2 = expData.freqList(1+expData.numFreqs : 2*expData.numFreqs) / 1E9;
            
    dataStack1 = expData.imgStack1;
    dataStack2 = expData.imgStack2;
    
    % check for 101 frequencies. File includes imgStack3
    if isfield(expData, 'imgStack3')      
        dataStack1a = expData.imgStack1;
        dataStack2a = expData.imgStack3;
        dataStack1b = expData.imgStack2;
        dataStack2b = expData.imgStack4;
        
        dataStack1 = [dataStack1a; dataStack1b];
        dataStack2 = [dataStack2a; dataStack2b];
    end
    
    data1 = zeros(expData.imgNumRows, expData.imgNumCols, expData.numFreqs);
    data2 = zeros(expData.imgNumRows, expData.imgNumCols, expData.numFreqs);
    
    % crop
    data1 = data1(SpanYTrans,SpanXTrans,:);
    data2 = data2(SpanYTrans,SpanXTrans,:);

    % reshape and transpose each image
    for y = 1:expData.numFreqs
        data1(:,:,y) = transpose(reshape(dataStack1(y, :), [expData.imgNumCols, expData.imgNumRows] ));
        data2(:,:,y) = transpose(reshape(dataStack2(y, :), [expData.imgNumCols, expData.imgNumRows] ));
    end
    
    % binning
    disp(['<>   binning data: binSize = ' num2str(binSize)])
    sizeXY = size(BinImage(data1(:,:,1),binSize));
    binData1 = zeros(sizeXY(1),sizeXY(2),length(Freqs1));
    binData2 = zeros(sizeXY(1),sizeXY(2),length(Freqs2)); 

    for y = 1:length(Freqs1)
        binData1(:,:,y) = BinImage(data1(:,:,y),binSize);
        binData2(:,:,y) = BinImage(data2(:,:,y),binSize);
    end
    
    sizeX = size(binData1,2); % binned image x-dimensions
    sizeY = size(binData1,1); % binned image y-dimensions
        
    % Correct for severely non-unity baseline by dividing pixelwise by
    % average of all frequency points
    
    binDataNorm1 = zeros(size(binData1));
    binDataNorm2 = zeros(size(binData2));
    NormalizationFactor1 = mean(binData1,3);    % compute average
    NormalizationFactor2 = mean(binData2,3);    % compute average
    
    for y = 1:length(Freqs1)
        binDataNorm1(:,:,y) = binData1(:,:,y) ./ NormalizationFactor1;
        binDataNorm2(:,:,y) = binData2(:,:,y) ./ NormalizationFactor2;
    end
    
    %% 1. GAUSSIAN BLUR
    % default = 0
    % gaussian filter on the guesses, can lead to some problems if high gradient
    if gaussianFilter ~= 0
        disp(['<>   smoothing data using gaussian blur: ' num2str(gaussianFilter)])
        gFilter = fspecial('gaussian',[20,20], gaussianFilter);
        binDataNorm1 = imfilter(binDataNorm1, gFilter, 'symmetric', 'conv');
        binDataNorm2 = imfilter(binDataNorm2, gFilter, 'symmetric', 'conv');         
    end

    disp('<>   Starting Parameter Estimation');

    %% global spectra subtraction
    binDataNorm1 = correct_global(binDataNorm1, globalFraction);
    binDataNorm2 = correct_global(binDataNorm2, globalFraction);
    
    %% first determine global guess
    meanData1 = squeeze(mean(mean(binDataNorm1,1),2));
    meanData2 = squeeze(mean(mean(binDataNorm2,1),2));
    
    %% GUESS INITIAL FIT PARAMETERS
    % start with global guess for all pixels
    guessRes1 = global_guess(binDataNorm1, Freqs1); % initial guess for GPUfit
    guessRes2 = global_guess(binDataNorm2, Freqs2); % initial guess for GPUfit

    %% local guess -> guess parameter for each pixel
    if type == 1
        disp('<>     Local Parameters')
        
        sizeX = size(binData1,2); % binned image x-dimensions
        sizeY = size(binData1,1); % binned image y-dimensions
                
        % iterate over all pixels with index (x,y)
        for x = 1:sizeX
            for y = 1:sizeY 
                pixel1 = squeeze(binDataNorm1(y,x,:)); 
                pixel2 = squeeze(binDataNorm2(y,x,:));
                
                %% Nuclear Spin polarization
                % default = 0
                if nucSpinPolFlg
                    guessRes1(y,x,:) = guessNucSpinPol(pixel1, Freqs1);    
                    guessRes2(y,x,:) = guessNucSpinPol(pixel2, Freqs2);

                %% no Nuclear Spin Polarization    
                else
                    % try getting peak positions from smoothed data
                    % LEFT
                    % NOTE: guess_peaks does not return the index anymore
                    % but Hz.
                    [pkVal1, pkRes1, fitFlg1] = guess_peaks(pixel1, meanData1, Freqs1, ...
                                                               'smoothDegree', smoothDegree, ...
                                                               'forceGuess', kwargs.forceGuess,...
                                                               'gaussianFit', kwargs.gaussianFit, ...
                                                               'pixel', [y x 1]);

                    % check if find peaks returned 3 peaks
                    % add them to badpixels
                    if ~strcmp(fitFlg1,'local')
                        badPixels1{end+1} =  {y x 2 fileNum, fitFlg1};
                    end
                    
                    if strcmp(fitFlg1,'gaussian')
                        % if it returns 3 -> replace the global guess with
                        % the local
                        resonance  = (pkRes1(1)+pkRes1(2)+pkRes1(3))/3;
                        width = 0.0005;
                        contrast = (mean(pixel1(1:10)) + pkVal1-1)';
                        baseline = mean(pixel1(1:10))-1;
                        guessRes1(y,x,:) = [resonance width contrast baseline];
                    end
                    
                    % RIGHT 
                    [pkVal2, pkLoc2, fitFlg2] = guess_peaks(pixel2, meanData2, Freqs2, ...
                                                               'smoothDegree', smoothDegree, ...
                                                               'forceGuess', kwargs.forceGuess,...
                                                               'gaussianFit', kwargs.gaussianFit, ...
                                                               'pixel', [y x 2]);
                    
                    % check if find peaks returned 3 peaks -> fitFlg is not
                    % 'local'
                    if ~strcmp(fitFlg2,'local')
                        badPixels2{end+1} =  {y x 2 fileNum, fitFlg2};
                    end
                    
                    if strcmp(fitFlg2,'gaussian')
                        resonance  = (pkLoc2(1)+pkLoc2(2)+pkLoc2(3))/3;
                        width = 0.0005;
                        contrast = (mean(pixel2(1:10)) +pkVal2-1)';
                        baseline = mean(pixel2(1:10))-1;
                        guessRes2(y,x,:) = [resonance width contrast baseline];
                    end
                end
            end
        end
    end
    
    %% Manual Guess: % manual guess redundant with forceGuess?
    if type == 3
        disp('<>    MANUAL Parameters');
        
        %Resonance 1:
        guess1 = zeros(sizeY, sizeX, 6);
        guess1(:,:,1) = 2.844728333333334;  % center location
        guess1(:,:,2) = 0.0005;             % width
        guess1(:,:,3) = 0.004083761711610;  % contrast 1
        guess1(:,:,4) = 0.004242795948007;  % contrast 2
        guess1(:,:,5) = 0.004431020178478;  % contrast 3
        guess1(:,:,6) = 0.001671817852821;  % baseline

        %Resonance 2:
        guess1 = zeros(sizeY, sizeX, 6);
        guess2(:,:,1) = 2.895892000000000;
        guess2(:,:,2) = 0.000500000000000;
        guess2(:,:,3) = 0.004529320655875;
        guess2(:,:,4) = 0.004631111505702;
        guess2(:,:,5) = 0.004663970581423;
        guess2(:,:,6) = 0.001912981765371;
    end
    
    res1 = mean2(guessRes1(:,:,1)); guessRes1(:,:,1) = guessRes1(:,:,1)-res1;
    res2 = mean2(guessRes2(:,:,1)); guessRes2(:,:,1) = guessRes2(:,:,1)-res2;
    
    % set-up guess parameters for gpu fit:
    imgpts = sizeX*sizeY;
    guess1 = reshape(guessRes1,[imgpts,6]);
    guess1 = transpose(guess1);
    guess1(1,:) = guess1(1,:)+res1; %updates to make p0 in gpu fit function

    %
    guess2 = reshape(guessRes2,[imgpts,6]);
    guess2 = transpose(guess2);
    guess2(1,:) = guess2(1,:)+res2; %updates to make p0 in gpu fit function
    
    nPixel = size(binDataNorm1,1) * size(binDataNorm1,2);
    
    fprintf('<>   INFO: initial parameter estimation complete in: %.1f s\n', toc(start)');
    
    if kwargs.gaussianFit
        disp(['<>     INFO: ' num2str(size(badPixels,1)) '/' num2str(nPixel*2) ' pixels was guessed from gaussian fit']);
    else
        disp(['<>     INFO: ' num2str(size(badPixels,1)) '/' num2str(nPixel*2) ' pixels had to be guessed from global']);
    end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RT GPU fit
    
    % data matricies:
    userinfo1 = Freqs1';
    userinfo2 = Freqs2';
    sweeplength = size(dataStack1,1);
    
    %binDataNorm1 used in fits: (y,x,freq)
    %gpu fit needs y-data and x-data seperatly...so need to make a matrix that
    %is 51 x pixels_total

    imgpts = sizeX*sizeY; % number of (x,y) pixels
    gpudata1 = reshape(binDataNorm1,[imgpts,sweeplength]); % make it into 2d matrix
    gpudata1 = transpose(gpudata1); %transpose to make it 51 x pixels
    %
    gpudata2 = reshape(binDataNorm2,[imgpts,sweeplength]); % make it into 2d matrix
    gpudata2 = transpose(gpudata2); %transpose to make it 51 x pixels
    
    
    %set-up inputs for gpu fits:
    data1 = single(gpudata1);
    userinfo1 = single(userinfo1);
    initial_parameters1 = single(guess1); %single(param1);
    %
    data2 = single(gpudata2);
    userinfo2 = single(userinfo2);
    initial_parameters2 = single(guess2); %single(param2);
    %
    tolerance = 1e-30;
    max_n_iterations = 100;
    model_id = ModelID.ESR3RT;
    %
    disp('<>   Starting GPU fits');
    % run Gpufit - Res 1
    [parameters1, states1, chi_squares1, n_iterations1, time1] = gpufit(data1, [], ...
        model_id, initial_parameters1, tolerance, max_n_iterations, [], EstimatorID.LSE, userinfo1);
    
    % run Gpufit - Res 2
    [parameters2, states2, chi_squares2, n_iterations2, time2] = gpufit(data2, [], ...
        model_id, initial_parameters2, tolerance, max_n_iterations, [], EstimatorID.LSE, userinfo2);
    disp('GPU fits Complete');
    
    %output is parameters and size is 10 x num_pixels
    %% Convergence check for each pixel
    % first consider states:
    convec1 = find(states1); % finds non-zero elements.
    con1 = ~isempty(convec1); % returns 1 if there are non-zero elements.
    if con1
        disp('GPU fits Res 1 - 2nd');
        max_n_iterations = 1000;
        upguess1 = parameters1;
        contot = size(convec1);
        for cnum = 1:contot(2)
            ncon = convec1(cnum);
            lown=max([ncon-10 1]);
            highn=min([ncon+10 size(parameters2,2)]);
            upguess2(1,ncon) = mean(parameters2(1,lown:highn));
            upguess2(2,ncon) = mean(parameters2(2,lown:highn));
            upguess2(3,ncon) = mean(parameters2(3,lown:highn));
            upguess2(4,ncon) = mean(parameters2(4,lown:highn));
            upguess2(5,ncon) = mean(parameters2(5,lown:highn));
            upguess2(6,ncon) = mean(parameters2(6,lown:highn));
        end
        initial_parameters1 = upguess1;
        [parameters1, states1, chi_squares1, n_iterations1, time1] = gpufit(data1, [], ...
            model_id, initial_parameters1, tolerance, max_n_iterations, [], EstimatorID.LSE, userinfo1);
        
    end
    
    convec2 = find(states2); % finds non-zero elements.
    con2 = ~isempty(convec2); % returns 1 if there are non-zero elements.
    if con2
        disp('GPU fits Res 2 - 2nd');
        max_n_iterations = 1000;
        upguess2 = parameters2;
        contot = size(convec2);
        for cnum = 1:contot(2)
            ncon = convec2(cnum);
            lown=max([ncon-10 1]);
            highn=min([ncon+10 size(parameters2,2)]);
            upguess2(1,ncon) = mean(parameters2(1,lown:highn));
            upguess2(2,ncon) = mean(parameters2(2,lown:highn));
            upguess2(3,ncon) = mean(parameters2(3,lown:highn));
            upguess2(4,ncon) = mean(parameters2(4,lown:highn));
            upguess2(5,ncon) = mean(parameters2(5,lown:highn));
            upguess2(6,ncon) = mean(parameters2(6,lown:highn));
        end
        initial_parameters2 = upguess2;
        [parameters2, states2, chi_squares2, n_iterations2, time2] = gpufit(data2, [], ...
            model_id, initial_parameters2, tolerance, max_n_iterations, [], EstimatorID.LSE, userinfo2);
        
        
    end
    
    
    disp('Convergence checks done');
    
    
    %% Output parameters into the proper matricies:
    
    %make parameters matrix into 3d matrix with x pixels, y pixels, and parameters
    output1 = zeros(6,sizeY,sizeX);
    output1(1,:,:) = reshape(parameters1(1,:),[sizeY,sizeX]);
    output1(2,:,:) = reshape(parameters1(2,:),[sizeY,sizeX]);
    output1(3,:,:) = reshape(parameters1(3,:),[sizeY,sizeX]);
    output1(4,:,:) = reshape(parameters1(4,:),[sizeY,sizeX]);
    output1(5,:,:) = reshape(parameters1(5,:),[sizeY,sizeX]);
    output1(6,:,:) = reshape(parameters1(6,:),[sizeY,sizeX]);
    
    % matricies with 2 dimensions for x and y pixels:
    Resonance1 = zeros(sizeY,sizeX);
    Resonance1(:,:) = output1(1,:,:);
    Width1 = zeros(sizeY,sizeX);
    Width1(:,:) = output1(2,:,:);
    ContrastA1 = zeros(sizeY,sizeX);
    ContrastA1(:,:) = output1(3,:,:);
    ContrastB1 = zeros(sizeY,sizeX);
    ContrastB1(:,:) = output1(4,:,:);
    ContrastC1 = zeros(sizeY,sizeX);
    ContrastC1(:,:) = output1(5,:,:);
    Baseline1 = zeros(sizeY,sizeX);
    Baseline1(:,:) = output1(6,:,:);
    Baseline1 = Baseline1+1;
    %
    Resonance1 = double(Resonance1);
    Width1 = double(Width1);
    ContrastA1 = double(ContrastA1);
    ContrastB1 = double(ContrastB1);
    ContrastC1 = double(ContrastC1);
    Baseline1 = double(Baseline1);
    %same for resonance 2:
    output2 = zeros(6,sizeY,sizeX);
    output2(1,:,:) = reshape(parameters2(1,:),[sizeY,sizeX]);
    output2(2,:,:) = reshape(parameters2(2,:),[sizeY,sizeX]);
    output2(3,:,:) = reshape(parameters2(3,:),[sizeY,sizeX]);
    output2(4,:,:) = reshape(parameters2(4,:),[sizeY,sizeX]);
    output2(5,:,:) = reshape(parameters2(5,:),[sizeY,sizeX]);
    output2(6,:,:) = reshape(parameters2(6,:),[sizeY,sizeX]);
    
    Resonance2 = zeros(sizeY,sizeX);
    Resonance2(:,:) = output2(1,:,:);
    Width2 = zeros(sizeY,sizeX);
    Width2(:,:) = output2(2,:,:);
    ContrastA2 = zeros(sizeY,sizeX);
    ContrastA2(:,:) = output2(3,:,:);
    ContrastB2 = zeros(sizeY,sizeX);
    ContrastB2(:,:) = output2(4,:,:);
    ContrastC2 = zeros(sizeY,sizeX);
    ContrastC2(:,:) = output2(5,:,:);
    Baseline2 = zeros(sizeY,sizeX);
    Baseline2(:,:) = 1+output2(6,:,:);
    %
    Resonance2 = double(Resonance2);
    Width2 = double(Width2);
    ContrastA2 = double(ContrastA2);
    ContrastB2 = double(ContrastB2);
    ContrastC2 = double(ContrastC2);
    Baseline2 = double(Baseline2);
    
    FitCvg = ones(sizeY,sizeX); %just to remove errors. This matrix is useless as is now.
    
    %% gpu model function
    
    modelgpu = @(p,x) 1+p(6)...
        -p(3)*p(2).^2./((x-p(1)+Ahyp).^2+p(2).^2)...
        -p(4)*p(2).^2./((x-p(1)).^2+p(2).^2)...
        -p(5)*p(2).^2./((x-p(1)-Ahyp).^2+p(2).^2);
    
    %% DEFINE FIT FUNCTION  (Old non-gpu fit), leaving here in case needed for comparision
    %
    % % % Define model function
    % % % Define without baseline offset parameter
    % % % Define with variable heights for 3 ESR peaks
    % % % Define with variable broadening for 3 ESR peaks, equal integrated area
    %
    % modelFun1 =  @(p,x) 1+p(6)...
    %     -p(3).*p(2).^2./((x-res1-p(1)+Ahyp).^2+p(2).^2)...
    %     -p(4).*p(2).^2./((x-res1-p(1)     ).^2+p(2).^2)...
    %     -p(5).*p(2).^2./((x-res1-p(1)-Ahyp).^2+p(2).^2);%For N-14 Lorentzian with variable peak amplitudes and variable baseline
    
    % modelFun2 =  @(p,x) 1+p(6)...
    %     -p(3).*p(2).^2./((x-res2-p(1)+Ahyp).^2+p(2).^2)...
    %     -p(4).*p(2).^2./((x-res2-p(1)     ).^2+p(2).^2)...
    %     -p(5).*p(2).^2./((x-res2-p(1)-Ahyp).^2+p(2).^2); %For N-14 Lorentzian with variable peak amplitudes and variable baseline
    %
    
    %% TAKE THE DIFFERENCE OF THE RESONANCES:
    ResDiff = (Resonance2 - Resonance1)/2;
    ResSum = (Resonance2 + Resonance1)/2;
    
    %HIGH PASS FOURIER-DOMAIN FILTERING
    %pixSize = 6 / Mag / BinSize; %in um
    %kCutOff = 30; % in um
    %xMax = size(ResDiff,2); yMax = size(ResDiff,1);
    %[X Y] = meshgrid(1:xMax, 1:yMax);
    %ResDiffFT = fftshift(fft2(ResDiff));
    %CircFiltHP = ((X-xMax/2).^2 / (xMax*pixSize/kCutOff)^2 +(Y-yMax/2).^2 / (yMax*pixSize/kCutOff)^2 > 1);
    %ResDiffFTFilt = CircFiltHP.*ResDiffFT;
    %ResDiffFilt = ifft2(ifftshift(ResDiffFTFilt));
    %dB = real(ResDiffFilt / gamma);
    
    if quadBGsubFlg
        %PARABOLIC BACKGROUND SUBTRACTION
        dB = QuadBGsub(ResDiff)/gamma;
        
        %PARABOLIC BG SUBTRACTION ON FIRST RESONANCE ONLY
        %   [xData, yData, zData] = prepareSurfaceData( X, Y, Resonance1 );
        %   [xData, yData, wData] = prepareSurfaceData( X, Y, FitCvg );
        %   [fitout gof] = fit([xData, yData],zData, 'poly22', 'Weight', wData);
        %   cvals = coeffvalues(fitout);
        %   fitFunction = cvals(1) + cvals(2)*X + cvals(3)*Y +...
        %     cvals(4)*X.*X + cvals(5)*X.*Y + cvals(6)*Y.*Y;
        %
        % %PARABOLIC BG SUBTRACTION ON SECOND RESONANCE ONLY
        %   if dualFitFlg
        %   [xData, yData, zData] = prepareSurfaceData( X, Y, Resonance2 );
        %   [xData, yData, wData] = prepareSurfaceData( X, Y, FitCvg );
        %   [fitout gof] = fit([xData, yData],zData, 'poly22', 'Weight', wData);
        %   cvals = coeffvalues(fitout);
        %   fitFunction = cvals(1) + cvals(2)*X + cvals(3)*Y +...
        %       cvals(4)*X.*X + cvals(5)*X.*Y + cvals(6)*Y.*Y;
        %   end
        
    else
        dB = ResDiff/gamma;
    end
    
    %EXTRA SMOOTHING AT THE END...
    % gFilterSize =4;
    % gFilter = fspecial('gaussian',[100,100],gFilterSize);
    % dB4 = imfilter(dB, gFilter, 'symmetric', 'conv');
    
    %% PLOTTING
    
    % f1 = figure; imagesc((Resonance1-zfs)/gamma);colorbar; colormap jet; axis equal tight; title 'Magnetic Field Resonance1 [Gauss]';  set(gca,'YDir','normal');
    % f2 = figure; imagesc((Resonance2-zfs)/gamma);colorbar; colormap jet; axis equal tight; title 'Magnetic Field Resonance2 [Gauss]'; set(gca,'YDir','normal');
    % f3 = figure; imagesc(dB); colorbar; colormap jet; axis equal tight; title(['Magnetic Field (Diff) Minus Gradient [Gauss].', 10, 'Std2(dB) = ' num2str(std2(dB))]); caxis([-.1 .1]);  set(gca,'YDir','normal');
    
%     if LEDimgFlg==1 f5 = figure; 
%         imagesc(transImgBin); 
%         axis equal tight; 
%         colormap gray; colorbar; 
%         title('LED image'); set(gca,'YDir','normal');    
%     end
    
    % f6 = figure; imagesc(ResDiff/gamma);colorbar; axis equal tight; title 'Magnetic Field (Diff) with Background [Gauss]';  set(gca,'YDir','normal');
%     
%     if diagplots
%         f7 = figure('units','normalized','outerposition',[0 0 1 1]);  set(gca,'YDir','normal'); suptitle('All freqs');
%         subplot(2,2,1); imagesc((Resonance1-zfs)/gamma); colorbar; colormap jet; axis equal tight; title 'Magnetic Field Resonance1 [gauss]';   set(gca,'YDir','normal');
%         subplot(2,2,2); imagesc((Resonance2-zfs)/gamma); colorbar; colormap jet; axis equal tight; title 'Magnetic Field Resonance2 [gauss]';  set(gca,'YDir','normal');
%         subplot(2,2,3); imagesc(ResDiff/gamma); colorbar; colormap jet; axis equal tight; title 'Res2 - Res1 [gauss]';  set(gca,'YDir','normal');
%         subplot(2,2,4); imagesc(dB); colorbar; colormap jet; axis equal tight; title(['Res2 - Res1 Minus Gradient [gauss]', 10, 'Std2(dB) = ' num2str(std2(dB))]); caxis([-.05 .05]);  set(gca,'YDir','normal');
%     end
    
    %strain figure
    %figure;imagesc(ResSum);colorbar; colormap jet; axis equal tight; title 'Strain [Ghz]';  set(gca,'YDir','normal');
    %
    strains = (Resonance1+Resonance2)/2-zfs;
    strains*1000;
    strains = strains - mean2(strains);
    % figure(); histogram(strains,1000000);
    % figure();imagesc(strains);colorbar; colormap jet; axis equal tight; title 'Res2 - Res1 [gauss]';  set(gca,'YDir','normal');
    
    
    %% SAVE FIT RESULTS%
    if kwargs.save
        fprintf('<>     INFO: saving data of %s\n',dataFile);
        save(fullfile(dataFolder, [dataFile, 'deltaBFit.mat']), 'dB', ...
            'Resonance1', 'Width1', 'ContrastA1', 'ContrastB1', 'ContrastC1', 'Baseline1', ...
            'Resonance2', 'Width2', 'ContrastA2', 'ContrastB2', 'ContrastC2', 'Baseline2', ...
            'binDataNorm1', 'Freqs1', 'chi_squares1', 'badPixels1',...
            'binDataNorm2', 'Freqs2', 'chi_squares2', 'badPixels2',...
            'FitCvg');
        
        if LEDimgFlg==1
            saveas(f5, [LEDimgFile 'CROPBIN.png'],'png');
            imwrite(transImgBinUint,gray, [LEDimgFile 'CROPBINPure.png'], 'png');
        end        
    end
    
    %% USEFUL TEST CODE TO PLOT FIT AT A PARTICULAR POINT I,J
    %(FOR LOWER RESONANCE):
    testplotflg = 0;
    if testplotflg
        %green is guess, red is real, cyan is smooth
        y = 440; x = 440;
        C = (squeeze(binDataNorm1(y,x,:)))';
        %coefEsts = nlinfit(Freq1, C, modelFun1, squeeze(guessRes1(y,x,:)));
        %coefEsts = [Resonance1(y,x)-res1,Width1(y,x),ContrastA1(y,x),ContrastB1(y,x),ContrastC1(y,x),Baseline1(y,x)-1] ;
        coefEsts = [Resonance1(y,x),Width1(y,x),ContrastA1(y,x),ContrastB1(y,x),ContrastC1(y,x),Baseline1(y,x)-1];
        %coefEsts = [2.8290401,0.0005000,0.0042796,0.0043262,0.0044002,0.0106340];
        %coefEsts = [2.838384628295898, 4.975722404196858e-04,0.218016400933266,1.115630149841309,  1.115103602409363,1.101928710937500];
        figure; plot(Freqs1, modelgpu(coefEsts,Freqs1),'r'); title(['x=' num2str(x) ' y=' num2str(y)]);
        hold on; plot(Freqs1, C, '.-b');plot(Freqs1, modelFun1(squeeze(guessRes1(y,x,:)),Freqs1),':g'); plot(Freqs1,smooth(C),'-c');
        %plot(Fres1, modelgpu(rtcoef,Fres1),'-k');
        hold off;
        
    
        y = 440; x = 440;
        C = (squeeze(binDataNorm2(y,x,:)))';
        % coefEsts = nlinfit(Freq2, C, modelFun2, squeeze(guessRes2(y,x,:)));
        coefEsts = [Resonance2(y,x)-res2,Width2(y,x),ContrastA2(y,x),ContrastB2(y,x),ContrastC2(y,x),Baseline2(y,x)-1] ;
        figure; plot(Freqs2, modelFun2(coefEsts,Freqs2),'r'); title(['x=' num2str(x) ' y=' num2str(y)]);
        %hold on; plot(Fres2, C, '.-b'); plot(Fres2, modelFun2(squeeze(guessRes2(y,x,:)),Fres2),':g'); plot(Fres2,smooth(C),'-c'); hold off;

        % Plot contrast
        figure('units','normalized','outerposition',[0 0 1 1]);  set(gca,'YDir','normal'); suptitle('contrast');
        subplot(3,2,1); imagesc(ContrastA1); colorbar; colormap jet; axis equal tight; title('ContrastA1'); caxis([0 .006]);  set(gca,'YDir','normal');
        subplot(3,2,2); imagesc(ContrastA2); colorbar; colormap jet; axis equal tight; title('ContrastA2'); caxis([0 .006]);  set(gca,'YDir','normal');
        subplot(3,2,3); imagesc(ContrastB1); colorbar; colormap jet; axis equal tight; title('ContrastB1'); caxis([0 .006]);  set(gca,'YDir','normal');
        subplot(3,2,4); imagesc(ContrastB2); colorbar; colormap jet; axis equal tight; title('ContrastB2'); caxis([0 .006]);  set(gca,'YDir','normal');
        subplot(3,2,5); imagesc(ContrastC1); colorbar; colormap jet; axis equal tight; title('ContrastC1'); caxis([0 .006]);  set(gca,'YDir','normal');
        subplot(3,2,6); imagesc(ContrastC2); colorbar; colormap jet; axis equal tight; title('ContrastC2'); caxis([0 .006]);  set(gca,'YDir','normal');
        
        %Plot Widths
        figure('units','normalized','outerposition',[0 0 1 1]);  set(gca,'YDir','normal'); suptitle('width');
        subplot(1,2,1); imagesc(Width1); colorbar; colormap jet; axis equal tight; title('Width1'); caxis([0 6.0E-4]);  set(gca,'YDir','normal');
        subplot(1,2,2); imagesc(Width2); colorbar; colormap jet; axis equal tight; title('Width2'); caxis([0 6.0E-4]);  set(gca,'YDir','normal');
    end
    
    
    if diagplots
        figure()
        plot(states1(1,:),'ok')
        figure()
        plot(states2(1,:),'ok')
        figure()
        plot(n_iterations1(1,:),'ok')
        figure()
        plot(n_iterations2(1,:),'ok')
        figure()
        plot(chi_squares1(1,:),'ok')
        figure()
        plot(chi_squares2(1,:),'ok')
    end
    
    
    
end