function fits = GPU_fit(dataFolder, binSize, kwargs)

% Parameters
% ----------
%     dataFolder: char
%         location of the data folder
%     binSize: int
%         number of pixels to be binned into one.
%         Uses 'BinImage' function

%     optional
%     ========
%     fielpolarity: int 
%         default: 0 
%         0: both polarities
%         1: first polarity only
%         2: second polarity only
%         4: Neg Pos Pos Neg
%     type: int
%         default: 0
%         0: ONLY global guess parameters
%         1: local guess parameters
%            use 'gaussianFit' to get peak positions if findpeaks fails
%         2: local guess using a gaussian pre fit
%         3: manual guess parameters %redundant???
%             todo:(set values in guess1 and guess2 matricies)
%     forceGuess: bool (false)
%     gaussianFit: bool (false)
%         Only used if the findpeaks function fails to find 3 peaks.
%         if true: uses a gaussian fit to estimate the center peak of the
%                  triplet. 
%         if false: uses the locations from the global estimation
%     gaussianFilter: numeric (0)
%         if 0: no filter is applied
%         if != 0: applies gaussian filter with a standard deviation of
%                 'gaussianFilter'. Previous versions of the code used 0.5.
%     checkPlot: bool (false)
%         display the fitted resonances. useful to establish the initial guesses for diagnostics
%     smoothDegree: int(2)
%     globalFraction: double (0.5)
%         amount of global illumination signal to subtract from data. 
%         Calls 'correct_global' function
%     save: bool (true)
%         if true the results are saved to 'dataFolder'

%     nucSpinPol: bool (false)
%         this is used for nuclear spin polarization -> NMR. Calls the
%         function guessNucSpinPol. Uses code in original state (< Nov 2020).
%     diamond: str (N14)
%         The type of diamond. Choses the type of fitting.
%
% Notes
% -----
% dataStack:     2 dimentional wwith xy pixels, frequencies
% data:          3 dimentional with y pixels, x pixels, frequencies
% binData:       3 dimentional with binned y pixels, binned x pixels, frequencies
% binDataNorm:   3 dimentional with normalized binned x pixels, normalized binned, y pixels, frequencies
%    if gaussianFilter: 3 dimentional with gaussian blurred normalized binned y pixels, gaussian blurred normalized binned, x pixels, frequencies
    
arguments
    dataFolder char
    binSize double
    kwargs.fieldPolarity (1,1) {mustBeMember(kwargs.fieldPolarity,[0,1,2, 4])} = 0
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 2
    kwargs.globalFraction (1,1) {mustBeNumeric} = 0.5
    kwargs.forceGuess (1,1) {mustBeMember(kwargs.forceGuess, [1, 0])} = 0
    kwargs.checkPlot (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0
    kwargs.gaussianFit (1,1) {mustBeMember(kwargs.gaussianFit, [1, 0])} = 0
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.nucSpinPol (1,1) {mustBeMember(kwargs.nucSpinPol, [1, 0])} = 0
    kwargs.save (1,1) {mustBeMember(kwargs.save, [1, 0])} = 1
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14'])} = 'N14'

end

tStart = tic;
fieldPolarity = kwargs.fieldPolarity;
gaussianFit = kwargs.gaussianFit;
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
polarities = {'Neg','Pos'};
sides = {'left' 'right'};
fits = struct();

%% GUESS PARAMETER ESTIMATION
for fileNum=startN:1:endN
    start = tic; % for timing 
    pol = polarities{fileNum};
    
    %%% select header and data file
    dataFile = dataFiles(fileNum).name;

    %%%
    LEDimgFile = 'laser.csv';
    
    loadStart = tic;
    fprintf('<>   loading data file:  %s\n', fullfile(dataFolder, dataFile));
    expData = load(fullfile(dataFolder, dataFile));

    fprintf('<>      loading of file %i/%i complete (%.1f s)\n', fileNum, size(startN:1:endN, 2), toc(loadStart));

    SpanXTrans = 1:expData.imgNumCols;
    SpanYTrans = 1:expData.imgNumRows;

    if LEDimgFlg==1
        transImg = load(fullfile(dataFolder,LEDimgFile));
        transImg = transImg(SpanYTrans,SpanXTrans);
        transImgBin = BinImage(transImg, binSize);
        transImgBinUint = im2uint8forExportDG(transImgBin,min(min(transImgBin)), max(max(transImgBin)) );
    end

    badPixels = struct();
    for nRes = 1:2
        side = sides{nRes};
            
        [Resfit, guess, badPixel] = fit_resonance(expData, binSize, nRes, ...
            'type',kwargs.type, 'globalFraction', kwargs.globalFraction, ...
            'diamond', kwargs.diamond,...
            'gaussianFit',gaussianFit, 'gaussianFilter', kwargs.gaussianFilter,...
            'smoothDegree', kwargs.smoothDegree, 'nucSpinPol', kwargs.nucSpinPol,...
            'checkPlot', kwargs.checkPlot);
        Resfit.fileName = fullfile(dataFolder, dataFile);
        fits.([side pol]) = Resfit;
        badPixels.([side pol]) = badPixel;
    end

    Resonance1 = fits.(['left' pol]).resonance; 
    Width1 = fits.(['left' pol]).width; 
    ContrastA1 = fits.(['left' pol]).contrastA; 
    ContrastB1 = fits.(['left' pol]).contrastB;
    
    if any(strcmp(fieldnames( fits.(['left' pol])), 'contrastC'))
        ContrastC1 = fits.(['left' pol]).contrastC;
    end
    
    Baseline1 = fits.(['left' pol]).baseline; 
    Freqs1 = fits.(['left' pol]).freq; 
    chiSquares1 = fits.(['left' pol]).chiSquares;
    p1 = fits.(['left' pol]).p;
    freq1 = fits.(['left' pol]).freq;
    
    Resonance2 = fits.(['right' pol]).resonance; 
    Width2 = fits.(['right' pol]).width; 
    ContrastA2 = fits.(['right' pol]).contrastA; 
    ContrastB2 = fits.(['right' pol]).contrastB;

    if any(strcmp(fieldnames( fits.(['right' pol])), 'contrastC'))
        ContrastC2 = fits.(['right' pol]).contrastC;
    end
    
    Baseline2 = fits.(['right' pol]).baseline;
    Freqs2 = fits.(['right' pol]).freq; 
    chiSquares2 = fits.(['right' pol]).chiSquares; 
    p2 = fits.(['right' pol]).p;
    freq2 = fits.(['right' pol]).freq;

    
    %% TAKE THE DIFFERENCE OF THE RESONANCES:
    ResDiff = (Resonance2 - Resonance1)/2;
    ResSum = (Resonance2 + Resonance1)/2;
      
    if quadBGsubFlg
        %PARABOLIC BACKGROUND SUBTRACTION
        dB = QuadBGsub(ResDiff)/gamma;     
    else
        dB = ResDiff/gamma;
    end
    
    %% SAVE FIT RESULTS%
    sizeX = size(Resonance1,2); sizeY = size(Resonance1,1); %Image dimensions
    FitCvg = ones(sizeY,sizeX); %just to remove errors. This matrix is useless as is now.
    if kwargs.save
        fprintf('<>      INFO: saving data of %s\n',dataFile);
        if strcmp(kwargs.diamond, 'N14')
            save(fullfile(dataFolder, [dataFile, 'deltaBFit.mat']), 'dB', ...
                'Resonance1', 'Width1', 'ContrastA1', 'ContrastB1', 'ContrastC1', 'Baseline1', ...
                'Freqs1', 'chiSquares1', 'p1','freq1',...
                'Resonance2', 'Width2', 'ContrastA2', 'ContrastB2', 'ContrastC2', 'Baseline2', ...
                'Freqs2', 'chiSquares2', 'p2','freq2',...
                'binSize','type','gaussianFit', 'FitCvg');
        elseif strcmp(kwargs.diamond, 'N15')
            save(fullfile(dataFolder, [dataFile, 'deltaBFit.mat']), 'dB', ...
            'Resonance1', 'Width1', 'ContrastA1', 'ContrastB1', 'Baseline1', ...
            'Freqs1', 'chiSquares1', 'p1','freq1',...
            'Resonance2', 'Width2', 'ContrastA2', 'ContrastB2', 'Baseline2', ...
            'Freqs2', 'chiSquares2', 'p2','freq2',...
            'binSize','type','gaussianFit', 'FitCvg');
        end
        
        if LEDimgFlg==1
            saveas(f5, [LEDimgFile 'CROPBIN.png'],'png');
            imwrite(transImgBinUint,gray, [LEDimgFile 'CROPBINPure.png'], 'png');
        end        
    end
end
fprintf('<>   INFO: all GPU fitting tasks completed in: %.1f s\n', toc(tStart)');