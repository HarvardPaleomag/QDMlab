function fits = GPU_fit(dataFolder, binSize, kwargs)
%[fits] = GPU_fit(dataFolder, binSize; 'fieldPolarity', 'type', 'globalFraction', 'quadBgSub', 'forceGuess', 'checkPlot', 'gaussianFit', 'gaussianFilter', 'smoothDegree', 'save', '['N15',', 'slopeCorrection', 'crop', 'fcrop')

% Parameters
% ----------
%     dataFolder: char
%         location of the data folder
%     binSize: int
%         number of pixels to be binned into one.
%         Uses 'BinImage' function
%     fielpolarity: int [0]
%         0: both polarities
%         1: first polarity only
%         2: second polarity only
%         4: Neg Pos Pos Neg
%     quadBgSub: bool [true]
%         Decides if a quadratic background is to be subtracted from the
%         data
%     type: int [0]
%         0: ONLY global guess parameters
%         1: local guess parameters
%            use 'gaussianFit' to get peak positions if findpeaks fails
%         2: local guess using a gaussian pre fit
%         3: manual guess parameters %redundant???
%             todo:(set values in guess1 and guess2 matricies)
%     forceGuess: bool [false]
%     gaussianFit: bool [false]
%         Only used if the findpeaks function fails to find 3 peaks.
%         if true: uses a gaussian fit to estimate the center peak of the
%                  triplet. 
%         if false: uses the locations from the global estimation
%     gaussianFilter: numeric [0]
%         if 0: no filter is applied
%         if != 0: applies gaussian filter with a standard deviation of
%                 'gaussianFilter'. Previous versions of the code used 0.5.
%     checkPlot: bool [false]
%         display the fitted resonances. useful to establish the initial guesses for diagnostics
%     smoothDegree: int[2]
%     globalFraction: double [0.5]
%         amount of global illumination signal to subtract from data. 
%         Calls 'correct_global' function
%     save: bool [true]
%         if true the results are saved to 'dataFolder'
%     diamond: str [N14]
%         The type of diamond. Choses the type of fitting.
%    kwargs.crop (1,1) {mustBeBoolean(kwargs.crop)} = 0
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
    kwargs.fieldPolarity (1,1) {mustBeMember(kwargs.fieldPolarity,[0,1,2,4])} = 0
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 2
    kwargs.globalFraction (1,1) {mustBeNumeric} = 0.5
    kwargs.quadBgSub (1,1) {mustBeBoolean(kwargs.quadBgSub)} = 1
    kwargs.forceGuess (1,1) {mustBeBoolean(kwargs.forceGuess)} = false;
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.gaussianFit (1,1) {mustBeBoolean(kwargs.gaussianFit)} = false;
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0;
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.save (1,1) {mustBeBoolean(kwargs.save)} = 1
    kwargs.diamond {mustBeMember(kwargs.diamond, ...
        ['N15', 'N14', 'DAC', 'singlet', 'doublet', 'triplet', 'gaussian'])} = 'N14'
    kwargs.slopeCorrection = false;
    kwargs.crop = 'none'
    kwargs.fcrop (1,1) {mustBeBoolean(kwargs.fcrop)}  = false

end

tStart = tic;
fieldPolarity = kwargs.fieldPolarity;
gaussianFit = kwargs.gaussianFit;
type = kwargs.type;

%% global variables
gamma = 0.0028;  % NV gyromagnetic ratio, in GHz / gauss
%zfs = 2.870;     % NV zero-field splitting, in GHz
%Ahyp = 0.002158; % longitudinal hyperfine for 14N

%% determine how many files need to be loaded
polarityMap = containers.Map([0,1,2,4],{[1 2];[1 1];[2 2];[1 4]});
startEnd = polarityMap(fieldPolarity);
startN = startEnd(1);
endN = startEnd(2);
%%

msg = sprintf('WORKING DIR: << %s >>', dataFolder);
logMsg('info',msg,1,0);

dataFiles = dir(fullfile(dataFolder,'run_0000*.mat'));
% filter files that are not run0000.mat
idx = ~cellfun('isempty', regexpi({dataFiles.name}, 'run_[0-9]{5}\.mat$','match'));

dataFiles = dataFiles(idx);
polarities = {'Neg','Pos'};
sides = {'left' 'right'};
fits = struct();

%% GUESS PARAMETER ESTIMATION
lowCropIdx = [0,0];
highCropIdx = [0,0];

for fileNum=startN:1:endN
    pol = polarities{fileNum};
    
    %%% select header and data file
    dataFile = dataFiles(fileNum).name;
    headerFile = strrep(dataFile,'.mat','_header.txt');
    
    loadStart = tic; % for timing 
    msg = ['loading data file: ', fullfile(dataFolder, dataFile)];
    logMsg('debug',msg,1,0);
    
    % read data & header
    expData = load(fullfile(dataFolder, dataFile));
    header = read_header(fullfile(dataFolder, headerFile));
    
    % add path to header struct
    header.dFile = fullfile(dataFolder, dataFile);
    header.headerFile = fullfile(dataFolder, headerFile);
    
    msg = sprintf('loading of file %i/%i complete (%.1f s)', fileNum, size(startN:1:endN, 2), toc(loadStart));
    logMsg('info',msg,1,1);

    pixelAlerts = struct();
    
    for nRes = 1:2
        side = sides{nRes};
        fRanges = get_franges(expData, header);
        % get the fcrop values for each side only one time
        if kwargs.fcrop
            n = expData.numFreqs;
            if nRes == 1 & all(lowCropIdx == [0,0])
                lowCropIdx = pick_fcrop(expData.disp1, fRanges{1});
            elseif nRes == 2 & all(highCropIdx == [0,0])
                highCropIdx = pick_fcrop(expData.disp2, fRanges{2});
            end
        end
        
        if kwargs.fcrop & nRes == 1
            kwargs.fcrop = lowCropIdx;
        elseif kwargs.fcrop & nRes == 2
            kwargs.fcrop = highCropIdx;
        end
        
        Resfit = fit_resonance(expData, header, binSize, nRes, ...
            'type',kwargs.type, ...
            'globalFraction', kwargs.globalFraction, ...
            'diamond', kwargs.diamond,...
            'slopeCorrection', kwargs.slopeCorrection,...
            'gaussianFit',gaussianFit, ...
            'gaussianFilter', kwargs.gaussianFilter,...
            'smoothDegree', kwargs.smoothDegree, ...
            'crop', kwargs.crop, ...
            'fcrop', kwargs.fcrop, ...
            'checkPlot', kwargs.checkPlot);
        
        Resfit.fileName = fullfile(dataFolder, dataFile);
        fits.([side pol]) = Resfit;
    end

    Resonance1 = fits.(['left' pol]).resonance; 
    Width1 = fits.(['left' pol]).width; 
    ContrastA1 = fits.(['left' pol]).contrastA;
    
    if any(strcmp(fieldnames( fits.(['left' pol])), 'contrastB'))
        ContrastB1 = fits.(['left' pol]).contrastB;
    end
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
    
    if any(strcmp(fieldnames( fits.(['right' pol])), 'contrastB'))
        ContrastB2 = fits.(['right' pol]).contrastB; 
    end

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
      
    if kwargs.quadBgSub
        %PARABOLIC BACKGROUND SUBTRACTION
        dB = QuadBGsub(ResDiff)/gamma;     
    else
        dB = ResDiff/gamma;
    end
    
    % fit convergance, if the fit failed for whatever reason, the value for this pixel is 1 will be
    pixelAlerts = fits.(['left' pol]).states ~= 0 | fits.(['right' pol]).states ~= 0;
    fits.(['pixelAlerts' pol]) = pixelAlerts;
    fits.(['header' pol]) = header;
    
    %% SAVE FIT RESULTS%
    if kwargs.save
        folderName = sprintf('%ix%iBinned', binSize, binSize);
        
        if isfolder(fullfile(dataFolder, folderName))
            msg = sprintf('folder < %s > already exists', folderName);
            logMsg('warn',msg,1,0);
        else
            msg = sprintf('creating folder < %s >', folderName);
            logMsg('info',msg,1,0);
            mkdir(fullfile(dataFolder, folderName));
        end
        
        msg = sprintf('saving data of %s',dataFile);
        logMsg('info',msg,1,0);        
        
        switch kwargs.diamond
            case {'N14', 'triplet'}
            save(fullfile(dataFolder, folderName, [dataFile, 'deltaBFit.mat']), 'dB', ...
                'Resonance1', 'Width1', 'ContrastA1', 'ContrastB1', 'ContrastC1', 'Baseline1', ...
                'Freqs1', 'chiSquares1', 'p1','freq1',...
                'Resonance2', 'Width2', 'ContrastA2', 'ContrastB2', 'ContrastC2', 'Baseline2', ...
                'Freqs2', 'chiSquares2', 'p2','freq2',...
                'binSize','type','gaussianFit', 'pixelAlerts');
            
            case {'N15', 'doublet'}
                save(fullfile(dataFolder, folderName, [dataFile, 'deltaBFit.mat']), 'dB', ...
                'Resonance1', 'Width1', 'ContrastA1', 'ContrastB1', 'Baseline1', ...
                'Freqs1', 'chiSquares1', 'p1','freq1',...
                'Resonance2', 'Width2', 'ContrastA2', 'ContrastB2', 'Baseline2', ...
                'Freqs2', 'chiSquares2', 'p2','freq2',...
                'binSize','type','gaussianFit', 'pixelAlerts');
        
            case {'DAC', 'singlet', 'gaussian'}
                save(fullfile(dataFolder, folderName, [dataFile, 'deltaBFit.mat']), 'dB', ...
                'Resonance1', 'Width1', 'ContrastA1', 'Baseline1', ...
                'Freqs1', 'chiSquares1', 'p1','freq1',...
                'Resonance2', 'Width2', 'ContrastA2', 'Baseline2', ...
                'Freqs2', 'chiSquares2', 'p2','freq2',...
                'binSize','type','gaussianFit', 'pixelAlerts');
        end
   
    end
end
msg = sprintf('all GPU fitting tasks completed in: %.1f s', toc(tStart));
logMsg('FINAL',msg,1,0);
end