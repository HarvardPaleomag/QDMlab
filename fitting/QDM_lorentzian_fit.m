function fits = QDM_lorentzian_fit(nFolders, binSizes, kwargs)
% :code:`QDM_lorentzian_fit` uses GPU_fit to calculate the field values for each pixel
% and then determines B111 field values from the different polarities.
%
% Parameters
% ----------
%   nFolders:
%   binSize:
%   fieldPolarity: (0)
%     0 = neg & pos,
%     1 = neg,
%     2 = pos,
%     4 = nppn
%   globalFraction: numeric (0.5)
%     Amount of global ilumination signal to be subtracted from the
%     measurements before fitting
%   type: [0,1,2] (0)
%     use global (0) or local (1) guesses for the fit
%   smoothDegree: int (2)
%     gaussian smoothing before fit
%   gaussianFit: bool (false)
%   gaussianFilter: numeric (0)
%   nucSpinPol: bool (false)
%   save: bool (true)
%   show diagnostics plots
%   plotGuessSpectra: bool (1)
%   forceGuess: bool (0)
%

arguments
    nFolders
    binSizes double
    % keyword arguments
    kwargs.fieldPolarity (1,1) {mustBeMember(kwargs.fieldPolarity,[0,1,2,4])} = 0
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 2
    kwargs.globalFraction (1,1) {mustBeNumeric} = 0.25
    kwargs.forceGuess (1,1) {mustBeBoolean(kwargs.forceGuess)} = 0
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)} = 0
    kwargs.plotGuessSpectra (1,1) {mustBeBoolean(kwargs.plotGuessSpectra)} = 0
    kwargs.gaussianFit (1,1) {mustBeBoolean(kwargs.gaussianFit)} = 0
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.nucSpinPol (1,1) {mustBeBoolean(kwargs.nucSpinPol)} = 0
    kwargs.save (1,1) {mustBeBoolean(kwargs.save)} = 1
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14'])} = 'N14'
end

% check if there is more than one folder
nFolders = correct_cell_shape(nFolders);

% check if there is one or more binSize
if isnumeric(binSizes)
    binSizes = [binSizes];
end

% select field polarity
fp = containers.Map({0 1 2 4},{'np  ' 'n   ' 'p   ' 'nppn'});
type = fp(kwargs.fieldPolarity);

for dataFolder = nFolders
    dataFolder = dataFolder{:};
    for n=1:size(binSizes,2)
        binSize=binSizes(n);
        %   GPU_fit_QDM(INFILE,polarities,bin,neighborguess,diagnostics)
        fits = GPU_fit(dataFolder, binSize,...
                        'fieldPolarity',kwargs.fieldPolarity, ...
                        'type', kwargs.type,...
                        'globalFraction', kwargs.globalFraction,...
                        'diamond', kwargs.diamond,...
                        'gaussianFit', kwargs.gaussianFit,...
                        'gaussianFilter', kwargs.gaussianFilter,...
                        'forceGuess', kwargs.forceGuess,...
                        'checkPlot', kwargs.checkPlot,...
                        'smoothDegree', kwargs.smoothDegree,...
                        'nucSpinPol', kwargs.nucSpinPol,...
                        'save', kwargs.save);
                    
        fits.kwargs = kwargs;
        fits.nFolder = dataFolder;
        
        folderName=[num2str(binSize) 'x' num2str(binSize) 'Binned'];
        mkdir(fullfile(dataFolder, folderName));
        fits = plotResults_CommLine(dataFolder, folderName, type, fits, binSize);
        save(fullfile(dataFolder, sprintf('final_fits_(%ix%i).mat', binSize, binSize)), '-struct', 'fits');

        %foldername=[num2str(bin) 'x' num2str(bin) 'Binned_' num2str(GF)];
        movefile(fullfile(dataFolder, 'run_00000.matdeltaBFit.mat'),fullfile(dataFolder, folderName))
        movefile(fullfile(dataFolder, 'run_00001.matdeltaBFit.mat'),fullfile(dataFolder, folderName))
        if strcmp(type, 'nppn')
        movefile(fullfile(dataFolder, 'run_00002.matdeltaBFit.mat'),fullfile(dataFolder, folderName))
        movefile(fullfile(dataFolder, 'run_00003.matdeltaBFit.mat'),fullfile(dataFolder, folderName))
        end
%         movefile(fullfile(dataFolder, 'B111dataToPlot.mat'),fullfile(dataFolder, folderName))
%         movefile(fullfile(dataFolder, 'negCurrent.png'),fullfile(dataFolder, folderName))
%         movefile(fullfile(dataFolder, 'posCurrent.png'),fullfile(dataFolder, folderName))
%         movefile(fullfile(dataFolder, 'ferromagImg.png'),fullfile(dataFolder, folderName))
%         movefile(fullfile(dataFolder, 'paramagImg.png'),fullfile(dataFolder, folderName))
%         copyfile(fullfile(dataFolder, 'ledImg.png'),fullfile(dataFolder, folderName))
%         movefile(fullfile(dataFolder, 'allPlots.png'),fullfile(dataFolder, folderName))
        copyfile(fullfile(dataFolder, 'laser.jpg'),fullfile(dataFolder, folderName))
    end
end

end
