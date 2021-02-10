function fits = QDM_lorentzian_fit(dataFolders, binSizes, kwargs)
% calculate_B111 uses GPU_fit to calculate the field values for each pixel
% and then determines B111 field values from the different polarities.
% 
% Parameters
% ----------
%     positional
%     ==========
%         dataFolder:
%         binSize:
% 
%     optional
%     ========
%         fieldPolarity: (0)
%             0 = neg & pos, 
%             1 = neg, 
%             2 = pos, 
%             4 = nppn
%         globalFraction: numeric (0.5)
%             Amount of global ilumination signal to be subtracted from the 
%             measurements before fitting
%         type: [0,1,2] (0)
%             use global (0) or local (1) guesses for the fit
%         diamond: str ('N14')
%             The type of diamond. Choses the type of fitting.
%         smoothDegree: int (2)
%             gaussian smoothing before fit
%         gaussianFit: bool (false)
%         gaussianFilter: numeric (0)
%         nucSpinPol: bool (false)
%         save: bool (true)
%         : bool (false)
%             show diagnostics plots
%         plotGuessSpectra: bool (1)
%         forceGuess: bool (0)

arguments
    dataFolders
    binSizes double
    % keyword arguments
    kwargs.fieldPolarity (1,1) {mustBeMember(kwargs.fieldPolarity,[0,1,2,4])} = 0
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 2
    kwargs.globalFraction (1,1) {mustBeNumeric} = 0.5
    kwargs.forceGuess (1,1) {mustBeMember(kwargs.forceGuess, [1, 0])} = 0
    kwargs.checkPlot (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0
    kwargs.plotGuessSpectra (1,1) {mustBeMember(kwargs.plotGuessSpectra, [1, 0])} = 0
    kwargs.gaussianFit (1,1) {mustBeMember(kwargs.gaussianFit, [1, 0])} = 0
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.nucSpinPol (1,1) {mustBeMember(kwargs.nucSpinPol, [1, 0])} = 0
    kwargs.save (1,1) {mustBeMember(kwargs.save, [1, 0])} = 1
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14'])} = 'N14'
end

% check if there is more than one folder
dataFolders = correct_cell_shape(dataFolders);

% check if there is one or more binSize
if isnumeric(binSizes)
    binSizes = [binSizes];
end

% select field polarity
if kwargs.fieldPolarity == 4
  type='nppn';
else
  type='np  ';
end

for dataFolder = dataFolders
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

      plotResults_CommLine(dataFolder,type)
      foldername=[num2str(binSize) 'x' num2str(binSize) 'Binned'];

      %foldername=[num2str(bin) 'x' num2str(bin) 'Binned_' num2str(GF)];
      mkdir(fullfile(dataFolder, foldername));
      movefile(fullfile(dataFolder, 'run_00000.matdeltaBFit.mat'),fullfile(dataFolder, foldername))
      movefile(fullfile(dataFolder, 'run_00001.matdeltaBFit.mat'),fullfile(dataFolder, foldername))
      if strcmp(type, 'nppn')
        movefile(fullfile(dataFolder, 'run_00002.matdeltaBFit.mat'),fullfile(dataFolder, foldername))
        movefile(fullfile(dataFolder, 'run_00003.matdeltaBFit.mat'),fullfile(dataFolder, foldername))
      end
      movefile(fullfile(dataFolder, 'B111dataToPlot.mat'),fullfile(dataFolder, foldername))
      movefile(fullfile(dataFolder, 'negCurrent.png'),fullfile(dataFolder, foldername))
      movefile(fullfile(dataFolder, 'posCurrent.png'),fullfile(dataFolder, foldername))
      movefile(fullfile(dataFolder, 'ferromagImg.png'),fullfile(dataFolder, foldername))
      movefile(fullfile(dataFolder, 'paramagImg.png'),fullfile(dataFolder, foldername))
      movefile(fullfile(dataFolder, 'ledImg.png'),fullfile(dataFolder, foldername))
      movefile(fullfile(dataFolder, 'allPlots.png'),fullfile(dataFolder, foldername))
      copyfile(fullfile(dataFolder, 'laser.jpg'),fullfile(dataFolder, foldername))
    end
end

end

