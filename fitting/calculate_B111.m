function calculate_B111(dataFolder, binSizes, varargin)
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
%         smoothDegree: int (2)
%             gaussian smoothing before fit
%         gaussianFit: bool (false)
%         gaussianFilter: numeric (0)
%         nucSpinPol: bool (false)
%         save: bool (true)
%         checkPlot: bool (false)
%             show diagnostics plots
%         plotGuessSpectra: bool (1)
%         forceGuess: bool (0)
inParse = inputParser;
addRequired(inParse, 'dataFolder');
addRequired(inParse, 'binSize');
addParameter(inParse, 'fieldPolarity', 0);
addParameter(inParse, 'type', 0);

addParameter(inParse, 'globalFraction', 0.5); %0.5;%0.237

% force initial guesses for FOVs with very strong fields
addParameter(inParse, 'forceGuess', 0);
% 
addParameter(inParse, 'checkPlot', false, @islogical);
addParameter(inParse, 'plotGuessSpectra', 1);

addParameter(inParse, 'gaussianFit', false, @islogical);
addParameter(inParse, 'gaussianFilter', 0, @isnumeric); %0.5;%0.237
addParameter(inParse, 'smoothDegree', 2);
addParameter(inParse, 'nucSpinPol', false, @islogical);
addParameter(inParse, 'save', true, @islogical);

parse(inParse, dataFolder, binSizes, varargin{:});

if inParse.Results.fieldPolarity == 4
  type='nppn';
else
  type='np  ';
end

for n=1:size(binSizes,2)
  binSize=binSizes(n);
%   GPU_fit_QDM(INFILE,polarities,bin,neighborguess,diagnostics)
  GPU_fit(dataFolder,binSize,...
      'fieldPolarity',inParse.Results.fieldPolarity, ...
      'type', inParse.Results.type,...
      'globalFraction', inParse.Results.globalFraction,...
      'gaussianFit', inParse.Results.gaussianFit...,
      'gaussianFilter', inParse.Results.gaussianFilter)
      
  plotResults_CommLine(dataFolder,type)
  foldername=[num2str(binSize) 'x' num2str(binSize) 'Binned'];

  %foldername=[num2str(bin) 'x' num2str(bin) 'Binned_' num2str(GF)];
  mkdir(fullfile(dataFolder, foldername));
  movefile(fullfile(dataFolder, 'run_00000.matdeltaBFit.mat'),fullfile(dataFolder, foldername))
  movefile(fullfile(dataFolder, 'run_00001.matdeltaBFit.mat'),fullfile(dataFolder, foldername))
  if type=='nppn';
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

