function fits = QDM_lorentzian_fit(kwargs)
%[fits] = QDM_lorentzian_fit('nFolders', 'binSizes', 'fieldPolarity', 'type', 'globalFraction', 'forceGuess', 'checkPlot', 'plotGuessSpectra', 'gaussianFit', 'gaussianFilter', 'smoothDegree', 'nucSpinPol', 'save', 'diamond', 'slopeCorrection')
% and then determines B111 field values from the different polarities.
%
% Parameters
% ----------
%   nFolders: cell ['none']
%       if 'none' lets you pick a folder
%   binSize: cell(int) [4]
%   fieldPolarity: (0)
%     0 = neg & pos,
%     1 = neg,
%     2 = pos,
%     4 = nppn
%   globalFraction: numeric (0.5)
%     Amount of global ilumination signal to be subtracted from the
%     measurements before fitting
%   type: [0,1,2] (2)
%     0: for global guess, 
%     1: local with gaussian fit (OLD)
%     2: local with GPU fit NEW
%   smoothDegree: int (2)
%     gaussian smoothing before fit
%   gaussianFit: bool (false)
%   gaussianFilter: numeric (0)
%   nucSpinPol: bool (false)
%   save: bool (true)
%   show diagnostics plots
%   plotGuessSpectra: bool (1)
%   forceGuess: bool (0)
%   slope_correction: bool [false]
%       uses a linear slope correction on the raw data to determine the
%       initial guess. 
%       Note: only works for type = 2

arguments
    kwargs.nFolders {foldersMustExist(kwargs.nFolders)} = 'none';
    kwargs.binSizes = 'none'
    % keyword arguments
    kwargs.fieldPolarity (1,1) {mustBeMember(kwargs.fieldPolarity,[0,1,2,4])} = 0
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 2
    kwargs.globalFraction = 'none';
    kwargs.forceGuess (1,1) {mustBeBoolean(kwargs.forceGuess)} = false;
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.plotGuessSpectra (1,1) {mustBeBoolean(kwargs.plotGuessSpectra)} = false;
    kwargs.gaussianFit (1,1) {mustBeBoolean(kwargs.gaussianFit)} = false;
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.save (1,1) {mustBeBoolean(kwargs.save)} = 1
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14'])} = 'N14'
    kwargs.slopeCorrection = false;
end

msg = sprintf('please use ODMR_to_B111. QDM_lorentzian_fit will be removed in  a later update');
logMsg('deprecated',msg,1,0);

kwargs = namedargs2cell(kwargs);
fits = ODMR_to_B111(kwargs{:});