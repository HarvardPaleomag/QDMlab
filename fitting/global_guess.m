function guess = global_guess(data, freq, kwargs)
%[guess] = global_guess(data, freq; 'forceGuess', 'checkPlot', 'smoothDegree', 'minPeakDistance')
% Returns a global guess for the given dataset
% 
% Returns
% -------
%     array with 6 x nPixel elements

arguments
    data double
    freq double
    % keyword arguments
    kwargs.forceGuess (1,1) {mustBeBoolean(kwargs.forceGuess)} = false
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)} = false
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.minPeakDistance (1,1) {mustBeNumeric} = 0
end

%%
msg = sprintf('generating initial guess from global resonance parameters');
logMsg('info',msg,1,0);

%% reshape (MxNxF) data
if ndims(data) == 3
    [~, ~, sweepLength] = size(data); % binned image x-dimensions
    data = transpose(reshape(data, [], sweepLength));
end

[~, nPixel] = size(data);

meanData = squeeze(mean(data, 2, 'omitnan'));

[pkVal, pkLoc, ~] = guess_peaks(meanData, freq, ...
                               'smoothDegree', kwargs.smoothDegree, ...
                               'forceGuess', kwargs.forceGuess,...
                               'checkPlot', kwargs.checkPlot);

Rguess = [(pkLoc(1)+pkLoc(2)+pkLoc(3))/3  0.0005  ( mean(meanData(1:10)) +pkVal-1)' mean(meanData(1:10))-1 ]; %resonance [GHz], Width [GHz], (contrast1; contrast2; contrast3)', baseline
guess = repmat(Rguess, nPixel, 1);
guess = single(guess');