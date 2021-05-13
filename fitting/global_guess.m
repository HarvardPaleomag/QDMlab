function guess = global_guess(data, freq, kwargs)
%[guess] = global_guess(data, freq; <checkPlot>, <forceGuess>, <minPeakDistance>, <smoothDegree>)
% Returns a global guess for the given dataset
% 
% Returns
% -------
%     array with 6 x X x Y elements, where X and Y are the indices of the
%     pixels

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

sizeX = size(data,2); % binned image x-dimensions
sizeY = size(data,1); % binned image y-dimensions
% meanData = squeeze(mean(mean(data,1),2));
% m1 = mean(data,1);
% m2 = nanmean(squeeze(m1),1);
meanData = squeeze(nanmean(data,[1,2]));

%% Resonance 1:
[pkVal, pkLoc, fitFlg] = guess_peaks(meanData, meanData, freq, ...
                               'smoothDegree', kwargs.smoothDegree, ...
                               'forceGuess', kwargs.forceGuess,...
                               'checkPlot', kwargs.checkPlot);

Rguess = [(pkLoc(1)+pkLoc(2)+pkLoc(3))/3  0.0005  ( mean(meanData(1:10)) +pkVal-1)' mean(meanData(1:10))-1 ]; %resonance [GHz], Width [GHz], (contrast1; contrast2; contrast3)', baseline
Rguess = Rguess';

% write guess
guess = zeros(sizeY, sizeX, 6);
guess(:,:,1) = Rguess(1);
guess(:,:,2) = Rguess(2);
guess(:,:,3) = Rguess(3);
guess(:,:,4) = Rguess(4);
guess(:,:,5) = Rguess(5);
guess(:,:,6) = Rguess(6);

guess = single(guess);