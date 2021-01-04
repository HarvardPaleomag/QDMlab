function guess = global_guess(data, freqs, varargin)
%{
Returns a global guess for the given dataset

Returns
-------
    array with 6 x X x Y elements, where X and Y are the indices of the
    pixels
%}
%% parser
inParse = inputParser;
addRequired(inParse, 'meanData');
addRequired(inParse, 'freqs');
addParameter(inParse, 'smoothDegree', 3);
addParameter(inParse, 'minPeakDistance', 0);
addParameter(inParse, 'forceGuess', 0);
addParameter(inParse, 'checkPlot', false, @islogical);

parse(inParse, data, freqs, varargin{:});

%%
disp('<>      generating initial guess from global resonance parameters');

sizeX = size(data,2); % binned image x-dimensions
sizeY = size(data,1); % binned image y-dimensions
meanData = squeeze(mean(mean(data,1),2));

%% Resonance 1:
[pkVal, pkLoc, fitFlg] = guess_peaks(meanData, meanData, freqs, ...
                               'smoothDegree', inParse.Results.smoothDegree, ...
                               'forceGuess', inParse.Results.forceGuess,...
                               'checkPlot', inParse.Results.checkPlot);

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