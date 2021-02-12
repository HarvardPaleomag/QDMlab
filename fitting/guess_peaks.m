function [peakValue, peakLocation, fitFlg] = guess_peaks(data, globalData, freqs, kwargs)
%{
helper function to determine an initial guess of a QDM spectra (i.e.
pixel).

Parameters
----------
    required
    ========
    data: double
        data array with numFreqs entries
    globalData: double
        mean of all pixels in image. Used only for plotting
    freqs: double
        frequency range

    optional
    ========
    pixel: 1x3
        x y location and left(1) right(2) spectrum 
        used to find pixels that were not guessed correct (i.e. global)
    forceGuess:
        HAVE NOT TESTED
        used to force a guess
    smoothDegree: int
        default: 3
        how much the data ius smoothed before finding peaks
    gaussianFit: bool
        if true: will use a simple gaussian fit to determine the location
        of the center peak and use the hyp. splitting to gues the peaks
        left and right to it.
        if false: will use the global to estimate the peak locations
    checkPlot: bool
        creates a plot to check for wrong guesses

Returns
-------
    peakValue: double
        three values for the peaks
    peakLocation: double
        three peakindeces
    fitFlg: [0 1 2]
        0: local guess worked
        1: gaussian (local failed)
        2: global (local failed)
%}
arguments
    data double
    globalData double
    freqs double

    kwargs.forceGuess (1,1) {mustBeMember(kwargs.forceGuess, [1, 0])} = 0
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)} = 0
    kwargs.gaussianFit (1,1) {mustBeBoolean(kwargs.gaussianFit)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.pixel  (1,3) {mustBeNumeric} = [nan nan nan]
end
smoothDegree = kwargs.smoothDegree;

pixel = kwargs.pixel;
peaksNotFound = 0;
fitFlg = 'local'; 

if kwargs.forceGuess
    peakValue = [.0041 .0044 .0040]' ;
    peakLocation = [15 23 32]' ;
    return
end

ahyp = 0.002158;
minpeakdistance = abs(round(ahyp / mean(diff(freqs))));

dataSmoothed = smooth(1-data,'sgolay', smoothDegree);

% try
    % try getting peak positions from smoothed data
    [peakValue, peakLocation] = findpeaks(dataSmoothed,...
                                'minpeakdistance', minpeakdistance/2, ...
                                'minpeakheight',0.5*(min(dataSmoothed)+max(dataSmoothed)),... %old version pre mike
                                'MinPeakProminence',0.0003, ...
                                'NPeaks',3,...
                                'sortstr','ascend');
    % index -> GHz values
	peakLocation = freqs(peakLocation).';
    
    if kwargs.checkPlot
        figure;
        hold on
        plot(freqs, data,'b.-','DisplayName','pixel data')
        plot(freqs, 1-dataSmoothed,'b--','DisplayName','pixel data smooth')
        plot(peakLocation, 1-peakValue, 'ob','DisplayName',sprintf('pixel peaks (%i)', size(peakLocation,1)))
        plot(freqs, globalData/max(globalData)*(1-max(dataSmoothed)), 'r.:','DisplayName','global data')
        legend()
    end
% catch
%     fprintf('<>     ERROR: findpeaks function failed for pixel (y,x): (%3i, %3i) Res. #%1i -> trying gaussian Fit on data',pixel);
%     peaksNotFound = true;
% end

%% First check
% if there are less than 3 peaks repeat findpeaks routine on unsmoothed
% data
if peaksNotFound || length(peakLocation) < 3
    if kwargs.gaussianFit
        fprintf('<>      INFO: less than 3 peaks found in pixel (y,x): (%3i, %3i) Res. #%1i -> trying gaussian Fit on data\n',pixel);
        y = data;
        y = reshape(y, size(freqs));
        y = y-max(data);
        f = fit(freqs.', y.','gauss1');
        
        centerLoc = f.b1;
        peakLocation = [centerLoc-ahyp, centerLoc, centerLoc+ahyp];
        
        if any(peakLocation < 0)
            fitFlg = 2;
            fprintf('<>      WARNING: peakLocation in pixel (y,x): (%3i, %3i) Res. #%1i out of range -> using global data', pixel);
            return %returns no value -> GPU_fit will use global value
        end
        peakValue = [1; 1; 1] * max(data)*f.a1; %data(peakLocation);
        
        % make sure peaks are sorted ascending
        [peakValue, order] = sort(peakValue, 'descend');
        peakLocation = peakLocation(order);
        
        if kwargs.checkPlot
            plot(freqs, f(freqs)+max(data),'k','DisplayName','gaussian fit')
            plot(peakLocation, peakValue+max(data), 'ok','DisplayName',['gaussian (+-' num2str(minpeakdistance) ') peaks'])
        end
        fitFlg = 1;

    else
        % if numPeaks is smaller than 3
        % return 0 0 -> GPU_fit then uses the global guess
        fitFlg = 2;
        fprintf('<>      WARNING: less than 3 peaks found in pixel (y,x): (%3i, %3i) Res. #%1i out of range -> using global data\n', pixel);
        peakValue = 0;
        peakLocation = 0;
    end
    
end





