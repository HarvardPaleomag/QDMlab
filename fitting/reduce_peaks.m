function [peakValue, peakLocation] = reduce_peaks(pkVals, pkLocs)
% REDUNDANT
%{ 
Function that selects 3 peaks from a set of peaks > 3. 
It returns the highest peaks and its closest neighbours to the left and right.

Parameters
----------
    pkVals: double
        list of peak values sorted by height ascending
    pkLocs: double
        list of peak locations sorted by height ascending
Returns
-------
    peakValue, peakLocation for exactly 3 peaks
%}

% get last element (i.e. maximum peak)
mx = pkLocs(end);
% sort indices of peaks
lst = sort(pkLocs);
% get index of the maximum peak in sorted list
idx = find(lst==mx);

% if left/righ most peak is highest take the second highest
if idx == 1 || idx == size(pkLocs,1)
    mx = pkLocs(end - 1);
    idx = find(lst==mx);
end

% in case it is again an edge case return no peaks -> global
if idx == 1 || idx == size(pkLocs,1)
    peakValue = [];
    peakLocation = [];
    return
else
    % get the left and right peak of the maximum
    new_pkLocs = lst(idx-1:idx+1);

    % get the original peak indices of the 3 peaks
    original_idx_peak = zeros([1,3]);
    for i = 1:size(new_pkLocs)
        elem = new_pkLocs(i);
        original_idx_peak(i) = find(pkLocs==elem);
    end
    original_idx_peak = sort(original_idx_peak);
end
% get the original peaks in the same order
peakValue = pkVals(original_idx_peak);
peakLocation = pkLocs(original_idx_peak);
end