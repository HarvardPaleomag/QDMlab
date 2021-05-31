function [idx, rIdx, cIdx] = return_bin_data(row, col, kwargs)
%[window] = return_bin_data(data, row, col, binSize)
% returns the data used for the bin binnedData(row,col,:)
arguments
    row
    col
    kwargs.binSize = 4;
    kwargs.nCols = 1920;
    kwargs.nRows = 1200;
    kwargs.type = 'binDataNorm';
end

idx = [];
binnedIdx = [];

rIdx = row * kwargs.binSize - 1:row * kwargs.binSize - 1 + kwargs.binSize - 1;
cIdx = col * kwargs.binSize - 1:col * kwargs.binSize - 1 + kwargs.binSize - 1;

for c = cIdx
    for r = rIdx
        i = xy2index(c, r, kwargs.nCols, 'type', kwargs.type);
        idx = [idx, i];
    end
end
