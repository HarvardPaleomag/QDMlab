function [idx, rIdx, cIdx] = return_bin_data(row, col, kwargs)
%[idx, rIdx, cIdx] = return_bin_data(row, col; 'binSize', 'sizeUnbinned', 'type')
% returns the data used for the bin binnedData(row,col,:)
% 
% Parameters
% ----------
%   row:
%   col:
%   binSize: (4)
%   sizeUnbinned: ([1200, 1920])
%   type: ('binDataNorm')
% 
% Returns
% ----------
%   idx:
%   rIdx:
%   cIdx:

arguments
    row
    col
    kwargs.binSize = 4;
    kwargs.sizeUnbinned = [1200, 1920];
    kwargs.type = 'binDataNorm';
end

idx = [];

rIdx = row * kwargs.binSize - 1:row * kwargs.binSize - 1 + kwargs.binSize - 1;
cIdx = col * kwargs.binSize - 1:col * kwargs.binSize - 1 + kwargs.binSize - 1;

for r = rIdx
    for c = cIdx
        i = xy2index(r, c, kwargs.sizeUnbinned, 'type', kwargs.type);
        idx = [idx, i];
    end
end
