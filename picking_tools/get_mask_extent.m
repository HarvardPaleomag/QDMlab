function [x, y, w, h] = get_mask_extent(mask)
%[x, y, w, h] = get_mask_extent(mask)
% Takes a mask and retuns the extent of the smallest possible rectangle 
% around it, where values ~= 0.

[row_index, col_index, ~] = find(~isnan(mask) & mask~= 0);

% get min/max indices of mask
x = int16(min(col_index));
w = int16(max(col_index)) - x;
y = int16(min(row_index));
h = int16(max(row_index)) - y;