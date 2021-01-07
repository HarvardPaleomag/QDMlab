function lMask = limit_mask(mask)
%Takes a mask and retuns the smallest possible rectangle where values ~= 0.

[x, y, w, h] = get_mask_extent(mask);

% cut around data
% disp('<>   cutting around mask')
lMask = mask(y:y+h, x:x+w);