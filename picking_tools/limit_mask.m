function lMask = limit_mask(mask, kwargs)
%Takes a mask and retuns the smallest possible rectangle where values ~= 0.
arguments
    mask
    kwargs.reference = 'none'
end

if strcmp(kwargs.reference,'none')
    kwargs.reference = mask;
end

[x, y, w, h] = get_mask_extent(kwargs.reference);

% cut around data
% disp('<>   cutting around mask')
lMask = mask(y:y+h, x:x+w);