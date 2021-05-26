function lMask = limit_mask(mask, kwargs)
%[lMask] = limit_mask(mask; 'reference')
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
msg = sprintf('cutting around mask lower left (x:%i, y:%i, w:%i, h:%i)', x,y,w,h);
logMsg('debug',msg,1,0);
lMask = mask(y:y+h, x:x+w);