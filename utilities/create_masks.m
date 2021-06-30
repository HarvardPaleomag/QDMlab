function [nMasks, nROI] = create_masks(data, selectionThreshold, kwargs)
%%
arguments
    data
    selectionThreshold
    kwargs.nROI = false;
    kwargs.freeHand  (1,1) {mustBeBoolean(kwargs.freeHand)} = false
    kwargs.freeHandFilter (1,1) {mustBeBoolean(kwargs.freeHandFilter)} = false
end
nROI = kwargs.nROI;

if iscell(nROI)
    % In the case that selections are passed to the function, we need to
    % check if they are the same size as the data (i.e. different binning)
    for iSel = 1:size(nROI, 2)
        nROI{iSel} = re_bin(nROI{iSel}, data);
    end
else
    % pick n areas from the QDM DATA for calculation
    msg = sprintf('pick your masks. press ESC to exit');
    logMsg('INPUT',msg,1,0);
    
    if kwargs.freeHand
        nROI = pick_area(data);
    else
        nROI = pick_box(data, 'closeFig', 1); % each Sel in the size of fixed data
    end
end
data(abs(data) >= 5) = nan;

%%
% CREATE MASK FROM SELECTIONS
%
% if freeHand and freeHandFilter is true then you can draw the mask
% directly in the image
if all([kwargs.freeHandFilter, kwargs.freeHand])
    nMasks = nROI;
% otherwise the mask will be calculated from the selection
else
    nMasks = {};
    for iSelect = 1: size(nROI, 2)
        % limit the data to only the selected region all other values set
        % to 0
        selData = data .* nROI{iSelect};
%         pcolor(selData);
%         axis xy; axis equal; axis tight; shading flat;

        % The masked data now gets filtered to create the final mask
        iMaskData = selData >= selectionThreshold * max(selData, [], 'all','omitnan');
        m = limit_mask(nROI{iSelect});
        msg = sprintf('creating mask #%i containing %i/%i pixel (%.2f %%)', iSelect, numel(nonzeros(iMaskData)));
        logMsg('debug',msg,1,0);

        % set mask
        nMasks{end+1} = iMaskData;
    end
end