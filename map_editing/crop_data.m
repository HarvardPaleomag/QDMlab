function croppedData = crop_data(data, mask)
%[croppedData] = crop_data(data, mask)
% Helper function to crop data according to a specified mask
% 
% Parameters
% ----------
%     positional
%     ----------
%     data: array
%         the data
%     mask: logical-array
%         the mask consisting of 0,1

    [x, y, w, h] = get_mask_extent(mask);

    % cut around data
    croppedData = data(y:y+h, x:x+w);
end