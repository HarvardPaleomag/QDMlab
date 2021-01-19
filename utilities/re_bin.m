function reBinned = re_bin(data, refData)
% RE_BIN corrects binning for not equally binned datasets.
% 
% Example
% -------
% 4 binned blank should be subtracted from 2 binned data. ->
%
% >> 2binned_data - 4binned_blank (will  cause an error)
% >> 2binned_data - re_bin(4binned_blank, 2binned_data)
% recalculates the correct binning and interpolates. Does not cause error
% 
% Note:
%     This should generally be avoided. Try to use the correct binning

if size(refData) ~= size(data)
    binning = (size(refData,1) / size(data,1));
    disp(['<>   WARNING: binning (' num2str(binning) ') difference detected -> correcting'])
    resize_binning = affine2d([binning, 0, 0; 0, binning, 0; 0, 0, 1]);
    reBinned = imwarp(data, resize_binning);
else
    reBinned = data;
end