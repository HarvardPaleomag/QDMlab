function [transForm, refFrame] = tform_bin_up(transForm, refFrame, bin)
%{
Transforms the binned data up to LED. Calls tform_led2data.

parameters:
    tform: affine2d
    binning: int

example:
    A binning of 2, the output will decrease the size of an
    aray from (600x960) -> (1200x1920)
%}
% Last change: April 21, 2020: Mike

binTransForm = [[1, 1, 1]; [1, 1, 1]; [bin, bin, 1]];
transForm.T = transForm.T .* binTransForm;

refFrame.YWorldLimits= refFrame.YWorldLimits * bin;
refFrame.XWorldLimits= refFrame.XWorldLimits * bin;
refFrame.ImageSize = refFrame.ImageSize * bin;
