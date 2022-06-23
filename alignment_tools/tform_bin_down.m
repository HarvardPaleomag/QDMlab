function [transForm, refFrame] = tform_bin_down(transForm, refFrame, bin)
%[transForm, refFrame] = tform_bin_down(transForm, refFrame, bin)
%
% Transforms the LED image down to binned data.
% 
% 
% Parameters
% ----------
%     tform: affine2d
%     binning: int
% 
% 
% Example
% -------
%     A binning of 2, the output will decrease the size of an
%     aray from (1200x1920) -> (600x960)

[transForm, refFrame] = tform_bin_up(transForm, refFrame, 1/bin);