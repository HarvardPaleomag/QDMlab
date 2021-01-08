function [transForm, refFrame] = tform_bin_down(transForm, refFrame, bin)
%{
Transforms the LED image down to binned data.

parameters:
    tform: affine2d
    binning: int

example:
    A binning of 2, the output will decrease the size of an
    aray from (1200x1920) -> (600x960)
%}
% Last change: April 21, 2020: Mike

[transForm, refFrame] = tform_bin_up(transForm, refFrame, 1/bin);