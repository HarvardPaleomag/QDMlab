function tform_data = tform_data(data, transForm, refFrame, kwargs)
%[tform_data] = tform_data(data, transForm, refFrame; 'binning')
%{
Convenience function to transforms data into a different reference frame
(e.g. 100G data -> NRM).

parameters:
    data: data to transform, may be QDM data or QDM LED image
    tform: affine2d
    rframe: reference frame for transformation
        %todo: 'none' if the data should just be transformed.
%}
% Last change: April 21, 2020: Mike
arguments
    data
    transForm
    refFrame
    kwargs.binning = 'none'
end

if ~isequal(kwargs.binning, false) & refFrame.ImageSize ~= size(data)
    if isequal(kwargs.binning,'none')
        kwargs.binning = detect_binning(data, 'refFrame', refFrame);
    end
    msg = ['binning (' num2str(kwargs.binning) ') detected correcting the tform'];
    logMsg('debug',msg,1,0);
    [transForm, refFrame] = tform_bin_down(transForm, refFrame, kwargs.binning);
end

tform_data = imwarp(data, transForm, 'OutputView', refFrame);
