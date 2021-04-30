function tform_data = tform_data(data, transForm, refFrame)
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
st = dbstack; funcName = st.name; % get functionName

if refFrame.ImageSize ~= size(data)
    binning = detect_binning(data, 'refFrame', refFrame);
    msg = ['binning (' num2str(binning) ') detected correcting the tform'];
    logMsg('info',funcName,msg,1,0);
    [transForm, refFrame] = tform_bin_down(transForm, refFrame, binning);
end

tform_data = imwarp(data, transForm, 'OutputView', refFrame);
