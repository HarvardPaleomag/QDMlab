function tform_data = tform_data(data, tform, rframe)
%{
Convenience function to transforms data into a different reference frame
(e.g. 100G data -> NRM).

parameters:
    data: data to transform, may be QDM data or DQM LED image
    tform: affine2d
    rframe: reference frame for transformation
        %todo: 'none' if the data should just be transformed.
%}
% Last change: April 21, 2020: Mike

if rframe.ImageSize ~= size(data)
    binning = (rframe.ImageSize / size(data));
    disp(['<>   Binning (' num2str(binning) ') detected correcting the tform'])
    resize_binning = affine2d([binning, 0, 0; 0, binning, 0; 0, 0, 1]);
    data = imwarp(data, resize_binning);
    rframe.ImageSize = rframe.ImageSize / binning;
end

tform_data = imwarp(data, tform, 'OutputView', rframe);
