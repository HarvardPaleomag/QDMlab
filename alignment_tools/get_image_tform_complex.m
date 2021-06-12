function [transForm, refFrame] = get_image_tform_complex(fixedData, movingData, kwargs)
%[transForm, refFrame] = get_image_tform_complex(fixedData, movingData, varargin; 'checkPlot', 'binning')
% Function lets you pick several points on a reference image and the target
% image. It calculates a transformation
% 
% Parameters
% ----------
%     fixedData: double QDM/LED data
%     movingData: double QDM/LED data
%         data to be matched to the refernce data
%     binning: int [1]
%       needed for correct transformation
%     checkPlot: bool [false]
%         Adds a plot to check alignment if true
arguments
    fixedData double
    movingData double
    kwargs.checkPlot {mustBeBoolean(kwargs.checkPlot)} = false 
    kwargs.binning {mustBeInteger(kwargs.binning)} = 1
end 

if max(fixedData, [], 'all') > 1
    fixedData = double(fixedData);
    fixedData = fixedData / max(fixedData, [], 'all');
end

if max(movingData, [], 'all') > 1
    movingData = double(movingData);
    movingData = movingData / max(movingData, [], 'all');
end

[mp, fp] = cpselect(movingData, fixedData, 'Wait', true);

transForm = fitgeotrans(mp, fp, 'similarity');

% create refence frame for fixed image
refFrame = imref2d(size(fixedData));

if isequal(kwargs.checkPlot, true)
    checkFigure = figure('Name', 'Align images');
    
    transformedData = imwarp(movingData, transForm, 'OutputView', refFrame,'binning', false);

    subplot(2, 1, 1)
    imshowpair(fixedData, movingData, 'Scaling', 'joint')
    title('uncorrected')
    axis xy; axis equal, axis tight

    subplot(2, 1, 2)
    imshowpair(fixedData, transformedData, 'Scaling', 'joint')
    title('corrected')
    axis xy; axis equal, axis tight

end