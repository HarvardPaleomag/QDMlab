function [transForm, refFrame] = get_image_tform(fixedData, movingData, varargin)
%{
takes reference data and calculate tform, rframe that tranforms target data
to match the reference in the reference frame

DEPRECATED - superseeded by get_image_tform2

parameters:
    fixedData: QDM/LED data
    movingData: QDM/LED data
        data to be matched to the refernce data

optional parameters:
    checkPlot: bool 
        default: false
        Adds a plot to check alignment if true
    title: char
        Adds a title to the checkPlot
%}
% Last change: April 21, 2020: Mike
p = inputParser;
addRequired(p, 'fixedData');
addRequired(p, 'movingData');
addParameter(p, 'checkPlot', false, @islogical);
addParameter(p, 'title', 'checkPlot alignment', @ischar);
parse(p, fixedData, movingData, varargin{:});

if fixedData == movingData
    disp('<> INFO: Transformation not needed. Same image detected')
    transForm = affine2d([[1 0 0]; [0 1 0]; [0 0 1]]);
    refFrame = imref2d(size(fixedData));
    return
end


disp('<> CALCULATING: Transformation, Used for imwarp(data, tform,''OutputView'', imref2d).  ')
transForm = imregtform(movingData, fixedData, 'similarity', optimizer, metric);
refFrame = imref2d(size(fixedData));

% do checkPlot
if p.Results.checkPlot
    transformedData = tform_data(movingData, transForm, refFrame); %imwarp(moving_data, tform,'OutputView', rframe);

    checkFigure = figure('Name', 'Align images');

    subplot(2, 1, 1)
    imshowpair(fixedData, movingData, 'Scaling', 'joint')
    title('uncorrected')
    axis xy; axis equal, axis tight

    subplot(2, 1, 2)
    imshowpair(fixedData, transformedData, 'Scaling', 'joint')
    title('corrected')
    axis xy; axis equal, axis tight
end