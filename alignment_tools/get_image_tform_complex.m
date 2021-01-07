function [transForm, refFrame] = get_image_tform_complex(fixedData, movingData, varargin)
% Function lets you pick several points on a reference image and the target
% image. It calculates a transformation
% 
% parameters:
%     fixedData: QDM/LED data
%     movingData: QDM/LED data
%         data to be matched to the refernce data
% 
% optional parameters:
%     checkPlot: bool 
%         default: false
%         Adds a plot to check alignment if true
p = inputParser;
addRequired(p, 'fixedData');
addRequired(p, 'movingData');
addParameter(p, 'checkPlot', false, @islogical);
parse(p, fixedData, movingData, varargin{:});

checkPlot = p.Results.checkPlot;

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

% create refence frame for fixed image //todo fix refFrame
refFrame = imref2d(size(fixedData));

if checkPlot == true
    checkFigure = figure('Name', 'Align images');
    
    transformedData = tform_data(movingData, transForm, refFrame); %imwarp(moving_data, tform,'OutputView', rframe);

    subplot(2, 1, 1)
    imshowpair(fixedData, movingData, 'Scaling', 'joint')
    title('uncorrected')
    axis xy; axis equal, axis tight

    subplot(2, 1, 2)
    imshowpair(fixedData, transformedData, 'Scaling', 'joint')
    title('corrected')
    axis xy; axis equal, axis tight
%     movedData = imwarp(movingData, transForm, 'OutputView', refFrame);
%     imshowpair(fixedData, movedData, 'blend');
end