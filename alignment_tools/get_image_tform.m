function [transForm, refFrame] = get_image_tform(fixedData, movingData, kwargs)
%[transForm, refFrame] = get_image_tform(fixedData, movingData; <checkPlot>, <sharpen>, <title>)
% takes reference data and calculate tform, rframe that tranforms target data
% to match the reference in the reference frame
% 
% Parameters
% ----------
%     fixedData: double
%         Reference QDM/LED data
%     movingData: double
%         QDM/LED data to be matched to the reference data
%     checkPlot: bool [false]
%         Adds a plot to check alignment if true
%     title: char ['checkPlot alignment']
%         Adds a title to the checkPlot

arguments
    fixedData
    movingData
    kwargs.checkPlot {mustBeBoolean(kwargs.checkPlot)}= false
    kwargs.sharpen = false;
    kwargs.title {ischar} = 'checkPlot alignment';
end

if fixedData == movingData
    msg = sprintf('Transformation not needed. Same image detected');
    logMsg('info',msg,1,0);
    transForm = affine2d([[1 0 0]; [0 1 0]; [0 0 1]]);
    refFrame = imref2d(size(fixedData));
    return
end
    
fixedData = uint8(255 * mat2gray(fixedData));
fixedData = imadjust(fixedData);
movingData = uint8(255 * mat2gray(movingData));
movingData = imadjust(movingData);

if kwargs.sharpen
    PSF = fspecial('gaussian',7,10);
    fixedData = deconvblind(fixedData,PSF);
    movingData = deconvblind(movingData,PSF);
end
    
msg = sprintf('detecting features in fixed and moving data');
logMsg('debug',msg,1,1);

ptsOriginal  = detectSURFFeatures(fixedData);
ptsDistorted = detectSURFFeatures(movingData);
[featuresOriginal,validPtsOriginal] = extractFeatures(fixedData,ptsOriginal);
[featuresDistorted,validPtsDistorted] = extractFeatures(movingData,ptsDistorted);

%% matching features
msg = sprintf('matching features in fixed and moving data');
logMsg('debug',msg,1,1);

index_pairs = matchFeatures(featuresOriginal,featuresDistorted);
matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));


%% remove outliers
[transForm,inlierPtsDistorted,inlierPtsOriginal] = ...
    estimateGeometricTransform(matchedPtsDistorted,matchedPtsOriginal,...
    'similarity');
refFrame = imref2d(size(fixedData));

if kwargs.checkPlot
    transformedData = imwarp(movingData,transForm, 'OutputView', refFrame);

    checkFigure = figure('Position', [610,150,378,835]);

    subplot(3,1, 1)
    imshowpair(fixedData, movingData, 'Scaling', 'joint')
    title('uncorrected')
    axis xy; axis equal, axis tight
    a = gca;

    subplot(3,1,2)
    showMatchedFeatures(fixedData,movingData, inlierPtsOriginal,inlierPtsDistorted);
    title('Matched features');
    axis xy; axis equal, axis tight
    b = gca;

    subplot(3,1, 3)
    imshowpair(fixedData, transformedData, 'Scaling', 'joint')
    title('corrected')
    axis xy; axis equal, axis tight
    c = gca;

    linkaxes([a b c])
end

    
    