function [transForm, refFrame] = get_image_tform2(fixedData, movingData, varargin)
% takes reference data and calculate tform, rframe that tranforms target data
% to match the reference in the reference frame
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
%     title: char
%         Adds a title to the checkPlot

inParse = inputParser;
addRequired(inParse, 'fixedData');
addRequired(inParse, 'movingData');
addParameter(inParse, 'checkPlot', false, @islogical);
addParameter(inParse, 'title', 'checkPlot alignment', @ischar);
parse(inParse, fixedData, movingData, varargin{:});

if fixedData == movingData
    disp('<> INFO: Transformation not needed. Same image detected')
    transForm = affine2d([[1 0 0]; [0 1 0]; [0 0 1]]);
    refFrame = imref2d(size(fixedData));
    return
end
    
    fixedData = uint8(255 * mat2gray(fixedData));
    movingData = uint8(255 * mat2gray(movingData));
    
    disp(['<> detecting features in fixed and moving data']) 
    ptsOriginal  = detectSURFFeatures(fixedData);
    ptsDistorted = detectSURFFeatures(movingData);
    [featuresOriginal,validPtsOriginal] = extractFeatures(fixedData,ptsOriginal);
    [featuresDistorted,validPtsDistorted] = extractFeatures(movingData,ptsDistorted);
    
    %% matching features
    disp(['<> matching features in fixed and moving data']) 
    index_pairs = matchFeatures(featuresOriginal,featuresDistorted);
    matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
    matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));

    
    %% remove outliers
    [transForm,inlierPtsDistorted,inlierPtsOriginal] = ...
        estimateGeometricTransform(matchedPtsDistorted,matchedPtsOriginal,...
        'similarity');
    refFrame = imref2d(size(fixedData));

    if inParse.Results.checkPlot
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
    
    
    