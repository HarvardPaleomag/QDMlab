function [results, files, nROI, nMasks] = estimate_coercivity(nFolders, kwargs)
% These codes (1) register the maps and (2) analizes a user selected magnetic
% pattern for changes from one map to the next.(folders, varargin)
% 
% Parameters
% ----------
%     positional
%     ----------
%     folders:list
%         list of absolute path to each data folder. First entry in the list
%         is used as reference image.
%         Note: calculations in the order the folders are given
% 
%     optional parameters
%     -------------------
%     fileName: str
%         name of the .mat file to be used.
%         default: 'Bz_uc0.mat'
%     transFormFile: str
%         default: 'none'
%         absolute path of a file where previously calculated transformations
%         are in. If file does not exist the file will be created so that
%         you dont have to do this every time.
%     removeHotPixel: bool
%         default: false
%         if True hotpixels will be removed
%         if False hotpixels will not be removed
%     includeHotPixel: bool (0)
%         if 1: hotpixel value is also used to calculate the new value that replaces the hot pixel
%         if 0: only pixels in window with winSize are used to calculate the new value that replaces the hot pixel
%     selectionThreshold: numeric
%         default: 0.5
%         defines the Threshold above which the mask is created.
%         Example:
%             maskSelection = [1 2 0; 1 1 2; 0 1 1]  selectionThreshold = 0.5 ->
%             mask          = [0 1 0; 1 0 1; 0 0 1]
%     upCont: cell
%         Cellarray of upward continuation distances for each of the files.
%         Needs to be in same order.
%     filter: str
%         Uses filter_hot_pixels to clean the maps before calculating the
%         results.
%         Note: 'includeHotPixel' is false by default
%     freeHand: bool
%         default: false
%         Instead of using a rectangular selection use freehand drawn
%     freeHandFilter: bool
%         default: false
%         If true: the selectionThreshold is used to create the mask
%                  from the maskSelection.
%         If false: the mask = maskSelection
%     fixedIdx: int
%         default: 1
%         index of the reference LED image. This will fixed while the other
%         image is transformed.
%     chi:
%         if true the chi2 value is used to filter and not the data itself
%         Note: the chi2 value is calculated from the sum(pos1, pos2, neg1,
%               neg2) chi values.
%     nROI: cell
%         cell with region of interest.
%         Note: the ROI and the selection are not necessarily the same. If
%         the ROI was selected with freeHand then ROI == selection. Otherwise
%         it is a rectangle around the selection.
%     reverse: bool ==> CURRENTLY NOT SUPPORTED
%         default: false
%         if true:  refernce - tform -> target
%         if false: target   - tform -> reference
%     bootStrapError: int, bool
%         default: 1
%         This uses a boot strapping approach to estimate errors. It will
%         shift the mask by 'bootStrapPixels' 'bootStrapError' times. The
%         result value is calculated from the mean of all calculations and an
%         error is estimated.
%         If 'bootStrapError' == 1 only one value is calculated -> no error
%         estimation.
%         Note: if bootStrapError == 1, 'bootStrapPixels' should be 0
%         otherwise the mask is shifted and the result will be wrong.
%     bootStrapPixels: int
%         default: 0
%         The number of pixels the mask will be randomly shifted to estimate
%         an error.
% 
% Returns
% -------
%     struct with 7 cells for each result (iFiles,
%     pPixels, pPixelRats, iPixelThreshs, 
%     sumDiffs, errs, maskSums, iMasks, transData).
%     Each is a cell with result{i,j} with ith mask and jth file:


arguments
    nFolders
    kwargs.fileName = 'Bz_uc0.mat'
    kwargs.transFormFile = 'none'
    kwargs.fixedIdx (1,1) {mustBePositive} = 1
    kwargs.upCont = 0
    kwargs.removeHotPixels = 0
    kwargs.includeHotPixel = 0
    kwargs.reverse  (1,1) {mustBeMember(kwargs.reverse, [1, 0])} = 0
    kwargs.freeHand  (1,1) {mustBeMember(kwargs.freeHand, [1, 0])} = 0
    kwargs.freeHandFilter (1,1) {mustBeMember(kwargs.freeHandFilter, [1, 0])} = 0
    kwargs.selectionThreshold (1,1) {mustBeNumeric} = 0.25
    kwargs.checkPlot  (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0
    kwargs.nROI = 0
    kwargs.chi = 0
    kwargs.winSize (1,1) = 4
    kwargs.bootStrapN = 1
    kwargs.bootStrapPixels = 4
end

% define optional function parameters
fileName = kwargs.fileName;
nROI = kwargs.nROI;

if ~endsWith(fileName, '.mat')
    fileName = [fileName, '.mat'];
end
%%
% fix shape for nFolders
nFolders = correct_cell_shape(nFolders);
%%

% ERROR estimation
nMonteCarlo = 1000; % number of shifts around the mask to estimate errors
pixelError = 3; % maximum number of pixel that the alignment could be wrong

% generate reference file name
fixedFile = [nFolders{kwargs.fixedIdx}, filesep, fileName];

%%
%% tranformation / filtering
[transformedData, nFiles] = get_transformed_maps(nFolders, ...
                  'fileName', kwargs.fileName, 'transFormFile', kwargs.transFormFile,...
                  'fixedIdx', kwargs.fixedIdx, 'removeHotPixels', kwargs.removeHotPixels,...
                  'includeHotPixel', kwargs.includeHotPixel,  'chi', kwargs.chi, ...
                  'winSize', kwargs.winSize, 'reverse', kwargs.reverse, ...
                  'checkPlot', kwargs.checkPlot);

% load the reference data
refFile = load(fixedFile);

% check for variable names
if contains(fileName, 'B111')
    fixedData = refFile.B111ferro;
    fixedLed = refFile.ledImg;
else
    fixedData = refFile.Bz;
    fixedLed = refFile.newLED;
end

if kwargs.removeHotPixels
    if kwargs.chi
        chi = refFile.chi2Pos1 + refFile.chi2Pos2 + refFile.chi2Neg1 + refFile.chi2Neg2;
    else
        chi = kwargs.chi;
    end

    fixedData = filter_hot_pixels(fixedData, 'cutOff', kwargs.removeHotPixels, ...
                'chi', chi, 'includeHotPixel',false, 'checkPlot', kwargs.checkPlot);
end

% detect binning
binning = detect_binning(refFile);

[nMasks, nROI] = create_masks(fixedData, kwargs.selectionThreshold,...
                      'nROI', nROI, 'freeHand', kwargs.freeHand,...
                      'freeHandFilter', kwargs.freeHandFilter);


% calculation of the mask for each file.
% Note: this is where the mask is transformed in case of 'reverse' == true
files = cell(numel(nFolders));

% preallocate the cells:
iFiles = {};
pPixels = [];
nPixels = [];
pPixelRats = [];
iPixelThreshs = [];
sumDiffs = [];
errs = [];
maskSums = [];
iMasks = {};
transDatas = {};
transLeds = {};

%% save checkPlot axes to link them later
axes = [];

% if checkPlot
%     checkfig = figure;
% end

for j = 1:size(nFiles, 2)
    iFile = nFiles{j};

    % load the data of this file
    iFileData = transformedData(iFile);
    files{j} = iFileData;

%     disp('<> ------------------------------------------------------------')
    % iterate over each mask
    for i = 1:size(nMasks, 2)
        iMask = nMasks{:, i};

        if kwargs.reverse
            disp(['<>    transforming mask to match << ...', iFile(end-40:end), ' >>'])
            iMask = tform_data(iMask, iFileData.transForm, iFileData.refFrame);
        end

        % create masked_data: mask is array with 0 where is should be
        % masked and 1 where it should not
        mData = iMask .* iFileData.transData;
        mData = mData - nanmedian(mData, 'all');

        % masked reference
        d0Select = crop_data(fixedData, nROI{i});
        d0 = iMask .* fixedData;
        d0Cut = limit_mask(d0);

        % cut around the mask
        mDataCut = crop_data(mData, iMask);
        d0Cut = crop_data(d0, iMask);

%         disp(['<> masking << ...', iFile(end-40:end), ' >>'])

        % predefine the variables
        mData = [];
        pPixel = [];
        pPixelRat = [];
        sumDiff = [];
        err = [];
        maskSum = [];
        iPixelThresh = [];


        if size(mDataCut) ~= size(d0Cut)
            disp('  WARNING mask too close to edge, skipping ... ')
            continue
        end
        
        dx = 0; dy = 0;
        % calculate parameters
        
        for n = 1:kwargs.bootStrapN
            % create masked_data: mask is array with 0 where is should be
            % masked and 1 where it should not
            if kwargs.bootStrapN
                dx = randi([-kwargs.bootStrapPixels, kwargs.bootStrapPixels]);
                dy = randi([-kwargs.bootStrapPixels, kwargs.bootStrapPixels]);
            end
            
            mask = shift_matrix(iMask, dx, dy);
            mData = iFileData.transData .* mask;
            mDataCut = limit_mask(mData);

            nPixel = numel(nonzeros(d0Cut));
            pPixel = [pPixel, numel(nonzeros(sign(nonzeros(mDataCut))+1))];
            pPixelRat = [pPixelRat, pPixel / nPixel];
            sumDiff = [sumDiff, sum(sum(mDataCut-d0Cut))];
            err = [err, immse(mDataCut, d0Cut)];
            maskSum = [maskSum, sum(mDataCut, 'all')];

            % Pixel above threshold
            thresh = mDataCut >= kwargs.selectionThreshold * max(mDataCut, [], 'all');
            iPixelThresh = [iPixelThresh, numel(nonzeros(thresh))];
        end
        
        % Pixel above threshold
        thresh = mDataCut >= mean(d0Select,'all') + kwargs.selectionThreshold * std(d0Select, 0, 'all');
        iPixelThresh = numel(nonzeros(thresh));

        % store everything
        iFiles{i, j} = iFile;
        % number of non 0 pixels in mask
        nPixels(i, j, 1) = mean(nPixel);
        nPixels(i, j, 2) = std(nPixel);
        % number of positive pixels in mask
        pPixels(i, j, 1) = mean(pPixel);
        pPixels(i, j, 2) = std(pPixel);
        % ratio of positive / total pixels in mask
        pPixelRats(i, j, 1) = mean(pPixelRat);
        pPixelRats(i, j, 2) = std(pPixelRat);
        % number of pixels above threshold
        iPixelThreshs(i, j, 1) = mean(iPixelThresh);
        iPixelThreshs(i, j, 2) = std(iPixelThresh);
        % sum of difference between mask and same mask in fixedData
        sumDiffs(i, j, 1) = mean(sumDiff);
        sumDiffs(i, j, 2) = std(sumDiff);
        errs(i, j, 1) = mean(err);
        errs(i, j, 2) = std(err);
        % Sum of the Mask
        maskSums(i, j, 1) = mean(maskSum);
        maskSums(i, j, 2) = std(maskSum);
        % mask itself
        iMasks{i} = iMask;
        % transformed data
        transDatas{j} = iFileData.transData;
        % transformed data
        transLeds{j} = iFileData.transLed;
    end

end

results = struct('nFiles', {iFiles}, 'pPixels', pPixels, 'pPixelRats', pPixelRats, ...
    'nPixelThreshs', iPixelThreshs, 'sumDiffs', sumDiffs, ...
    'errs', errs, 'maskSums', maskSums, 'nPixels', nPixels, ...
    'nMasks', {iMasks}, 'nROI', {nROI},...
    'transDatas', {transDatas}, 'fixedData', fixedData, 'transLeds', {transLeds});

fprintf('<>   INFO: coercivity estimation complete. Output: (%i x %i x 2) = (ROI, file, (value, std)\n', size(nROI,2), size(iFiles,2));

if kwargs.checkPlot
    coercivity_result_plot(results)
end