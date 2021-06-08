function [results, files, nROI, nMasks] = demag_behavior(nFolders, kwargs, selection, error, filter)
%[results, files, nROI, nMasks] = demag_behavior(nFolders; 'fileName', 'transFormFile', 'fixedIdx', 'upCont', 'checkPlot', 'reverse', 'nROI', 'freeHand', 'freeHandSelection', 'selectionThreshold', 'bootStrapN', 'pixelShift', 'filterProps')
% These codes (1) register the maps and (2) analizes a user selected magnetic
% pattern for changes from one map to the next.(folders, varargin)
%
% Parameters
% ----------
%     folders: list
%         list of absolute path to each data folder. First entry in the list
%         is used as reference image.
%         Note: calculations in the order the folders are given
%     fileName: str ['Bz_uc0.mat']
%         name of the .mat file to be used.
%     transFormFile: str ['none']
%         absolute path of a file where previously calculated transformations
%         are in. If file does not exist the file will be created so that
%         you dont have to do this every time.
%     selectionThreshold: numeric [0.25]
%         defines the Threshold above which the mask is created.
%         **Example:** selectionThreshold = 0.5
%         :code:`maskSelection = [1 2 0; 1 1 2; 0 1 1]`  ->
%         :code:`mask          = [0 1 0; 1 0 1; 0 0 1]`
%     upCont: cell
%         Cellarray of upward continuation distances for each of the files.
%         Needs to be in same order.
%     freeHand: bool [0]
%         Instead of using a rectangular selection use freehand drawn
%     freeHandSelection: bool [0]
%         If true: the selectionThreshold is used to create the mask from the maskSelection.
%         If false: the mask = maskSelection
%     fixedIdx: int [1]
%         index of the reference LED image. This will fixed while the other
%         image is transformed.
%     chi:
%         if true the chi2 value is used to filter and not the data itself
%         **Note:** the chi2 value is calculated from the sum(pos1, pos2, neg1,neg2) chi values.
%     nROI: cell
%         cell with region of interest.
%         **Note:** the ROI and the selection are not necessarily the same.|
%         If the ROI was selected with freeHand then ROI == selection.
%         Otherwise it is a rectangle around the selection.
%     bootStrapN: int, bool [1]
%         This uses a boot strapping approach to estimate errors. It will
%         shift the mask by 'pixelError' 'bootStrapN' times. The
%         result value is calculated from the mean of all calculations and an
%         error is estimated.
%         If 'bootStrapError' == 1 only one value is calculated -> no error
%         estimation.
%     pixelError: int [0]
%         Used for bootstrapping an error estimate.
%         The number of pixels the mask will be randomly shifted to estimate
%         an error.
%
% Returns
% -------
%     struct with 7 cells for each result (iFiles,
%     pPixels, pPixelRats, iPixelThreshs,
%     sumDiffs, errs, maskSums, iMasks, transData).
%     Each is a cell with result{i,j} with ith mask and jth file:
%
% hint
% ----
%     **Example:** Here we have a cell of folders, containing QDM images ('B111dataToPlot', 'Bz_uc0').
%     literal blocks::
%         >> nFolders = {'home/NRM, 'home/10mT','home/20mT','home/30mT'};
%
%     The most simple way to estimate the coercivity is to call:
%     literal blocks::
%         >> demag_behavior(nFolders)
%
%     This calls the function with the default parameters (the values in square brackets *[]*).     


arguments
    nFolders cell {foldersMustExist(nFolders)}
    kwargs.fileName char {fileMustExistInFolder(kwargs.fileName, nFolders)} = 'Bz_uc0'
    kwargs.transFormFile = 'none'
    kwargs.fixedIdx (1,1) {mustBePositive} = 1
    kwargs.upCont = false;
    
    kwargs.checkPlot  (1,1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.reverse  (1,1) {mustBeBoolean(kwargs.reverse)} = false;
    kwargs.nROI = false
    
    selection.freeHand  (1,1) {mustBeBoolean(selection.freeHand)} = false
    selection.freeHandSelection (1,1) {mustBeBoolean(selection.freeHandSelection)} = false
    selection.selectionThreshold (1,1) {mustBeNumeric} = 0.25
    
    error.bootStrapN (1,1) {mustBeNumeric} = 1
    error.pixelShift (1,1) {mustBeNumeric} = 4
    
    filter.filterProps struct = struct();

end
%%
% fix shape for nFolders
nFolders = correct_cell_shape(nFolders);

% define optional function parameters
fileName = kwargs.fileName;
fileName = check_suffix(fileName);

nROI = kwargs.nROI;

if ~kwargs.upCont
    kwargs.upCont = num2cell(zeros(1, size(nFolders,2)));
end

% generate reference file name
fixedFile = [nFolders{kwargs.fixedIdx}, filesep, fileName];

%% load the reference data
refFile = load(fixedFile);
[~, dataName, ledName] = is_B111(refFile);
% read data and threshold to 5
fixedData = refFile.(dataName);
% read LED
fixedLed = refFile.(ledName);

%%
if ~all( structfun(@isempty, filter.filterProps))
    filterProps = namedargs2cell(filter.filterProps);
    fixedData = filter_hot_pixels(fixedData, filterProps{:});
end

%% tranformation / filtering
[transformedData, nFiles] = get_transformed_maps(nFolders, ...
                  'fileName', kwargs.fileName, 'transFormFile', kwargs.transFormFile,...
                  'fixedIdx', kwargs.fixedIdx, 'reverse', kwargs.reverse, ...
                  'upCont', kwargs.upCont, 'checkPlot', kwargs.checkPlot, ...
                  'filterProps', filter.filterProps);

[nMasks, nROI] = create_masks(fixedData, selection.selectionThreshold,...
                      'nROI', nROI, 'freeHand', selection.freeHand,...
                      'freeHandSelection', selection.freeHandSelection);


% calculation of the mask for each file.
% Note: this is where the mask is transformed in case of 'reverse' == true
files = cell(numel(nFolders));

% preallocate the cells:
iFiles = {};
pPixels = zeros(size(nMasks, 2), size(nFiles, 2), 2);
nPixels = zeros(size(nMasks, 2), size(nFiles, 2), 2);
pPixelRats = zeros(size(nMasks, 2), size(nFiles, 2), 2);
iPixelThreshs = zeros(size(nMasks, 2), size(nFiles, 2), 2);
sumDiffs = zeros(size(nMasks, 2), size(nFiles, 2), 2);
errs = zeros(size(nMasks, 2), size(nFiles, 2), 2);
maskSums = zeros(size(nMasks, 2), size(nFiles, 2), 2);
iFiles = cell(size(nMasks, 2), size(nFiles, 2));
iMasks = cell(1, size(nMasks, 2));
transDatas = cell(1, size(nFiles, 2));
transLeds = cell(1, size(nFiles, 2));

% if checkPlot
%     checkfig = figure;
% end

for j = 1:size(nFiles, 2)
    iFile = nFiles{j};

    % load the data of this file
    iFileData = transformedData(iFile);
    files{j} = iFileData;

    % iterate over each mask
    for i = 1:size(nMasks, 2)
        iMask = nMasks{:, i};

        if kwargs.reverse
            msg = ['transforming mask to match << ...', iFile(end-40:end), ' >>'];
            logMsg('info','dipole_fit',msg,1,0);
            iMask = tform_data(iMask, iFileData.transForm, iFileData.refFrame);
        end

        % create masked_data: mask is array with 0 where is should be
        % masked and 1 where it should not
        mData = iMask .* iFileData.transData;
%         mData = mData - nanmedian(mData, 'all');

        % masked reference
        d0ROI = crop_data(fixedData, nROI{i});
        d0 = iMask .* fixedData;

        % cut around the ROI/mask
        mDataCut = crop_data(mData, iMask);
        d0Cut = crop_data(d0, iMask);

%         disp(['<> masking << ...', iFile(end-40:end), ' >>'])

        % predefine the variables
        nPixel = zeros(error.bootStrapN,1);
        pPixel = zeros(error.bootStrapN,1);
        pPixelRat = zeros(error.bootStrapN,1);
        sumDiff = zeros(error.bootStrapN,1);
        err = zeros(error.bootStrapN,1);
        maskSum = zeros(error.bootStrapN,1);
        iPixelThresh = zeros(error.bootStrapN,1);


        if size(mDataCut) ~= size(d0Cut)
            msg = sprintf('mask too close to edge, skipping ... ');
            logMsg('warn',msg,1,0);
            continue
        end

        dx = 0; dy = 0;
        % calculate parameters

        for n = 1:error.bootStrapN
            % create masked_data: mask is array with 0 where is should be
            % masked and 1 where it should not
            if error.bootStrapN > 1
                dx = randi([-error.pixelShift, error.pixelShift]);
                dy = randi([-error.pixelShift, error.pixelShift]);
            end
            
            if dx ~= 0 && dy ~= 0
                mask = shift_matrix(iMask, dx, dy);
            else
                mask = iMask;
            end
            
            mData = iFileData.transData .* mask;
            mDataCut = crop_data(mData, mask);

            nPixel(n) = numel(nonzeros(d0Cut));
            pPixel(n) = numel(nonzeros(mDataCut>0));
            pPixelRat(n) = pPixel(n) / nPixel(n);
            sumDiff(n) = sum(sum(mDataCut-d0Cut));
            err(n) = immse(mDataCut, d0Cut);
            maskSum(n) = sum(mDataCut, 'all');

            % Pixel above threshold
            thresh = mDataCut >= selection.selectionThreshold * max(mDataCut, [], 'all');
            iPixelThresh(n) = numel(nonzeros(thresh));
        end

        % Pixel above threshold
        thresh = mDataCut >= mean(d0ROI,'all') + selection.selectionThreshold * std(d0ROI, 0, 'all');
        iPixelThresh = numel(nonzeros(thresh));

        % store everything
        iFiles{i, j} = iFile;
        % number of non 0 pixels in mask
        nPixels(i, j, :) = [mean(nPixel) std(nPixel)];
        % number of positive pixels in mask
        pPixels(i, j, :) = [mean(pPixel) std(pPixel)];
        % ratio of positive / total pixels in mask
        pPixelRats(i, j, :) = [ mean(pPixelRat)  std(pPixelRat)];
        % number of pixels above threshold
        iPixelThreshs(i, j, :) = [mean(iPixelThresh) std(iPixelThresh)];
        % sum of difference between mask and same mask in fixedData
        sumDiffs(i, j, :) = [mean(sumDiff)  std(sumDiff)];
        errs(i, j, :) = [mean(err) std(err)];
        % Sum of the Mask
        maskSums(i, j, :) = [mean(maskSum) std(maskSum)];
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

msg = sprintf('coercivity estimation complete. Output: (%i x %i x 2) = (ROI, file, (value, std)', size(nROI,2), size(iFiles,2));
logMsg('finished',msg,1,0);

if kwargs.checkPlot
    demag_behavior_plot(results)
end
