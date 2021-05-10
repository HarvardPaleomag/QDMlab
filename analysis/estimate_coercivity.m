function [results, files, nROI, nMasks] = estimate_coercivity(nFolders, kwargs, selection, error, filter)
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
%     removeHotPixel: bool [0]
%         if True hotpixels will be removed
%         if False hotpixels will not be removed
%     includeHotPixel: bool [0]
%         | if 1: hotpixel value is also used to calculate the new value that replaces the hot pixel
%         | if 0: only pixels in window with winSize are used to calculate the new value that replaces the hot pixel
%     selectionThreshold: numeric [0.25]
%         | defines the Threshold above which the mask is created.
%         | **Example:** selectionThreshold = 0.5
%         | :code:`maskSelection = [1 2 0; 1 1 2; 0 1 1]`  ->
%         | :code:`mask          = [0 1 0; 1 0 1; 0 0 1]`
%     upCont: cell
%         Cellarray of upward continuation distances for each of the files.
%         Needs to be in same order.
%     filter: str
%         Uses filter_hot_pixels to clean the maps before calculating the
%         results.
%         Note: 'includeHotPixel' is false by default
%     freeHand: bool [0]
%         Instead of using a rectangular selection use freehand drawn
%     freeHandSelection: bool [0]
%         | If true: the selectionThreshold is used to create the mask from the maskSelection.
%         | If false: the mask = maskSelection
%     fixedIdx: int [1]
%         index of the reference LED image. This will fixed while the other
%         image is transformed.
%     chi:
%         | if true the chi2 value is used to filter and not the data itself
%         | **Note:** the chi2 value is calculated from the sum(pos1, pos2, neg1,neg2) chi values.
%     nROI: cell
%         | cell with region of interest.
%         | **Note:** the ROI and the selection are not necessarily the same.|
%           If the ROI was selected with freeHand then ROI == selection.
%           Otherwise it is a rectangle around the selection.
%     reverse: bool ==> CURRENTLY NOT SUPPORTED
%         default: false
%         if true:  refernce - tform -> target
%         if false: target   - tform -> reference
%     bootStrapError: int, bool [1]
%         This uses a boot strapping approach to estimate errors. It will
%         shift the mask by 'pixelError' 'bootStrapError' times. The
%         result value is calculated from the mean of all calculations and an
%         error is estimated.
%         If 'bootStrapError' == 1 only one value is calculated -> no error
%         estimation.
%         Note: if bootStrapError == 1, 'pixelError' should be 0
%         otherwise the mask is shifted and the result will be wrong.
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
%         >> estimate_coercivity(nFolders)
%
%     This calls the function with the default parameters (the values in square brackets *[]*).     


arguments
    nFolders cell {foldersMustExist(nFolders)}
    kwargs.fileName char {fileMustExistInFolder(kwargs.fileName, nFolders)} = 'Bz_uc0'
    kwargs.transFormFile = 'none'
    kwargs.fixedIdx (1,1) {mustBePositive} = 1
    kwargs.upCont = num2cell(zeros(size(nFolders,1)))
    
    kwargs.checkPlot  (1,1) {mustBeBoolean(kwargs.checkPlot)} = 0
    kwargs.reverse  (1,1) {mustBeBoolean(kwargs.reverse)} = 0
    kwargs.nROI (1,1) {mustBeBoolean(kwargs.nROI)} = false
    
    selection.freeHand  (1,1) {mustBeBoolean(selection.freeHand)} = false
    selection.freeHandSelection (1,1) {mustBeBoolean(selection.freeHandSelection)} = false
    selection.selectionThreshold (1,1) {mustBeNumeric} = 0.25
    
    error.bootStrapN = 1
    error.pixelError = 4
    
    filter.removeHotPixels = 0
    filter.threshold = 5
    filter.includeHotPixel = 0
    filter.chi = 0
    filter.winSize (1,1) = 4
    


end

% define optional function parameters
fileName = kwargs.fileName;
fileName = check_suffix(fileName);

nROI = kwargs.nROI;

%%
% fix shape for nFolders
nFolders = correct_cell_shape(nFolders);
%%

% generate reference file name
fixedFile = [nFolders{kwargs.fixedIdx}, filesep, fileName];

%%
%% tranformation / filtering
[transformedData, nFiles] = get_transformed_maps(nFolders, ...
                  'fileName', kwargs.fileName, 'transFormFile', kwargs.transFormFile,...
                  'fixedIdx', kwargs.fixedIdx, 'reverse', kwargs.reverse, ...
                  'upCont', kwargs.upCont, 'checkPlot', kwargs.checkPlot, ...
                  'removeHotPixels', filter.removeHotPixels, 'winSize', filter.winSize, ...
                  'includeHotPixel', filter.includeHotPixel, 'chi', filter.chi, ...
                  'threshold', filter.threshold);

% load the reference data
refFile = load(fixedFile);
[bool, dataName, ledName] = is_B111(refFile);

% read data and threshold to 5
fixedData = refFile.(dataName);
fixedData = filter_hot_pixels(fixedData, 'threshold', filter.threshold);
% read LED
fixedLed = refFile.(ledName);

if filter.removeHotPixels
    if filter.chi
        chi = refFile.chi2Pos1 + refFile.chi2Pos2 + refFile.chi2Neg1 + refFile.chi2Neg2;
    else
        chi = filter.chi;
    end

    fixedData = filter_hot_pixels(fixedData, 'cutOff', filter.removeHotPixels, ...
                'chi', chi, 'includeHotPixel',false, 'checkPlot', filter.checkPlot);
end



[nMasks, nROI] = create_masks(fixedData, selection.selectionThreshold,...
                      'nROI', nROI, 'freeHand', selection.freeHand,...
                      'freeHandSelection', selection.freeHandSelection);


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
            msg = ['transforming mask to match << ...', iFile(end-40:end), ' >>'];
            logMsg('info','dipole_fit',msg,1,0);
            iMask = tform_data(iMask, iFileData.transForm, iFileData.refFrame);
        end

        % create masked_data: mask is array with 0 where is should be
        % masked and 1 where it should not
        mData = iMask .* iFileData.transData;
        mData = mData - nanmedian(mData, 'all');

        % masked reference
        d0Select = crop_data(fixedData, nROI{i});
        d0 = iMask .* fixedData;

        % cut around the mask
        mDataCut = crop_data(mData, iMask);
        d0Cut = crop_data(d0, iMask);

%         disp(['<> masking << ...', iFile(end-40:end), ' >>'])

        % predefine the variables
        pPixel = [];
        pPixelRat = [];
        sumDiff = [];
        err = [];
        maskSum = [];
        iPixelThresh = [];


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
            if error.bootStrapN
                dx = randi([-error.pixelError, error.pixelError]);
                dy = randi([-error.pixelError, error.pixelError]);
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
            thresh = mDataCut >= selection.selectionThreshold * max(mDataCut, [], 'all');
            iPixelThresh = [iPixelThresh, numel(nonzeros(thresh))];
        end

        % Pixel above threshold
        thresh = mDataCut >= mean(d0Select,'all') + selection.selectionThreshold * std(d0Select, 0, 'all');
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

msg = sprintf('coercivity estimation complete. Output: (%i x %i x 2) = (ROI, file, (value, std)', size(nROI,2), size(iFiles,2));
logMsg('finished',msg,1,0);

if kwargs.checkPlot
    coercivity_result_plot(results)
end
