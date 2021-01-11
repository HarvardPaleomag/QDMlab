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
    kwargs.reverse  (1,1) {mustBeMember(kwargs.reverse, [1, 0])} = 0
    kwargs.freeHand  (1,1) {mustBeMember(kwargs.freeHand, [1, 0])} = 0
    kwargs.freeHandFilter (1,1) {mustBeMember(kwargs.freeHandFilter, [1, 0])} = 0
    kwargs.selectionThreshold (1,1) {mustBeNumeric} = 0.25
    kwargs.checkPlot  (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0
    kwargs.nROI = 0
    kwargs.chi = 0
    kwargs.winSize (1,1) {mustBeNumeric} = 4
    kwargs.bootStrapError = 1
    kwargs.bootStrapPixels = 0
end

% define optional function parameters
fileName = kwargs.fileName;

if ~endsWith(fileName, '.mat')
    fileName = [fileName, '.mat'];
end

transFormFile = kwargs.transFormFile;
fixedIdx = kwargs.fixedIdx;
reverse = kwargs.reverse;
checkPlot = kwargs.checkPlot;
selectionThreshold = kwargs.selectionThreshold;
nROI = kwargs.nROI;
removeHotPixels = kwargs.removeHotPixels;
upCont = kwargs.upCont;
chi = kwargs.chi;
winSize = kwargs.winSize;
bootStrapError = kwargs.bootStrapError;
bootStrapPixels = kwargs.bootStrapPixels;
%%
% fix shape for nFolders
nFolders = correct_cell_shape(nFolders);
%%
% define QDM parameters
pixelsize = 4.68e-6;

% ERROR estimation
nMonteCarlo = 1000; % number of shifts around the mask to estimate errors
pixelError = 3; % maximum number of pixel that the alignment could be wrong

% generate reference file name
fixedFile = [nFolders{fixedIdx}, filesep, fileName];

%%
% get transformations and rframes
[nTransForms, nRefFrames] = get_tform_multi(fixedFile, nFolders, ...
    'transFormFile', transFormFile, 'reverse', reverse, 'checkPlot', checkPlot);

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

if removeHotPixels
    if chi
        chi = refFile.chi2Pos1 + refFile.chi2Pos2 + refFile.chi2Neg1 + refFile.chi2Neg2;
    end

    fixedData = filter_hot_pixels(fixedData, 'cutOff',removeHotPixels, ...
        'includeHotPixel',false, 'checkPlot', checkPlot);
end

% detect binning
binning = size(fixedLed, 1) / size(fixedData, 1);
resizeBinning = affine2d([1/binning, 0, 0; 0, 1/binning, 0; 0, 0, 1]);

%%
% The transformation was done on the LED image. Therefore, if there is
% binning involved, the tform 2daffine object needs to be corrected for
% that.
if binning ~= 1
    disp(['<> WARNING: differently sized data/led detected. ', ...
        'Changing all transForms to match data: ', num2str(binning)])
    fnames = keys(nTransForms);
    for i = 1:size(nTransForms)
        fname = fnames{i};
        [tf, rf] = tform_bin_down(nTransForms(fname), nRefFrames(fname), binning);
        nTransForms(fname) = tf;
        nRefFrames(fname) = rf;
    end
end

%%
% transForms and refFrames are now in data coordinates
if iscell(nROI)
    % In the case that selections are passed to the function, we need to
    % check if they are the same size as the data (i.e. different binning)
    for iSel = 1:size(nROI, 2)
        if size(nROI{iSel}, 1) ~= size(fixedData, 1)
            binSel = size(fixedData, 1)/size(nROI{iSel}, 1);
            resizeBinningSel = affine2d([binSel, 0, 0; 0, binSel, 0; 0, 0, 1]);
            newSel = imwarp(nROI{iSel}, resizeBinningSel);
            nROI{iSel} = newSel;
        end
    end
else
    % pick n areas from the QDM DATA for calculation
    disp('<> pick masks')
    if kwargs.freeHand
        nROI = pick_area(fixedData);
    else
        nROI = pick_box(fixedData); % each Sel in the size of fixed data
    end
end
%%
% CREATE MASK FROM SELECTIONS
%
% if freeHand and freeHandFilter is true then you can draw the mask
% directly in the image
if all([kwargs.freeHandFilter, kwargs.freeHand])
    nMasks = nROI;
% otherwise the mask will be calculated from the selection
else
    nMasks = {};
    for iSelect = 1: size(nROI, 2)
        % limit the data to only the selected region all other values set
        % to 0
        selFixedData = fixedData .* nROI{iSelect};

        
        % The masked data now gets filtered to create the final mask
        iMaskData = selFixedData >= selectionThreshold * max(selFixedData, [], 'all', 'omitnan');
        disp(['<> creating mask #', num2str(i), ' containing ', num2str(numel(nonzeros(iMaskData))), ' Pixel'])

        % set mask
        nMasks{end+1} = iMaskData;
    end
end

%% tranformation / filtering
transformedData = containers.Map;
nFiles = {};
% cycle through all folders
for i = 1:size(nFolders, 2)
    % create filename
    iFolder = nFolders{i};
    iFile = fullfile(iFolder, filesep, fileName);
    iFile = check_suffix(iFile);
    nFiles{end+1} = iFile;
    
    disp(['<> loading: target file for transformation: '])
    disp(['      << ', iFile '>>'])
    target = load(iFile);

    if contains(fileName, 'B111')
        targetData = target.B111ferro;
        targetLed = target.ledImg;
    else
        targetData = target.Bz;
        targetLed = target.newLED;
    end

    %% upward cont.
    if iscell(upCont)
        h = upCont{i};
        if std(abs(targetData(:))) > mean(abs(targetData(:)))
            disp('<>   HOT pixels found: filtering before upward continuation')
            targetData = filter_hot_pixels(targetData, 'cutOff', 10);
        end

        if h ~= 0
            disp(['<>   calculating upward continuation (' num2str(h) ') micron'])
            targetData = UpCont(targetData, h*1e-6, 1/pixelsize);
        end
    end

    %% filtering
    if removeHotPixels
        if isnumeric(chi)
            chi = target.chi2Pos1 + target.chi2Pos2 + target.chi2Neg1 + target.chi2Neg2;
        end
        disp(['<>   filtering: ...' iFile(end-40:end-20)  '... .mat'])

        targetData = filter_hot_pixels(targetData, 'cutOff',removeHotPixels, ...
            'includeHotPixel',false, 'winSize', winSize, 'checkPlot', checkPlot,'chi', chi);
    end

    iTransForm = nTransForms(iFile);
    iRefFrame = nRefFrames(iFile);

    %% reverse
    % in the case of reverse, tform and rframe are the ref -> target
    % tranformation. Therefore, the data/led of the target does not need to
    % be transformed. However, the mask itself needs to be transformed form
    % the reference coordinates to the target coordinates later.

    if reverse
        transData = targetData;
        transLed = targetLed;
    else
        disp(['<>    transforming: target data/LED'])
        disp(['         << ...', fileName, ' >>'])
        transData = tform_data(targetData, iTransForm, iRefFrame);
        transLed = tform_data(targetLed, iTransForm, iRefFrame);
    end

    if iscell(upCont)
        h = upCont{i};
        disp(['<> calculating upward continuation (', num2str(h), ') micron'])
        transData = UpCont(transData, h*1e-6, 1/pixelsize);
    end
    % create struct for the the transformed data of this file
    fileTransForm = struct;
    fileTransForm.transData = transData;
    fileTransForm.transLed = transLed;
    fileTransForm.transForm = iTransForm;
    fileTransForm.refFrame = iRefFrame;

    % save the result in the trans_data container for later use
    transformedData(iFile) = fileTransForm;
end

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

    for i = 1:size(nMasks, 2)
        iMask = nMasks{:, i};

        if reverse
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
        
        % calculate parameters
        for n = 1:bootStrapError
            % create masked_data: mask is array with 0 where is should be
            % masked and 1 where it should not
            dx = randi([-bootStrapPixels, bootStrapPixels]);
            dy = randi([-bootStrapPixels, bootStrapPixels]);
            mask = shift_matrix(iMask, dx, dy);
            mData = mask .* iFileData.transData;
            mDataCut = limit_mask(mData);

            nPixel = numel(nonzeros(d0Cut));
            pPixel = [pPixel, numel(nonzeros(sign(nonzeros(mDataCut))+1))];
            pPixelRat = [pPixelRat, pPixel / nPixel];
            sumDiff = [sumDiff, sum(sum(mDataCut-d0Cut))];
            err = [err, immse(mDataCut, d0Cut)];
            maskSum = [maskSum, sum(mDataCut, 'all')];

            % Pixel above threshold
            thresh = mDataCut >= selectionThreshold * max(mDataCut, [], 'all');
            iPixelThresh = [iPixelThresh, numel(nonzeros(thresh))];
        end
        
        % Pixel above threshold
        thresh = mDataCut >= mean(d0Select,'all') + selectionThreshold * std(d0Select, 0, 'all');
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

if checkPlot
    coercivity_result_plot(results)
end