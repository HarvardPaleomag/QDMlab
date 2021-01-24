function [transformedData, nFiles] = get_transformed_maps(nFolders, varargin)
% These codes (1) register the maps and (2) analizes a user selected magnetic
% pattern for changes from one map to the next.(folders, varargin)
% 
% Parameters
% ----------
%     positional
%     ----------
%     nFolders:list
%         list of absolute path to each data folder. First entry in the list
%         is used as reference image.
%         Note: calculations in the order the folders are given
% 
%     optional parameters
%     -------------------
%     fileName: str
%         name of the .mat file to be used. 
%         default: 'Bz_uc0'
%     transFormFile: str
%         default: 'none'
%         absolute path of a file where previously calculated transformations
%         are in. If file does not exist the file will be created so that
%         you dont have to do this every time.
%     removeHotPixel: double
%         default: bool, false
%         Uses 'filter_hot_pixel' function. All filtering is done pre
%         tranfornation.
%         if false: hotpixels will not be removed
%         else: removeHotPixel determines how many standard deviations 
%               have to be exceeded for the pixel to be filtered
%     includeHotPixel: bool 
%         default: false
%           uses 'filter_hot_pixel' function
%         if true: hotpixels is included in the calculation of the mean
%         if false: hotpixels is NOT included in the calculation of the mean
%     winSize: int
%         default: 4
%         Window size around the hot pixel for mean calculation (e.g. 4 ->
%         4x4 window)
%     chi: bool
%         if true: chi^2 values are used to filter (sum of target chi^2)
%         if false: data is filtered by data values
%     fixedIdx: int
%         default: 1
%         index of the reference LED image. This will fixed while the other
%         image is transformed.
%     reverse: bool - NOT IMPLEMENTED YET
%         default: false
%         if true:  refernce - tform -> target
%         if false: target   - tform -> reference
%
% See also, filter_hot_pixels
inParse = inputParser;
str_or_char = @(x) isstring(x) | ischar(x);

addRequired(inParse, 'nFolders', @iscell);
addParameter(inParse, 'fileName', 'Bz_uc0', str_or_char);
addParameter(inParse, 'transFormFile', 'none', str_or_char);
addParameter(inParse, 'fixedIdx', 1, @isnumeric);
addParameter(inParse, 'upCont', 'none', @iscell);
% filter related
addParameter(inParse, 'removeHotPixels', false, @isnumeric);
addParameter(inParse, 'includeHotPixel', false);
addParameter(inParse, 'chi', 0);
addParameter(inParse, 'winSize', 4, @isnumeric);
% other
addParameter(inParse, 'reverse', 0);
addParameter(inParse, 'checkPlot', 0);

parse(inParse, nFolders, varargin{:});
fileName = inParse.Results.fileName;
transFormFile = inParse.Results.transFormFile;
fixedIdx = inParse.Results.fixedIdx;

nFolders = correct_cell_shape(nFolders);

%% tranformation / filtering
transformedData = containers.Map;
nFiles = {};

% define QDM parameters
pixelsize = 4.68e-6;

% check if filename has '.mat'
fileName = check_suffix(fileName);

% generate reference file name
refFile = [nFolders{fixedIdx}, filesep, fileName];
refFileData = load(refFile);

% get transformations and rframes
if inParse.Results.checkPlot
    inParse.Results.checkPlot
    [nTransForms, nRefFrames] = align_images(nFolders, 0, ...
        'fileName', fileName, 'fixedIdx', fixedIdx);
else
    [nTransForms, nRefFrames] = get_tform_multi(refFile, nFolders, ...
    'transFormFile', transFormFile, 'reverse', inParse.Results.reverse);
end

if contains(fileName, 'B111')
    refData = refFileData.B111ferro;
    refLed = refFileData.ledImg;
else
    refData = refFileData.Bz;
    refLed = refFileData.newLED;
end

% cycle through all folders
for i = 1:size(nFolders, 2)
    % create filename
    iFolder = nFolders{i};  
    iFile = fullfile(iFolder, filesep, fileName);
    
    fprintf('<> loading << %s >> target file for transformation\n', iFile(end-50:end))
    
    nFiles{end+1} = iFile;
    
    target = load(iFile);

    if contains(fileName, 'B111')
        targetData = target.B111ferro;
        targetLed = target.ledImg;
    else
        targetData = target.Bz;
        targetLed = target.newLED;
    end
    
    % pre filtering
    targetData = filter_hot_pixels(targetData);

    %% upward cont.
    if iscell(inParse.Results.upCont)
        h = inParse.Results.upCont{i};
        if h ~= 0
            disp(['<>   calculating upward continuation (' num2str(h) ') micron'])
            targetData = UpCont(targetData, h*1e-6, 1/pixelsize);
        end
    end

    %% filtering
    if inParse.Results.removeHotPixels
        if inParse.Results.chi
            chi = target.chi2Pos1 + target.chi2Pos2 + target.chi2Neg1 + target.chi2Neg2;
        else
            chi = inParse.Results.chi;
        end
        disp(['<>   filtering: ...' iFile(end-40:end-20)  '... .mat'])

        targetData = filter_hot_pixels(targetData, ...
            'cutOff',inParse.Results.removeHotPixels, ...
            'includeHotPixel',inParse.Results.includeHotPixel, ...
            'winSize', inParse.Results.winSize, ...
            'checkPlot', inParse.Results.checkPlot,'chi', chi);
    end

    iTransForm = nTransForms(iFile);
    iRefFrame = nRefFrames(iFile);

    %% reverse
    % in the case of reverse, tform and rframe are the ref -> target
    % tranformation. Therefore, the data/led of the target does not need to
    % be transformed. However, the mask itself needs to be transformed form
    % the reference coordinates to the target coordinates later.

    if inParse.Results.reverse
        disp(['<>   WARNING: reverse doesnt work, yet!'])
        transData = targetData;
        transLed = targetLed;
    else
        fprintf('<>   transforming: target data & LED  << ... %s >>\n', iFile(end-50:end))
        transData = tform_data(targetData, iTransForm, iRefFrame);
        transLed = tform_data(targetLed, iTransForm, iRefFrame);
    end

    % create struct for the the transformed data of this file
    fileTransForm = struct;
    
    fileTransForm.refFile = refFile;
    fileTransForm.refData = refData;
    fileTransForm.refLed = refLed;
    
    fileTransForm.fileName = iFile;
    fileTransForm.targetLed = targetLed;
    fileTransForm.targetData = targetData;

    fileTransForm.transData = transData;
    fileTransForm.transLed = transLed;
    fileTransForm.transForm = iTransForm;
    fileTransForm.refFrame = iRefFrame;
    fileTransForm.reverse = inParse.Results.reverse;
    fileTransForm.binning = detect_binning(target);

    % save the result in the trans_data container for later use
    transformedData(iFile) = fileTransForm;
    
    %%% create checkplots
    if inParse.Results.checkPlot
        if i ~= fixedIdx
            check_plot(fileTransForm);
        end
    end
end
end

function check_plot(fileTransForm)
    f = figure('units','normalized','outerposition',[0.2 0.4 0.5 0.5],'NumberTitle', 'off', 'Name', fileTransForm.fileName);
    ax1 = subplot(2,3,1);
    ref = fileTransForm.refData;
    ref = filter_hot_pixels(ref);
    pcolor(ref);
    axis xy; axis equal; axis tight; shading flat;
    title('reference Data')

    ax2 = subplot(2,3,2);
    target = fileTransForm.targetData;
    target = filter_hot_pixels(target);
    pcolor(target);
    axis xy; axis equal; axis tight;shading flat;
    title('target Data')

    ax3 = subplot(2,3,3);
    trans = fileTransForm.transData;
    trans = filter_hot_pixels(trans);
    pcolor(trans);
    axis xy; axis equal; axis tight;shading flat;
    title('transformed Data')

    ax4 = subplot(2,3,4);
    pcolor(fileTransForm.refLed);
    axis xy; axis equal; axis tight; shading flat;
    colormap(ax4, bone)
    title('reference LED')

    ax5 = subplot(2,3,5);
    pcolor(fileTransForm.targetLed);
    axis xy; axis equal; axis tight;shading flat;
    colormap(ax5, bone)
    title('target LED')

    ax6 = subplot(2,3,6);
    pcolor(fileTransForm.transLed);

    axis xy; axis equal; axis tight;shading flat;
    title('transformed LED')
    colormap(ax6, bone)

    linkaxes([ax1 ax2 ax3]);
    linkaxes([ax4 ax5 ax6]);
end