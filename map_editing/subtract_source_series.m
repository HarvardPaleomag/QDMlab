function results = subtract_source_series(nFolders, varargin)
% pick_sources_and_fit is used to bulk analyze datasets it 
% (1) registers the maps with respect to the first file in nFolders
% (2) lets you pick the sources (can be passed using 'nROI' parameter)
% (3) calculates the upwards continuation (UC) for each source (passed through
%     'upCont' parameter
% (4) fits a dipole to each source in each map at every UC distance
%     -> output structure is {source, map, UC}
% 
% Parameters
% ----------
%     positional
%     ----------
%     folders:list
%         list of absolute path to each data folder. First entry in the list
%         is used as reference image.
% 
%     optional keywords
%     -----------------
%     fileName: str
%         default: 'Bz_uc0.mat'
%         name of the .mat file to be used.
%     transFormFile: str
%         default: 'none'
%         absolute path of a file where previously calculated transformations
%         are in. If file does not exist the file will be created so that
%         you dont have to do this every time.
%     upCont: cell
%         default: false
%         A cell of numbers for upward continuation
%     refIdx: int
%         default: 1
%         index of the reference LED image. This will be used as reference
%         for the transformation calculations
%     reverse: bool -> NOT SUPPORTED, YET
%         default: false
%         if true:  refernce - tform -> target
%         if false: target   - tform -> reference
%     outputTrue: bool
%         default: false
%         if true `dipoleinversions.txt` will be written to disk
%         if false `dipoleinversions.txt` will not be written to disk
%     nROI:
%         selections
%     IMAGEFOLDER: char
%         name of folder to store fit images
%     save: bool (false)
%         if true: individual fits are saved
%         if false: fits are not saved
%
% Returns
% -------
%     struct with cells for:
%            - nFiles
%             - moments{i,j,k}
%             - incs{i,j,k}
%             - decs{i,j,k}
%             - heights{i,j,k}
%             - residuals{i,j,k}
%             - xMin{i,j,k}
%             - xMax{i,j,k}
%             - yMin{i,j,k}
%             - yMax{i,j,k}
%     Each is a cell with result{i,j,k} with:
%           ith ROI, jth file, kth UC value:


nParams = inputParser;
str_or_char = @(x) isstring(x) | ischar(x);

addRequired(nParams, 'nFolders', @iscell);
addParameter(nParams, 'fileName', 'Bz_uc0.mat', str_or_char);
addParameter(nParams, 'transFormFile', 'none', str_or_char);
addParameter(nParams, 'refIdx', 1, @isnumeric);
addParameter(nParams, 'checkPlot', false, @islogical);
addParameter(nParams, 'outputTrue', false, @islogical);
addParameter(nParams, 'save', false, @islogical);
addParameter(nParams, 'upCont', {0}, @iscell);
addParameter(nParams, 'nROI', false, @iscell);
addParameter(nParams, 'imageFolder', false, @ischar);
addParameter(nParams, 'unTransform', true, @islogical);
addParameter(nParams, 'saveAllSub', false, @islogical);
parse(nParams, nFolders, varargin{:});

% define optional parameters
fileName = nParams.Results.fileName;
tformFile = nParams.Results.transFormFile;
refIdx = nParams.Results.refIdx;
checkPlot = nParams.Results.checkPlot;
outputTrue = nParams.Results.outputTrue;
upCont = nParams.Results.upCont;
imageFolder=nParams.Results.imageFolder;
unTransform=nParams.Results.unTransform;

% define QDM parameters
pixelsize = 4.68e-6;

% generate reference file name
refFile = [nFolders{refIdx}, filesep, fileName];

% get transformations and rframes
[nTransForms, nRefFrames] = get_tform_multi(refFile, nFolders, ...
    'transFormFile', tformFile, 'check', checkPlot);

% pick source on data
fixed = load(refFile);

if contains(fileName, 'B111')
    msg = sprintf('B111 data NOT SUPPORTED, YET');
    logMsg('error',msg,1,0);
    fixedData = fixed.B111ferro;
    fixedLed = fixed.ledImg;
else
    fixedData = fixed.Bz;
    fixedLed = fixed.newLED;
end

fixedData = filter_hot_pixels(fixedData, 'cutOff', 12);

if ~islogical(nParams.Results.nROI)
    nROI = nParams.Results.nROI;
else
    nROI = pick_box(fixedData, 'led', false, ...
        'returnCoordinates', true, 'closeFig', true);
end

allResults = containers.Map;

numberoffolders=size(nFolders,1);

% preallocate the cells:
iFiles = {};
moments = zeros(size(nROI,2),numberoffolders);
inclinations = zeros(size(nROI,2),numberoffolders);
declinations = zeros(size(nROI,2),numberoffolders);
heights = zeros(size(nROI,2),numberoffolders);
residuals = zeros(size(nROI,2),numberoffolders);
xMin = zeros(size(nROI,2),numberoffolders);
xMax = zeros(size(nROI,2),numberoffolders);
yMin = zeros(size(nROI,2),numberoffolders);
yMax = zeros(size(nROI,2),numberoffolders);

if nParams.Results.saveAllSub
    BzSub = cell(size(nROI,2),numberoffolders);
else
    BzSub = cell(1,numberoffolders);
end

% list to record failed fits
errList = [];


%%
% cycle through all folders
for j = 1 : numberoffolders
    if contains(fileName,'.mat')
        iFile = fullfile(nFolders{j}, filesep, fileName);
    else
        iFile = fullfile(nFolders{j}, filesep, [fileName '.mat']);
    end
    iData = load(iFile);
    binning = round(detect_binning(iData));
    
    % get transformation for that file
    tForm = nTransForms(iFile);
    % get refFrame for that file
    rframe = nRefFrames(iFile);

    % claculate transformed data
    msg = sprintf('transforming LED');
    logMsg('info',msg,1,0);
    movedLed = tform_data(iData.newLED, tForm, rframe);
    msg = sprintf('transforming Bz');
    logMsg('info',msg,1,0);
    movedData = tform_data(iData.Bz, tForm, rframe);
    
    % todo add filter_hot_pixel
    % movedData = filter_hot_pixels(movedData, 'cutOff', 12);
    
    % replace data/LED with transformed data/LED
    iData.newLED = movedLed;
    iData.Bz = movedData;

    % copy the data so that Bz is not overwritten when doing the
    % UC calculations later
    transDataUC = iData;

    fileResults = {};
    disp(j)
    
    
    % cycle through all rectangles (i.e. sources)
    for i = 1:size(nROI, 2)

       
        iRect = nROI{i};

        xLim = round([iRect(1), iRect(1) + iRect(3)]);
        yLim = round([iRect(2), iRect(2) + iRect(4)]);

        %% actual fitting
        SOURCENAME=['Source' num2str(i) '_Step' num2str(j) ];
        try
            % Dipole... returns a struct('dfile', 'm', 'inc', 'dec', 'h', 'res');
            iResult = dipole_fit('filePath', iFile, 'fitOrder', 1, ...
                'cropFactor', 20, 'save', nParams.Results.save, ...
                'xy', iRect(1:2), 'dx', iRect(3), 'dy', iRect(4), ...
                'expData', transDataUC, 'method', 0, ...
                'imageFolder',imageFolder,'sourceName',SOURCENAME);

            % reassign subtracted region
            ic = 1;
            for is = xLim(1):xLim(2)
                jc = 1;
                for js = yLim(1):yLim(2)
                    transDataUC.Bz(js, is) = -iResult.residuals(jc, ic);%-residualMap(jc, ic);
                    jc = jc + 1;
                end
                ic = ic + 1;
            end

            %% results
            iResult.xLims = xLim;
            iResult.yLims = yLim;
            iResult.iSource = i;
            fileResults{end+1} = iResult;

            iFiles{end+1} = iFile;
            moments(i,j) = iResult.m;
            inclinations(i,j) = iResult.inc;
            declinations(i,j) = iResult.dec;
            heights(i,j) = iResult.h;
            residuals(i,j) = iResult.res;
            xMin(i,j) = xLim(1);
            xMax(i,j) = xLim(2);
            yMin(i,j) = yLim(1);
            yMax(i,j) = yLim(2);
        
            if unTransform && nParams.Results.saveAllSub
                BzSub{i,j} = tform_data(transDataUC.Bz, invert(tForm), rframe);
            elseif nParams.Results.saveAllSub
                BzSub{i,j} = transDataUC.Bz;
            end

            figure();imagesc(transDataUC.Bz);
            axis xy, axis equal, axis tight;
            colorbar;
            colormap(jet);
            caxis([-1,1]*8*10^(-7));
        catch
            errList(end+1,:) = [i,j];
            % just don't want it to stop every time it can't fit something!
        end
        
        close all


    end
    
    
    if unTransform && ~nParams.Results.saveAllSub
        BzSub{1,j} = tform_data(transDataUC.Bz, invert(tForm), rframe);
    elseif ~nParams.Results.saveAllSub
        BzSub{1,j} = transDataUC.Bz;
    end
    
    allResults(iFile) = fileResults;
end


results = struct('nFiles', {iFiles}, 'moments', moments, 'incs', inclinations, ...
                 'decs', declinations, 'heights', heights, 'res', residuals, ...
                 'xMin', xMin, 'xMax', xMax, 'yMin', yMin, 'yMax', yMax, ...
                 'BzSub', {BzSub},'errList',errList);           
	
             
