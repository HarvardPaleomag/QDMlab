function results = dipole_fit_series(nFiles, kwargs, filter)
%[results] = dipole_fit_series(nFiles; 'refIdx', 'pixelSize', 'upCont', 'nROI', 'checkPlot', 'imageFolder', 'transFormFile', 'outputTrue', 'save', 'filterProps')
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
%     kwargs.fileName: str
%         default: 'Bz_uc0.mat'
%         name of the .mat file to be used.
%     transFormFile: str
%         default: 'none'
%         absolute path of a file where previously calculated transformations
%         are in. If file does not exist the file will be created so that
%         you dont have to do this every time.
%     kwargs.upCont: cell
%         default: false
%         A cell of numbers for upward continuation
%     kwargs.refIdx: int
%         default: 1
%         index of the reference LED image. This will be used as reference
%         for the transformation calculations
%     reverse: bool -> NOT SUPPORTED, YET
%         default: false
%         if true:  refernce - tform -> target
%         if false: target   - tform -> reference
%     kwargs.outputTrue: bool
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

arguments
    nFiles
    kwargs.refIdx = 1;
    kwargs.pixelSize = 4.68e-6;
    kwargs.upCont = {0};
    kwargs.nROI = false;

    kwargs.checkPlot = false;
    
    kwargs.imageFolder = false;
    kwargs.transFormFile = false;
    kwargs.outputTrue =  false;
    kwargs.save = false;

    filter.filterProps struct = struct();
end

%% filtering
if ~all( structfun(@isempty, filter.filterProps))
    filterProps = namedargs2cell(filter.filterProps);
    filtered = true;
else
    filtered = false;
end

% reference file name
refFile = nFiles{kwargs.refIdx};

% get transformations and rframes
[nTransForms, nRefFrames] = get_tform_multi(refFile, nFolders, ...
                            'transFormFile', kwargs.transFormFile, 'check', kwargs.checkPlot);

% pick source on data
fixed = load(refFile);

[bool, dataName, ledName] = is_B111(fixed);

if bool
    msg = sprintf('B111 data NOT SUPPORTED, YET');
    logMsg('error',msg,1,0);
end

% read data and threshold to 5
fixedData = fixed.(dataName);

if filtered 
    if isfield(filterProps, 'chi')
        filterProps.chi = fixedData.chi2Pos1 + fixedData.chi2Pos2 + fixedData.chi2Neg1 + fixedData.chi2Neg2;
    end
    fixedData = filter_hot_pixels(fixedData, filterProps{:});
end



% read LED
fixedLed = fixed.(ledName);

if ~islogical(kwargs.nROI)
    nROI = kwargs.nROI;
else
    [~, nROI] = pick_box(fixedData, 'led', false, 'closeFig', true);
end

allResults = containers.Map;

numberoffolders=size(nFolders,1);

% preallocate the cells:
iFiles = {};
moments = zeros(size(nROI,2),numberoffolders,size({0,10},2));
inclinations = zeros(size(nROI,2),numberoffolders,size({0,10},2));
declinations = zeros(size(nROI,2),numberoffolders,size({0,10},2));
heights = zeros(size(nROI,2),numberoffolders,size({0,10},2));
residuals = zeros(size(nROI,2),numberoffolders,size({0,10},2));
xMin = zeros(size(nROI,2),numberoffolders,size({0,10},2));
xMax = zeros(size(nROI,2),numberoffolders,size({0,10},2));
yMin = zeros(size(nROI,2),numberoffolders,size({0,10},2));
yMax = zeros(size(nROI,2),numberoffolders,size({0,10},2));

%%
% cycle through all folders
for j = 1 : numberoffolders
    if contains(kwargs.fileName,'.mat')
        iFile = fullfile(nFolders{j}, filesep, kwargs.fileName);
    else
        iFile = fullfile(nFolders{j}, filesep, [kwargs.fileName '.mat']);
    end
    iData = load(iFile);
    
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
    if filtered 
        if isfield(filterProps, 'chi')
            filterProps.chi = movedData.chi2Pos1 + movedData.chi2Pos2 + movedData.chi2Neg1 + movedData.chi2Neg2;
        end
        fixedData = filter_hot_pixels(fixedData, filterProps{:});
    end

    % replace data/LED with transformed data/LED
    iData.newLED = movedLed;
    iData.Bz = movedData;

    % copy the data so that Bz is not overwritten when doing the
    % UC calculations later
    transDataUC = iData;

    fileResults = {};

    % cycle through all rectangles (i.e. sources)
    for i = 1:size(nROI, 2)

        for k = 1:size(kwargs.upCont, 2)
            h = kwargs.upCont{k};

            if h > 0
                % calculate the UC
                msg = sprintf('calculating upward continuation (%i) micron', h);
                logMsg('info',msg,1,0);
                
                % replace the last Bz data with UC data of non UC iData
                transDataUC.Bz = UpCont(iData.Bz, h*1e-6, 1/pixelSize);
            end

            iRect = nROI{i};

            xLim = round([iRect(1), iRect(1) + iRect(3)]);
            yLim = round([iRect(2), iRect(2) + iRect(4)]);
            
            %% actual fitting
            SOURCENAME=['Source' num2str(i) '_Step' num2str(j) ];
            
            % Dipole... returns a struct('dfile', 'm', 'inc', 'dec', 'h', 'res');
            iResult = dipole_fit('filePath', iFile, 'fitOrder', 1, ...
                'cropFactor', 20, 'save', kwargs.save, ...
                'xy', iRect(1:2), 'dx', iRect(3), 'dy', iRect(4), ...
                'expData', transDataUC, 'checkPlot', kwargs.checkPlot, ...
                'imageFolder',kwargs.imageFolder,'sourceName', SOURCENAME);
            
            %% results
            iResult.xLims = xLim;
            iResult.yLims = yLim;
            iResult.iSource = i;
            fileResults{end+1} = iResult;

            iFiles{end+1} = iFile;
            moments(i,j,k) = iResult.m;
            inclinations(i,j,k) = iResult.inc;
            declinations(i,j,k) = iResult.dec;
            heights(i,j,k) = iResult.h;
            residuals(i,j,k) = iResult.res;
            xMin(i,j,k) = xLim(1);
            xMax(i,j,k) = xLim(2);
            yMin(i,j,k) = yLim(1);
            yMax(i,j,k) = yLim(2);
            
%             close all
        end
    end
    allResults(iFile) = fileResults;
end


results = struct('nFiles', {iFiles}, 'moments', moments, 'incs', inclinations, ...
                 'decs', declinations, 'heights', heights, 'res', residuals, ...
                 'xMin', xMin, 'xMax', xMax, 'yMin', yMin, 'yMax', yMax);

