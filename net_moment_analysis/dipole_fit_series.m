function results = dipole_fit_series(nFolders, kwargs, constrains, filter)
%[results] = dipole_fit_series(nFolders; 'transFormFile', 'refIdx', 'checkPlot', 'outputTrue', 'save', 'upCont', 'nROI', 'imageFolder', 'fileName')
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
%     nFolders: list
%         list of absolute path to each data folder. First entry in the list
%         is used as reference image.
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
%     constraints
%     -----------
%     m0: double [1e-12]
%       initial moment guess
%     hguess: double [2.5e-5]
%       initial height guess
%     minheight: double [0]
%       minimum height constraint for fit
%     maxheight: double [100e-6]
%       maximum height constraint for fit
%     boxwidth: double [100e-6]
%       x/y constraint for fit
%     filterProps: struct [struct()]
%       passed to filter_hot_pixels
%
% Returns
% -------
%     struct with cells for:
%             - nFiles{i,j,k}
%             - moments{i,j,k}
%             - incs{i,j,k}
%             - decs{i,j,k}
%             - heights{i,j,k}
%             - dipolarity{i,j,k}
%             - x{i,j,k}
%             - y{i,j,k}
%             - data{i,j,k}
%             - model{i,j,k}
%             - res{i,j,k}
%             - xMin{i,j,k}
%             - xMax{i,j,k}
%             - yMin{i,j,k}
%             - yMax{i,j,k}
%             - individualResults{i,j,k} results for each ROI/file/UC 
%
%     Each is a cell with result{i,j,k} with:
%           ith ROI, jth file, kth UC value:

arguments
    nFolders
    kwargs.refIdx = 1;
    kwargs.upCont = {0};
    kwargs.nROI = false;
    kwargs.pixelSize = 4.68e-6;
    kwargs.nRuns = 10;
    
    kwargs.outputTrue = false;
    kwargs.save = false;

    kwargs.checkPlot = false;
    kwargs.transFormFile = false;

    kwargs.imageFolder = false;
    kwargs.fileName = 'Bz_uc0.mat';

    constrains.constrained {mustBeBoolean(constrains.constrained)} = false;
    constrains.m0 = 1e-12;
    constrains.hguess = 2.5e-5;
    constrains.minheight = 0;
    constrains.maxheight = 100e-6;
    constrains.boxwidth = 100e-6;

    filter.filterProps struct = struct();
end

if ~all(structfun(@isempty, filter.filterProps))
    filterProps = namedargs2cell(filter.filterProps);
    filtered = true;
else
    filtered = false;
end

nFolders = correct_cell_shape(nFolders);
% generate reference file name
refFile = [nFolders{kwargs.refIdx}, filesep, kwargs.fileName];

% get transformations and rframes
[nTransForms, nRefFrames] = get_tform_multi(refFile, nFolders, ...
    'transFormFile', kwargs.transFormFile, 'check', kwargs.checkPlot);

% pick source on data
fixed = load(refFile);

[bool, dataName, ledName] = is_B111(fixed);

if bool
    msg = sprintf('B111 data NOT SUPPORTED, YET');
    logMsg('error', msg, 1, 0);
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

numberoffolders = size(nFolders, 2);

% preallocate the cells:
preAllocatedArray = zeros(size(nROI, 2), numberoffolders, size(kwargs.upCont, 2));
iFiles = cell(size(preAllocatedArray));
moments = preAllocatedArray;
inclinations = preAllocatedArray;
declinations = preAllocatedArray;
heights = preAllocatedArray;
datas = cell(size(preAllocatedArray));
models = cell(size(preAllocatedArray));
residuals = preAllocatedArray;
dipolarity = preAllocatedArray;
xMin = preAllocatedArray;
xMax = preAllocatedArray;
yMin = preAllocatedArray;
yMax = preAllocatedArray;
xloc = preAllocatedArray;
yloc = preAllocatedArray;
fileResults = num2cell(preAllocatedArray);

%%
% cycle through all folders
for j = 1:numberoffolders
    if contains(kwargs.fileName, '.mat')
        iFile = fullfile(nFolders{j}, filesep, kwargs.fileName);
    else
        iFile = fullfile(nFolders{j}, filesep, [kwargs.fileName, '.mat']);
    end
    iData = load(iFile);

    % get transformation for that file
    tForm = nTransForms(iFile);
    % get refFrame for that file
    rframe = nRefFrames(iFile);

    % claculate transformed data
    msg = sprintf('transforming LED');
    logMsg('info', msg, 1, 0);
    movedLed = tform_data(iData.newLED, tForm, rframe);
    msg = sprintf('transforming Bz');
    logMsg('info', msg, 1, 0);
    movedData = tform_data(iData.Bz, tForm, rframe);

    %% filtering
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
                logMsg('info', msg, 1, 0);

                % replace the last Bz data with UC data of non UC iData
                transDataUC.Bz = UpCont(iData.Bz, h*1e-6, 1/kwargs.pixelSize);
            end

            iRect = nROI{i};

            xLim = round([iRect(1), iRect(1) + iRect(3)]);
            yLim = round([iRect(2), iRect(2) + iRect(4)]);

            %% actual fitting
            SOURCENAME = ['Source', num2str(i), '_Step', num2str(j)];

            % Dipole... returns a struct('dfile', 'm', 'inc', 'dec', 'h', 'res');
            constrainArgs = namedargs2cell(constrains);
            iResult = dipole_fit('filePath', iFile, 'fitOrder', 1, 'nRuns', kwargs.nRuns, ...
                'cropFactor', 20, 'save', kwargs.save, ...
                'xy', iRect(1:2), 'dx', iRect(3), 'dy', iRect(4), ...
                'expData', transDataUC, 'checkPlot', kwargs.checkPlot, ...
                constrainArgs{:}, ...
                'imageFolder', kwargs.imageFolder, 'sourceName', SOURCENAME);

            %% results
            iResult.xLims = xLim;
            iResult.yLims = yLim;
            iResult.iSource = i;
            iResult.sourceName = SOURCENAME;
            
            fileResults{i,j,k} = iResult;
            iFiles{i, j, k} = iFile;
            % fit results
            moments(i, j, k) = iResult.m;
            inclinations(i, j, k) = iResult.inc;
            declinations(i, j, k) = iResult.dec;
            heights(i, j, k) = iResult.h;
            dipolarity(i, j, k) = iResult.dipolarity;
            xloc(i, j, k) = iResult.x;
            yloc(i, j, k) = iResult.y;
            
            % data
            datas{i, j, k} = iResult.data;
            models{i, j, k} = iResult.model;
            residuals(i, j, k) = iResult.res;
            
            % croplimits
            xMin(i, j, k) = xLim(1);
            xMax(i, j, k) = xLim(2);
            yMin(i, j, k) = yLim(1);
            yMax(i, j, k) = yLim(2);
            %             close all
        end
    end
    
end


results = struct();
results.nFiles = iFiles;

results.moments = moments;
results.decs = declinations;
results.incs = inclinations;                
results.heights = heights;
results.dipolarity = dipolarity;
results.x = xloc;
results.y = yloc;

results.data = datas;
results.model = models;
results.res = residuals;
                 
results.xMin = xMin;
results.xMax = xMax; 
results.yMin = yMin; 
results.yMax = yMax;

results.individualResults = fileResults;
msg = sprintf('returning structure where result{i,j,k} is: ith ROI, jth file, kth UC value');
logMsg('result',msg,1,0);
