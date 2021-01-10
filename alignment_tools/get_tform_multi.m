function [nTransForms, nRefFrames] = get_tform_multi(fixedFile, nMovingFolders, varargin)
% parameters:
%     fixedFile: str
%         Path to the reference image. Needs to be a file not folder
%     nMovingFolders: cell(str)
%         A cell with all folders that contain the data to be aligned.
% 
% optional parameters:
%         tform_file: str
%             default: 'none'
%             Location for the file the transformation map is stored in for
%             later use
%         binning: int
%             default: 2
%             binning used in fit
%         reverse: bool
%             if true:  refernce - tform -> target
%             if false: target   - tform -> reference
%         checkPlot: bool
%             default: false
%             If true creates an overlay image of reference/target for each
%             file in target_fodlers matching the filename of reference_file
% 
% Returns:
%     a map with the filname as the key and the transform as the value. Can
%     be accessed like `nTransForms(fname) = iTransForm`

inParams = inputParser;
str_or_char = @(x) isstring(x) | ischar(x);

addRequired(inParams, 'fixedFile', str_or_char);
addRequired(inParams, 'nMovingFolders', @iscell);
addParameter(inParams, 'transFormFile', false);
addParameter(inParams, 'checkPlot', false, @islogical);
addParameter(inParams, 'reverse', false, @islogical);
addParameter(inParams, 'binning', 2, @islogical);
addParameter(inParams, 'laser', 0, @islogical);

parse(inParams, fixedFile, nMovingFolders, varargin{:});

transFormFile = inParams.Results.transFormFile;
checkPlot = inParams.Results.checkPlot;
reverse = inParams.Results.reverse;
binning = inParams.Results.binning;
laser = inParams.Results.laser;

nMovingFolders = correct_cell_shape(nMovingFolders);

% define resizing affine transformation
LED2data = affine2d([1 / binning, 0, 0; 0, 1 / binning, 0; 0, 0, 1]);
data2LED = affine2d([binning, 0, 0; 0, binning, 0; 0, 0, 1]);

[fixedPath, refFileName, refExtension] = fileparts(fixedFile);

if laser
    refExtension = 'laser.jpg';
    fixedData = imread([fixedPath filesep refExtension]);
else
    fixedData = load(fixedFile);
end


if laser
    fixedLed = fixedData;
%     fixedLed = adapthisteq(fixedLed);
%     fixedLed = (fixedLed / max(fixedLed, [], 'all')) * 512;
    
elseif contains(refFileName, 'B111')
    fixedLed = fixedData.ledImg;
else
    fixedLed = fixedData.newLED;
end

% resize led to match data binning
% todo see if it is not better to not do that but in the functiopn that
% calculates the transformation
nTransForms = containers.Map;
nRefFrames = containers.Map;

% check if file exists an/or should be created
if transFormFile ~= 0
    if isfile(transFormFile)
        newFile = input('transformation file already exists, overwrite? (y/n)? ', 's');
        if strcmp(newFile, 'y')
            nTransForms = containers.Map;
            nRefFrames = containers.Map;
        else
            load(transFormFile);
            return
        end
    end
end

for iFolder = nMovingFolders

    iFolder = iFolder{:};

    % load moving data
    movingPath = [iFolder, filesep, refFileName];

    if laser
        moving = imread([iFolder filesep refExtension]);
%         moving = adapthisteq(moving);
    else
        moving = load(movingPath);
    end
    % check for differences in LED naming
    if isfield(moving, 'ledImg')
        movingLed = moving.ledImg;
    elseif isfield(moving, 'newLED')
        movingLed = moving.newLED;
    elseif laser
        movingLed = moving;
        movingLed = movingLed - min(movingLed, [], 'all');
    end
    
    if reverse
        [tForm, refFrame] = get_image_tform2(movingLed, fixedLed, 'checkPlot', checkPlot);
        disp(['<>   ', fixedFile, '->'])
        disp(['<>   ', movingPath])
    else
        [tForm, refFrame] = get_image_tform2(fixedLed, movingLed, 'checkPlot', checkPlot);
        disp(['<>   ', movingPath, '->'])
        disp(['<>   ', fixedFile])
    end
    
    nTransForms([movingPath, '.mat']) = tForm;
    nRefFrames([movingPath, '.mat']) = refFrame;

end

if transFormFile == 0
    disp('<>   returning: Map(tForms), Map(refFrames)')
    return
else
    disp(['<>   saving... to ', transFormFile])
    save(transFormFile, 'nTransForms', 'nRefFrames')
end
