function [nTransForms, nRefFrames] = get_tform_multi(fixedFile, nMovingFolders, kwargs)
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

arguments
    fixedFile
    nMovingFolders
    kwargs.transFormFile = 'none'
    kwargs.checkPlot  (1,1) {mustBeBoolean(kwargs.checkPlot)} = 0
    kwargs.reverse  (1,1) {mustBeBoolean(kwargs.reverse)} = 0
	kwargs.binning (1,1) {mustBePositive} = 2;
    kwargs.laser  (1,1) {mustBeBoolean(kwargs.laser)} = 0
end

transFormFile = kwargs.transFormFile;
checkPlot = kwargs.checkPlot;
reverse = kwargs.reverse;
binning = kwargs.binning;
laser = kwargs.laser;

nMovingFolders = correct_cell_shape(nMovingFolders);

fixedFile = check_suffix(fixedFile);
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
    movingFile = [iFolder, filesep, refFileName];
    movingFile = check_suffix(movingFile);
    
    if laser
        moving = imread([iFolder filesep refExtension]);
%         moving = adapthisteq(moving);
    else
        moving = load(movingFile);
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
    if movingLed == fixedLed
        
    end
    
    if reverse
        [tForm, refFrame] = get_image_tform(movingLed, fixedLed, 'checkPlot', checkPlot);
        disp(['<>   ', fixedFile, '->'])
        disp(['<>   ', movingFile])
    else
        [tForm, refFrame] = get_image_tform(fixedLed, movingLed, 'checkPlot', checkPlot);
        disp(['<>   ', movingFile, '->'])
        disp(['<>   ', fixedFile])
    end
    
    nTransForms(movingFile) = tForm;
    nRefFrames(movingFile) = refFrame;

end

if transFormFile == 0
    msg = sprintf('returning: Map(tForms), Map(refFrames)');
    logMsg('info',msg,1,0);
    return
else
    msg = sprintf('saving... to ', transFormFile');
    logMsg('info',msg,1,0);
%     disp(['<>   saving... to ', transFormFile])
    save(transFormFile, 'nTransForms', 'nRefFrames')
end
