function [nTransForms, nRefFrames] = align_images(nFolders, kwargs)
% Function to aling a set of images. Function will automatically align the
% images first and you can check if it is ok. If not the complex alignment
% will be called.
%
% Parameters
% ----------
%   nFolders: cell, char
%
%   transFormFile: path
%   fixedIdx: int [1]
%   checkPlot: bool [1]
%   fileName: str ['Bz_uc0.mat']
%   sequence: bool [0]
%   reverse: bool [0]
%   laser: uses laser image for alignment

arguments
    nFolders cell {foldersMustExist(nFolders)}
    kwargs.transFormFile char = 'none';
    kwargs.fixedIdx int16 = 1
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)} = 1
    kwargs.fileName char {mustBeMember(kwargs.fileName, ['Bz_uc0', 'B111dataToPlot', 'Bz_uc0.mat', 'B111dataToPlot.mat']), ...
                          fileMustExistInFolder(kwargs.fileName, nFolders)} = 'Bz_uc0'
    kwargs.sequence (1,1) {mustBeBoolean(kwargs.sequence)} = 0
    kwargs.reverse (1,1) {mustBeBoolean(kwargs.reverse)} = 0
    kwargs.laser (1,1) {mustBeBoolean(kwargs.laser)} = false
end
nFolders = correct_cell_shape(nFolders);

fileName = check_suffix(kwargs.fileName);

fixedIdx = kwargs.fixedIdx;
sequence = kwargs.sequence;
reverse = kwargs.reverse;

% generate reference file name
fixedFile = [nFolders{fixedIdx}, filesep, fileName];
fixedData = load(fixedFile);


if contains(fileName, 'B111')
    fixedLed = fixedData.ledImg;
else
    fixedLed = fixedData.newLED;
end

if kwargs.laser
    fixedLed = get_laser(nFolders{fixedIdx}, 'data', fixedData);
end
% resize led to match data binning
% todo see if it is not better to not do that but in the functiopn that
% calculates the transformation
% fixed_led = imwarp(fixed_led, resize_binning);
nTransForms = containers.Map;
nRefFrames = containers.Map;

% check if file exists an/or should be created
if kwargs.transFormFile ~= 0 | ~strcmp(kwargs.transFormFile, 'none')
    if isfile(kwargs.transFormFile)
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

for iFolder = nFolders

    iFolder = iFolder{:};

    % load moving data
    movingPath = fullfile(iFolder, filesep, fileName);
    movingData = load(movingPath);

    % check for differences in LED naming
    if kwargs.laser
        movingLed = get_laser(movingPath, 'data', movingData);
    elseif isfield(movingData, 'ledImg')
        movingLed = movingData.ledImg;
    elseif isfield(movingData, 'newLED')
        movingLed = movingData.newLED;
    end
    
    msg = sprintf('aligning images:');
    logMsg('info',msg,1,0);
    if reverse
        msg = sprintf('reversed alignment (fixed -> moving)');
        logMsg('info',msg,1,0);
%         disp(['<>   INFO: reversed alignment (fixed -> moving)'])
        logMsg('info',[fixedFile, '->'],1,1);
        logMsg('info',[movingPath],1,1);
        [tForm, refFrame] = get_image_tform(movingLed, fixedLed, 'check', true);
    else
        logMsg('info',[movingPath, '->'],1,1);
        logMsg('info',[fixedFile],1,1);
        [tForm, refFrame] = get_image_tform(fixedLed, movingLed, 'check', true);
    end
    % ask only if the images are different

    if ~strcmp(movingPath,fixedFile)
        ok = questdlg('Is the allignment ok?', 'Yes', 'No');
        switch ok
            case 'Yes'
                close
            case 'No'
                [tForm, refFrame] = get_image_tform_complex(fixedLed, movingLed, 'check', true);
                pause(5)
                close
        end
    end

    nTransForms([movingPath]) = tForm;
    nRefFrames([movingPath]) = refFrame;

    if sequence
        fixedLed = movingLed;
        fixedFile = movingPath;
    end

end

if kwargs.transFormFile == 0 | strcmp(kwargs.transFormFile, 'none')
    msg = sprintf('returning: Map(tForms), Map(refFrames)');
    logMsg('info',msg,1,0);
%     disp('<>   returning: Map(tForms), Map(refFrames)')
    return
else
    msg = sprintf('saving... to << %s >>', transFormFile');
    logMsg('info',msg,1,0);
%     fprintf('<>   saving... to ''%s''', transFormFile)
    save(transFormFile, 'nTransForms', 'nRefFrames', 'reverse')
end
