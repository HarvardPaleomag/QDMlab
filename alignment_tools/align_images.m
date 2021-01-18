function [nTransForms, nRefFrames] = align_images(nFolders, transFormFile, kwargs)
% Function to aling a set of images. Function will automatically align the
% images first and you can check if it is ok. If not the complex alignment
% will be called.

arguments
    nFolders cell
    transFormFile
    kwargs.fixedIdx = 1
    kwargs.checkPlot = 1
    kwargs.fileName = 'Bz_uc0.mat'
    kwargs.sequence = 0
    kwargs.reverse = 0
end

nFolders = correct_cell_shape(nFolders);

fileName = kwargs.fileName;

if ~endsWith(fileName, '.mat')
    fileName = [fileName, '.mat'];
end

checkPlot = kwargs.checkPlot;
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

% resize led to match data binning
% todo see if it is not better to not do that but in the functiopn that
% calculates the transformation
% fixed_led = imwarp(fixed_led, resize_binning);
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

for iFolder = nFolders

    iFolder = iFolder{:};

    % load moving data
    movingPath = fullfile(iFolder, filesep, fileName);
    moving = load(movingPath);

    % check for differences in LED naming
    if isfield(moving, 'ledImg')
        movingLed = moving.ledImg;
    elseif isfield(moving, 'newLED')
        movingLed = moving.newLED;
    end
    
    if reverse
        disp(['<>   INFO: reversed alignment (fixed -> moving)'])
        disp(['<>   ', fixedFile, '->'])
        disp(['<>   ', movingPath])
        [tForm, refFrame] = get_image_tform2(movingLed, fixedLed, 'check', true);        
    else
        disp(['<>   ', movingPath, '->'])
        disp(['<>   ', fixedFile])
        [tForm, refFrame] = get_image_tform2(fixedLed, movingLed, 'check', true);
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

if transFormFile == 0
    disp('<>   returning: Map(tForms), Map(refFrames)')
    return
else
    fprintf('<>   saving... to ''%s''', transFormFile)
    save(transFormFile, 'nTransForms', 'nRefFrames', 'reverse')
end