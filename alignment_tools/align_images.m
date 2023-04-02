function [nTransForms, nRefFrames] = align_images(nFolders, kwargs)
% Align a set of images to a fixed reference image.
%
% INPUTS:
% - nFolders: A cell array of folder paths containing the image files to
%   be aligned. Each folder should contain one or more image files with the
%   same prefix and file type (e.g., "img_*.png").
%
% - kwargs: A struct of optional parameters for the alignment. The following
%   fields are recognized:
%   - fixedIdx: The index of the fixed reference image in the folderPaths
%     array. Default is 1.
%   - checkPlot: A boolean value indicating whether to show a plot to check
%     the alignment of each image. Default is true.
%   - reference: A string indicating the type of reference image to use for
%     the alignment. Possible values are "led" (default) or "laser".
%   - laser: A boolean value indicating whether to use a laser image for
%     alignment instead of the LED image. Default is false.
%
% OUTPUTS:
% - tForms: A cell array of transformation matrices that align each image to
%   the fixed reference image. Each matrix is a 3x3 affine transformation
%   matrix that maps the coordinates of the input image to the coordinates
%   of the reference image.
%
% - refFrames: A cell array of reference frames for each image. Each frame
%   is a 2D matrix that represents the intensity values of the reference
%   image at the transformed coordinates of the input image.
%
% EXAMPLE USAGE:
% folderPaths = {'/path/to/folder1', '/path/to/folder2', ...};
% options = struct('fixedIdx', 1, 'checkPlot', true, 'reference', 'led', 'laser', false);
% [tForms, refFrames] = align_images(folderPaths, options);
%
% NOTE:
% - This function assumes that all input images have the same size and
%   aspect ratio.
% - The alignment algorithm uses a simple cross-correlation method to
%   estimate the initial transformation and a manual check to refine the
%   alignment if necessary. For best results, the input images should have
%   high contrast and low noise.
% - The function currently only supports grayscale images. To align
%   multi-channel or color images, you need to align each channel separately
%   and then combine the aligned channels into a single image.
% - The function does not currently support aligning images with different
%   resolutions or orientations. To align such images, you need to apply a
%   suitable transformation to the input images before aligning them.
% See also
% --------
%  'get_image_transform', 'get_image_tform_complex'

arguments
    nFolders cell {foldersMustExist(nFolders)}
    kwargs.transFormFile char = 'none';
    kwargs.fixedIdx int16 = 1
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = true;
    kwargs.fileName char {fileMustExistInFolder(kwargs.fileName, nFolders)} = 'Bz_uc0'
    kwargs.sequence (1,1) {mustBeBoolean(kwargs.sequence)} = 0
    kwargs.reverse (1,1) {mustBeBoolean(kwargs.reverse)} = 0
    kwargs.reference = 'led'
    kwargs.laser (1,1) {mustBeBoolean(kwargs.laser)} = false
end
nFolders = correct_cell_shape(nFolders);

fileName = check_suffix(kwargs.fileName);

fixedIdx = kwargs.fixedIdx;
sequence = kwargs.sequence;
reverse = kwargs.reverse;

% generate reference file name
fixedFile = fullfile(nFolders{fixedIdx}, filesep, fileName);
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
% todo see if it is not better to not do that but in the function that
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
