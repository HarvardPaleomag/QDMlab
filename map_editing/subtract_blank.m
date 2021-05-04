function subtractedData = subtract_blank(nFolders, blankFolder, kwargs)
% function subtractedData = subtract_blank(kwargs)
% Subtracts a blank map from the Data
%
% Parameters
% ----------
%   nFolders ['none']
%       A path (char) or cell of paths that contain the data and are to be
%       corrected.
%   blankFolder: path ['none']
%       The path of the blank B111 map.
%   checkPlot: bool [false]
%       if true: creates a plot to check if the subtraction worked as
%       expected
%   save: bool [true]
%       if true: a new file with the name: `B111BlankSub.mat` will be
%       created. Otherwise the results will only be returned.
%
% Note
% ----
%   If no :code:`nFolders` or :code:`blankFolder` is passed to the function,
%   you will be prompted to select them.

arguments
    nFolders
    blankFolder
    kwargs.checkPlot {mustBeBoolean(kwargs.checkPlot)} = false
    kwargs.save {mustBeBoolean(kwargs.save)} = true
end

close all

checkPlot = kwargs.checkPlot;

fileName = 'B111dataToPlot.mat';
laserFileName = 'laser.jpg';

%% manual
% if nFoilders and blankData uses default values i.e. false

if strcmp(nFolders, 'none')
    %Load the correct "B111dataToPlot.mat" file
    f = helpdlg('Pick a B111 file');
    pause(2)

    if ishghandle(f)
        close(f)
    end

    [longfilename, pathname] = uigetfile('*.mat', 'Pick a B111 file');
    fullfilename = [pathname longfilename];
    nFolders = {pathname};
end

if strcmp(blankFolder, 'none')
    f = helpdlg('Pick a blank file');
    pause(2)

    if ishghandle(f)
        close(f)
    end

    [longfilenameBLANK, pathnameBLANK] = uigetfile('*.mat', 'Pick a Blank file');
    fullfilenameBLANK=[pathnameBLANK longfilenameBLANK];
    blankFolder = fullfilenameBLANK;
end

%% automatic subtraction for all folders
% checks if none of the default arguments is used
nFolders = correct_cell_shape(nFolders);

disp(['<> loading blank file: <<' blankFolder '>>'])
blankFile = fullfile(blankFolder, fileName);
blankData = load(blankFile);

% [nTransForms, nRefFrames] = get_tform_multi(blankFolder, nFolders, ...
%                             'reverse', true, ...
%                             'laser',true, 'checkPlot', checkPlot);

movingData = imread(fullfile(blankFolder, laserFileName));
subtractedData = containers.Map();

for i = 1 : size(nFolders, 2)
    iFolder = nFolders{i};
    iFile = fullfile(iFolder, filesep, fileName);

    disp(['<> reading: << ' iFile ' >> and blankData'])
    fixedData = imread(fullfile(iFolder, laserFileName));
    fileData = load(iFile);
    
    [transForm, refFrame] = get_image_tform(fixedData, movingData,...
        'checkPlot', kwargs.checkPlot, 'title', 'laser alignment');

    % transform the blank
    binning = detect_binning(fileData);
    B111ferroTransformed = tform_data(blankData.B111ferro, transForm, refFrame, binning);
    B111paraTransformed = tform_data(blankData.B111para, transForm, refFrame, binning);

    % crop the FOV and subtract blank
    [x, y, w, h] = get_mask_extent(B111ferroTransformed);
    fileB111ferro = fileData.B111ferro(y:y+h, x:x+w);
    fileB111para = fileData.B111para(y:y+h, x:x+w);
    
    fileData.B111ferro = fileB111ferro- B111ferroTransformed(y:y+h, x:x+w);
    fileData.B111para = fileB111para - B111paraTransformed(y:y+h, x:x+w);
%     B111ferro = fileData.B111ferro - B111ferroTransformed;
%     B111para = fileData.B111para - B111paraTransformed;
    subtractedData(iFolder) = fileData;
    
    if kwargs.save
        saveFilePath = fullfile(iFolder, 'B111BlankSub.mat');
        save(saveFilePath, '-struct', 'fileData');
    end

    if kwargs.checkPlot
        figure
        sp1 = subplot(2,2,1);
        imagesc(fileData.B111ferro);
        title('Original');
        sp2 = subplot(2,2,2);
        imagesc(blankData.B111ferro);
        title('Blank');
        sp3 = subplot(2,2,3);
        blank = re_bin(blankData.B111ferro, fileData.B111ferro);
        imagesc(fileData.B111ferro - blank);
        title('Unaligned subtraction');
        sp4 = subplot(2,2,4);
        imagesc(B111ferro);
        title('Aligned subtraction');
        linkaxes([sp1, sp2, sp3, sp4])
    end
end
