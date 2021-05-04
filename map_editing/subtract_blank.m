function subtractedData = subtract_blank(kwargs)
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
    kwargs.nFolders = 'none'
    kwargs.blankFolder = 'none'
    kwargs.checkPlot {mustBeBoolean(kwargs.checkPlot)} = false
    kwargs.save {mustBeBoolean(kwargs.save)} = true
end

close all

fileName = 'B111dataToPlot.mat';
laserFileName = 'laser.jpg';

%% manual
% if nFoilders and blankData uses default values i.e. false

nFolders = automatic_input_ui__(kwargs.nFolders, 'title', 'Select measurement folder');
blankFolder = automatic_input_ui__(kwargs.blankFolder, 'title', 'Select blank folder');
blankFolder = blankFolder{:};


%% automatic subtraction for all folders
% checks if none of the default arguments is used
nFolders = correct_cell_shape(nFolders);

msg = sprintf('loading blank file: << %s >>', blankFolder);
logMsg('info',msg,1,0);

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

    msg = sprintf('loading laser & magnetic data: << %s >>', iFolder);
    logMsg('info',msg,1,0);
    
    fixedData = imread(fullfile(iFolder, laserFileName));
    fileData = load(iFile);
    
    [transForm, refFrame] = get_image_tform(fixedData, movingData,...
        'checkPlot', kwargs.checkPlot, 'title', 'laser alignment');

    % transform the blank
    B111ferroTransformed = tform_data(blankData.B111ferro, transForm, refFrame);
    B111paraTransformed = tform_data(blankData.B111para, transForm, refFrame);

    % crop the FOV and subtract blank
    [x, y, w, h] = get_mask_extent(B111ferroTransformed);
    fileB111ferro = fileData.B111ferro(y:y+h, x:x+w);
    fileB111para = fileData.B111para(y:y+h, x:x+w);
    
    fileData.B111ferro = fileB111ferro- B111ferroTransformed(y:y+h, x:x+w);
    fileData.B111para = fileB111para - B111paraTransformed(y:y+h, x:x+w);
%     B111para = fileData.B111para - B111paraTransformed;
    subtractedData(iFolder) = fileData;
    
    if kwargs.save
        saveFilePath = fullfile(iFolder, 'B111BlankSub.mat');
        save(saveFilePath, '-struct', 'fileData');
    end

    if kwargs.checkPlot
        med = abs(median(fileB111ferro, 'all', 'omitnan'));
        st = std(fileB111ferro, [], 'all', 'omitnan');
        
        fig = figure('Units', 'normalized', ...
                 'Position',[0.1 0.1 0.8 0.8], 'Name', 'Blank Subtraction');
        
        sp1 = subplot(2,2,1);
        imagesc(fileB111ferro);
        title('Original');
        
        sp2 = subplot(2,2,2);
        imagesc(blankData.B111ferro);
        title('Blank');
        
        sp3 = subplot(2,2,3);
        blank = re_bin(blankData.B111ferro, fileData.B111ferro);
        imagesc(fileB111ferro - blank);
        title('Unaligned subtraction');
        
        sp4 = subplot(2,2,4);
        imagesc(fileData.B111ferro);
        title('Aligned subtraction');
        ax = [sp1, sp2, sp3, sp4];
        linkaxes(ax)
        
        for a = ax
           colorbar(a)
           set(ax, 'CLim', [-1, 1]*(med + 4 * st));
        end
    end
end
