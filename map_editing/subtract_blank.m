function subtractedData = subtract_blank(varargin)
%{

%}
close all

inParse = inputParser;
addParameter(inParse, 'nFolders', false);
addParameter(inParse, 'blankFile', false);
addParameter(inParse, 'checkPlot', false);
addParameter(inParse, 'save', true);

parse(inParse, varargin{:});

nFolders = inParse.Results.nFolders;
blankFile = inParse.Results.blankFile;
checkPlot = inParse.Results.checkPlot;

fileName = 'B111dataToPlot.mat';
laserFileName = 'laser.jpg';

%% manual
% if nFoilders and blankData uses default values i.e. false

if find(strcmp(inParse.UsingDefaults, 'nFolders'))
    %Load the correct "B111dataToPlot.mat" file
    f = helpdlg('Pick a B111 file');
    pause(2)

    if ishghandle(f)
        close(f)
    end

    [longfilename, pathname] = uigetfile('*.mat', 'Pick a B111 file');
    fullfilename = [pathname longfilename]
    nFolders = {pathname};
end

if find(strcmp(inParse.UsingDefaults, 'blankFile'))
    f = helpdlg('Pick a blank file');
    pause(2)

    if ishghandle(f)
        close(f)
    end

    [longfilenameBLANK, pathnameBLANK] = uigetfile('*.mat', 'Pick a Blank file');
    fullfilenameBLANK=[pathnameBLANK longfilenameBLANK];
    blankFile = fullfilenameBLANK;
end

%% automatic subtraction for all folders
% checks if none of the default arguments is used
nFolders = correct_cell_shape(nFolders);

disp(['<> loading blank file: <<' blankFile '>>'])
blankData = load(blankFile);
[nTransForms, nRefFrames] = get_tform_multi(blankFile, nFolders, ...
                            'reverse', true, ...
                            'laser',true, 'checkPlot', checkPlot);

for i = 1 : size(nFolders, 2)
    iFile = nFolders{i};
    iFile = [iFile filesep fileName];
    disp(['<> reading: << ' iFile ' >> and blankData'])

    [filepath,name,ext] = fileparts(iFile);
    fileData = load(iFile);

    % transform the blank
    B111ferroTransformed = tform_data(blankData.B111ferro, nTransForms(iFile), nRefFrames(iFile));
    B111paraTransformed = tform_data(blankData.B111para, nTransForms(iFile), nRefFrames(iFile));

    % crop the FOV and subtract blank
    [x, y, w, h] = get_mask_extent(B111ferroTransformed);
    fileB111ferro = fileData.B111ferro(y:y+h, x:x+w);
    fileB111para = fileData.B111para(y:y+h, x:x+w);
    B111ferro = fileB111ferro- B111ferroTransformed(y:y+h, x:x+w);
    B111para = fileB111para - B111paraTransformed(y:y+h, x:x+w);

%     FitCvgBoth=fileData.FitCvgBoth;
    ledImg=fileData.ledImg(y:y+h, x:x+w);
    negDiff=fileData.ledImg(y:y+h, x:x+w);
    posDiff=fileData.ledImg(y:y+h, x:x+w);

    if inParse.Results.save
        save([filepath 'B111BlankSub.mat'],'negDiff','posDiff', 'B111ferro', 'B111para', 'ledImg');
    end

    if inParse.Results.checkPlot
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
