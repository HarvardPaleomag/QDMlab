% newQDM_makeB111file.m -- MATLAB script for taking B111 diff and sum files
% and combining with laser/led images to make traditional B111toplot.mat

% Open a dialog box for the user to pick a folder
folderPath = uigetdir(pwd, 'Select folder containing QDM data');
bintype = 1;

%%%%%%%%%%%%%%%%%%%% Read the map data files

% If the user did not cancel
if folderPath ~= 0
    % Build full path to B111_diff.mat
    matFilediff = fullfile(folderPath, 'B111_diff.mat');
    matFilesum = fullfile(folderPath, 'B111_sum.mat');
    matFiles1 = fullfile(folderPath, 'B111_s1.mat');
    matFiles2 = fullfile(folderPath, 'B111_s2.mat');
    snapshot = fullfile(folderPath, 'snapshot.csv');
    
    % Check that the file exists
    if exist(matFilediff, 'file')
        % Load only the variable B111_diff from the MAT-file
        S = load(matFilediff, 'B111_diff');
        % Create new variable B111_ferro from B111_diff
        B111ferro = S.B111_diff;
    else
        error('File B111_diff.mat not found in the selected folder.');
    end
    % Check that the file exists
    if exist(matFilesum, 'file')
        % Load only the variable B111_diff from the MAT-file
        S = load(matFilesum, 'B111_sum');
        % Create new variable B111_ferro from B111_diff
        B111para = S.B111_sum;
    else
        error('File B111_sum.mat not found in the selected folder.');
    end
    % Check that the file exists
    if exist(matFiles1, 'file')
        % Load only the variable B111_diff from the MAT-file
        S = load(matFiles1, 'B111_s1');
        % Create new variable B111_ferro from B111_diff
        B111_s1 = S.B111_s1;
    else
        error('File B111_s1.mat not found in the selected folder.');
    end
    % Check that the file exists
    if exist(matFiles2, 'file')
        % Load only the variable B111_diff from the MAT-file
        S = load(matFiles2, 'B111_s2');
        % Create new variable B111_ferro from B111_diff
        B111_s2 = S.B111_s2;
    else
        error('File B111_s2.mat not found in the selected folder.');
    end
    % Check that the file exists
    if exist(snapshot, 'file')
        % Load only the variable B111_diff from the MAT-file
        ledImg = readmatrix(snapshot);
    else
        error('File B111_s2.mat not found in the selected folder.');
    end
else
    disp('Folder selection canceled by user.');
end


%%%%%%%%%%%%%%%%%%%% Read the yaml file
yamlFile = fullfile(folderPath, 'parameters.yaml');;   % YAML file

% Read all lines (R2020b+). For older versions, use fopen/fgetl loop.
lines = readlines(yamlFile);

S = struct();  % output struct

for i = 1:numel(lines)
    line = strtrim(lines(i));

    % Skip empty or comment lines
    if line == "" || startsWith(line, "#")
        continue;
    end

    % Find first colon; skip line if none
    idx = strfind(line, ':');
    if isempty(idx)
        continue;   % << skip line with no ":"
    end

    % Split into key and value
    key    = strtrim(extractBefore(line, idx(1)));
    valStr = strtrim(extractAfter(line, idx(1)));

    % Convert value string to numeric
    valNum = str2double(valStr);
    if isnan(valNum)
        % If not numeric, skip (or handle as you wish)
        continue;
    end

    % Ensure valid MATLAB field name
    key = matlab.lang.makeValidName(key);

    % Assign to struct
    S.(key) = valNum;
end

nsweeps = S.TotalNumberOfSweeps;
timepersweep = S.NumberOfFrequencyPointsPerBlock * S.NumberOfBlocks / S.FrameRate_Hz_;
totaltime = nsweeps * timepersweep;


%%%%%%%%%%%%%%%%%%%% Cropping and saving
if ismac
    sep = '/';
elseif ispc
    sep = '\';
else
    % Fallback for other systems (Linux, etc.)
    sep = '/';
end

%figure out binning
if size(B111ferro,2) == 1936
    bintype = 1;
elseif size(B111ferro,2) == 967
    bintype = 2;
elseif size(B111ferro,2) == 483
    bintype = 3;
end

% crop to remove the weird hot pixels in the corner of the new sensor
if bintype == 1
    B111ferro = B111ferro(1:1200, 1:1920);
    B111para = B111para(1:1200, 1:1920);
    B111_s1 = B111_s1(1:1200, 1:1920);
    B111_s2 = B111_s2(1:1200, 1:1920);
    ledImg = ledImg(1:1200, 1:1920);
elseif bintype == 2
    B111ferro = B111ferro(1:600, 1:960);
    B111para = B111para(1:600, 1:960);
    B111_s1 = B111_s1(1:600, 1:960);
    B111_s2 = B111_s2(1:600, 1:960);
    ledImg = ledImg(1:1200, 1:1920);
elseif bintype == 3
    B111ferro = B111ferro(1:300, 1:480);
    B111para = B111para(1:300, 1:480);
    B111_s1 = B111_s1(1:300, 1:480);
    B111_s2 = B111_s2(1:300, 1:480);
    ledImg = ledImg(1:1200, 1:1920);
end

% flip up-down to maintain sample north
B111ferro = flipud(B111ferro);
B111para = flipud(B111para);
B111_s1 = flipud(B111_s1);
B111_s2 = flipud(B111_s2);
ledImg = flipud(ledImg);

save([folderPath sep 'B111data.mat'], 'B111ferro', 'B111para', 'B111_s1','B111_s2','ledImg','nsweeps','timepersweep','totaltime');








