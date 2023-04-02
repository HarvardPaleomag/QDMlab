function data = read_file(runFile)
% Reads a data file and returns the data in a structured format.
%
% Parameters
% ----------
%   runFile: str
%       Path to the data file or 'none' to select the file using a UI.
% 
% Returns
% ----------
%   data: struct
%       Struct containing the data read from the file.

arguments 
    runFile = 'none'
end

% If runFile is 'none', use UI to select the file
runFile = automatic_input_ui__(runFile, 'type', 'file', 'single', true, 'title', 'select run_00*.mat file');

% If runFile is not a struct, assume it is a path to the file and load the file metadata
if ~isstruct(runFile)
    runFile = dir(runFile);
end

% Load the data file
dataFull = load(fullfile(runFile.folder, runFile.name));

% If there are 4 image stacks, assume we have low and high frequency data
if isfield(dataFull, 'imgStack3')
    data.lowData = [dataFull.imgStack1; dataFull.imgStack2];
    data.highData = [dataFull.imgStack3; dataFull.imgStack4];
else
    data.lowData = dataFull.imgStack1;
    data.highData = dataFull.imgStack2;
end

% Add mean spectra
data.lowMean = mean_without_nan(data.lowData);
data.highMean = mean_without_nan(data.highData);

% Get the number of frequency points and the image dimensions
data.numFreqs = dataFull.numFreqs;
data.imgNumCols = dataFull.imgNumCols;
data.imgNumRows = dataFull.imgNumRows;

% Get the frequency axis (in GHz)
data.fLow = dataFull.freqList(1:size(data.lowData))'/1e9;
data.fHigh = dataFull.freqList(size(data.highData)+1:end)'/1e9;

% Log a message indicating that the file was successfully loaded
msg = sprintf('File << %s >> in %s successfully loaded', runFile.name, runFile.folder);
logMsg('debug',msg,1,0);
end
