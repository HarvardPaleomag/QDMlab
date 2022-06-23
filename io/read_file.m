function data = read_file(runFile)
%[data] = read_file('runFile')
arguments 
    runFile = 'none'
end

runFile = automatic_input_ui__(runFile, 'type', 'file', 'single', true, 'title', 'select run_00*.mat file');

if ~isstruct(runFile)
    runFile = dir(runFile);
end

dataFull = load(fullfile(runFile.folder, ... 
        runFile.name));

if isfield(dataFull, 'imgStack3')
    data.lowData = [dataFull.imgStack1; dataFull.imgStack2];
    data.highData = [dataFull.imgStack3; dataFull.imgStack4];
else
    data.lowData = dataFull.imgStack1;
    data.highData = dataFull.imgStack2;
end

% add mean spectra
data.lowMean = mean_without_nan(data.lowData);
data.highMean = mean_without_nan(data.highData);

data.numFreqs = dataFull.numFreqs;
data.imgNumCols = dataFull.imgNumCols;
data.imgNumRows = dataFull.imgNumRows;

%% frequencies in GHz
data.fLow = dataFull.freqList(1:size(data.lowData))'/1e9;
data.fHigh = dataFull.freqList(size(data.highData)+1:end)'/1e9;

msg = sprintf('File << %s >> in %s successfully loaded', runFile.name, runFile.folder);
logMsg('debug',msg,1,0);
end