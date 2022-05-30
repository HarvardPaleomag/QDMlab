function data = read_file(runFile)
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
    data.left = [dataFull.imgStack1; dataFull.imgStack2];
    data.right = [dataFull.imgStack3; dataFull.imgStack4];
else
    data.left = dataFull.imgStack1;
    data.right = dataFull.imgStack2;
end

data.leftMean = mean_without_nan(data.left);
data.rightMean = mean_without_nan(data.right);

%% frequencies in GHz
data.leftF = dataFull.freqList(1:size(data.left))/1e9;
data.rightF = dataFull.freqList(size(data.left)+1:end)/1e9;

msg = sprintf('File << %s >> in %s successfully loaded', runFile.name, runFile.folder);
logMsg('debug',msg,1,0);
end