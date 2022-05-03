function runFiles = load_ODMR_data(folder)
%[runFiles] = load_ODMR_data(folder)

runFiles = dir(fullfile(folder, 'run_*.mat'));

for i = 1:size(runFiles)
    try
        data = load(fullfile(runFiles(i).folder, ... 
            runFiles(i).name), 'imgStack1', 'imgStack2', 'imgStack3', 'imgStack4', 'freqList');
    catch
        data = load(fullfile(runFiles(i).folder, ... 
            runFiles(i).name), 'imgStack1', 'imgStack2', 'freqList');
    end

    if isfield(data, 'imgStack3')
        runFiles(i).left = [data.imgStack1; data.imgStack2];
        runFiles(i).right = [data.imgStack3; data.imgStack4];
    else
        runFiles(i).left = data.imgStack1;
        runFiles(i).right = data.imgStack2;
    end

    runFiles(i).leftMean = mean_without_nan(runFiles(i).left);
    runFiles(i).rightMean = mean_without_nan(runFiles(i).right);

    %% frequencies in GHz
    runFiles(i).leftF = data.freqList(1:size(runFiles(i).left))/1e9;
    runFiles(i).rightF = data.freqList(size(runFiles(i).left)+1:end)/1e9;
end
end