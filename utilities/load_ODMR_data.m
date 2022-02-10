function runFiles = load_ODMR_data(folder)

runFiles = dir(fullfile(folder, 'run_*.mat'));

for i = 1:size(runFiles)
    data = load(fullfile(runFiles(i).folder, ... 
        runFiles(i).name), 'imgStack1', 'imgStack2', 'imgStack3', 'imgStack4', 'freqList');
    
    if isfield(data, 'imgStack3')
        runFiles(i).left = [data.imgStack1; data.imgStack2];
        runFiles(i).right = [data.imgStack3; data.imgStack4];
    else
        runFiles(i).left = data.imgStack1;
        runFiles(i).right = data.imgStack2;
    end
    %% frequencies in GHz
    runFiles(i).leftF = data.freqList(1:size(runFiles(i).left))/1e9;
    runFiles(i).rightF = data.freqList(size(runFiles(i).left)+1:end)/1e9;
end
end