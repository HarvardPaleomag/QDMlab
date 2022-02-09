function [meanData freq] = load_mean_spectra(folder)


runFiles = dir(fullfile(folder, 'run_*.mat'));

for i = 1:size(runFiles)
    runFiles(i).data = load(fullfile(runFiles(i).folder, ... 
        runFiles(i).name), 'disp1', 'disp2', 'freqList');
end

meanData = [];
freq = [];
for i = 1:size(runFiles)
    meanData = [meanData; [runFiles(i).data.disp1, runFiles(i).data.disp2]];
    freq = [freq; runFiles(i).data.freqList];
end