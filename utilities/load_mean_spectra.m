function [meanData freq struct] = load_mean_spectra(folder)


runFiles = dir(fullfile(folder, 'run_*.mat'));

for i = 1:size(runFiles)
    runFiles(i).data = load(fullfile(runFiles(i).folder, ... 
        runFiles(i).name), 'disp1', 'disp2', 'freqList');
end
struct = [];
meanData = [];
freq = [];
for i = 1:size(runFiles)
    struct(i).left = runFiles(i).data.disp1;
    struct(i).right = runFiles(i).data.disp2;
    struct(i).leftF = runFiles(i).data.freqList(1:size(runFiles(i).data.disp1,2));
    struct(i).rightF = runFiles(i).data.freqList(size(runFiles(i).data.disp1,2)+1:end);
    
    meanData = [meanData; [runFiles(i).data.disp1, runFiles(i).data.disp2]];
    freq = [freq; runFiles(i).data.freqList];
end