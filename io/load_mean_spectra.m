function [meanData freq struct] = load_mean_spectra(folder)
runFiles = dir(fullfile(folder, 'run_*.mat'));
headerFiles = dir(fullfile(folder, 'run_*header.txt'));
headerNames = {};
struct = [];

for i = 1:size(headerFiles)
   headerNames{end+1} = headerFiles(i).name; 
end

for i = 1:size(runFiles)
    runFiles(i).data = load(fullfile(runFiles(i).folder, ... 
        runFiles(i).name), 'disp1', 'disp2', 'freqList');
    runName = strsplit(runFiles(i).name,'.');
    runName = runName{1};
    idx = find(contains(headerNames, runName));
    
    headerFile = fullfile(headerFiles(idx).folder, headerFiles(idx).name);
    
    struct(i).header = read_header(headerFile);
    struct(i).header.headerFile = headerFile;
    
end

meanData = [];
freq = [];
for i = 1:size(runFiles)
    struct(i).left = runFiles(i).data.disp1;
    struct(i).right = runFiles(i).data.disp2;
    meanData = [meanData; [runFiles(i).data.disp1, runFiles(i).data.disp2]];
    
    % read frequencies, sometime LabView does not write them -> catch
    % exception
    try
        struct(i).leftF = runFiles(i).data.freqList(1:size(runFiles(i).data.disp1,2));
        struct(i).rightF = runFiles(i).data.freqList(size(runFiles(i).data.disp1,2)+1:end);
    catch
        struct(i).leftF = linspace(struct(i).header.startfreq1_Hz,struct(i).header.endfreq1_Hz, struct(i).header.numfreqs);
        struct(i).rightF = linspace(struct(i).header.startfreq2_Hz,struct(i).header.endfreq2_Hz, struct(i).header.numfreqs);
    end
    freq = [freq; [struct(i).leftF struct(i).rightF]];


end