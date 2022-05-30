function runFiles = load_data(folder, nRes)
% loads the data from a specific folder
%

arguments
    folder = 'none'
    nRes = 'all'
end

folder = automatic_input_ui__(folder, 'type', 'dir', 'single', true);

runFiles = dir(fullfile(folder, 'run_*.mat'));
headerFiles = dir(fullfile(folder, 'run_*_header.txt'));

switch nRes
    case 'all'
        idx = 1:size(runFiles);
    case 1
        idx = 1;
    case 2
        idx = 2;
end

for i = idx
    runFiles(i).data = read_file(runFiles(i));
    runFiles(i).header = read_header(fullfile(headerFiles(i).folder, headerFiles(i).name));
end

msg = sprintf('%i file(s) in %s successfully loaded', size(idx,2), folder);
logMsg('debug',msg,1,0);
end
