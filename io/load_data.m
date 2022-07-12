function runFiles = load_data(folder, polarity)
%[runFiles] = load_data('folder', 'polarity')
% loads the data from a specific folder
%
% Parameters
% ----------
%   folder:
%   polarity:
%   folder: ('none')
%   polarity: ('all')
% 
% Returns
% ----------
%   runFiles:

arguments
    folder = 'none'
    polarity = 'all'
end

folder = automatic_input_ui__(folder, 'type', 'dir', 'single', true);

runFiles = dir(fullfile(folder, 'run_*.mat'));
headerFiles = dir(fullfile(folder, 'run_*_header.txt'));

if isequal(polarity, 'all')
    idx = 1:size(runFiles);
else 
    idx = polarity;
end

for i = idx
    runFiles(i).data = read_file(runFiles(i));
    runFiles(i).header = read_header(fullfile(headerFiles(i).folder, headerFiles(i).name));
end

end
