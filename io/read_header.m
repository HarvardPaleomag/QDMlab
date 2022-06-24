function headerFile = read_header(headerFile)
%[headerFile] = read_header('headerFile')
%
% Parameters
% ----------
%   headerFile:
%   headerFile: ('none')
% 
% Returns
% ----------
%   headerFile:

arguments 
    headerFile = 'none'
end

headerFile = automatic_input_ui__(headerFile, 'type', 'file', 'single', true, 'title', 'select run_00*_header.txt file', 'filter', '*.txt');

if ~isstruct(headerFile)
    headerFile = dir(headerFile);
end


data = importdata(fullfile(headerFile.folder, headerFile.name));

for i=1:size(data.rowheaders)
    key = data.rowheaders{i};
    key = strrep(key, ' ', '_');
    key = strrep(key, '(', '');
    key = strrep(key, ')', '');
    key = strrep(key, '?', '');

    headerFile.(key) = data.data(i);
end

msg = sprintf('File << %s >> in %s successfully loaded', headerFile.name, headerFile.folder);
logMsg('debug',msg,1,0);

end