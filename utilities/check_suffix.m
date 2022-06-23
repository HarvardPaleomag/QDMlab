function fileName = check_suffix(fileName)
%[fileName] = check_suffix(fileName)
% adds '.mat' as a suffix if the file does not have one
if ~endsWith(fileName, '.mat')
    fileName = char([fileName, '.mat']);
end