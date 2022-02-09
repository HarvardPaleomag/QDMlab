function fileName = check_suffix(fileName)
%[fileName] = check_suffix(fileName)

if ~endsWith(fileName, '.mat')
    fileName = char([fileName, '.mat']);
end