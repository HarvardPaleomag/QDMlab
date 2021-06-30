function fileName = check_suffix(fileName)

if ~endsWith(fileName, '.mat')
    fileName = [fileName, '.mat'];
end