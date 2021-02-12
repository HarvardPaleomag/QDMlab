% Custom validation function
function fileMustExistInFolder(fileName, folders)
    fileName = check_suffix(fileName);
    
    % construct folder
    folders = cellfun(@(S) fullfile(S, fileName), folders, 'Uniform', 0);
    
    % Test for equal size
    if ~all(cellfun(@exist, folders))
        eid = 'function:inputError';
        msg = sprintf('At least one folder does not contain << %s >>! Check data\n', fileName);
        throwAsCaller(MException(eid,msg))
    end
end
