% Custom validation function
function folderMustExist(folder)
    % Test for equal size
    if strcmp(folder, 'none')
        return
    end
        
    if ~exist(folder, 'dir')
        eid = 'function:inputError';
        msg = sprintf('folder << %s >> does not exist!', strrep(folder, '\', '\\'));
        throwAsCaller(MException(eid,msg))
    end
end
