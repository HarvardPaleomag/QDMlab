% Custom validation function
function foldersMustExist(folders)
    % Test for equal size
    folders = correct_cell_shape(folders);
    if ~all(cellfun(@exist, folders))
        eid = 'function:inputError';
        msg = sprintf('Some folders do not exist! Check cell\n');
        throwAsCaller(MException(eid,msg))
    end
end
