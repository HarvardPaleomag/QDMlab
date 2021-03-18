% Custom validation function
function foldersMustExist(folders)
    % Test for equal size
%    folders = correct_cell_shape(folders); Mike: I don't think this
%    function works- it's transposing it even when the cell array is OK
    if ~all(cellfun(@exist, folders))
        eid = 'function:inputError';
        msg = sprintf('Some folders do not exist! Check cell\n');
        throwAsCaller(MException(eid,msg))
    end
end
