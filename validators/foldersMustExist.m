% Custom validation function
function foldersMustExist(folders)
    % Test for equal size
    folders = correct_cell_shape(folders);
    for f = folders
        f = f{:};
        folderMustExist(f);
    end
end
