function nFolders = correct_cell_shape(nFolders)
% Takes a cell and converts its shape to (1,n) if the shape is (n,1)

% fix shape for nFolders
folderShape = size(nFolders);
if folderShape(1) > folderShape(2)
    nFolders = transpose(nFolders);
end