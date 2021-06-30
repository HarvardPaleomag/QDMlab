function elements = correct_cell_shape(elements)
% Takes a cell and converts its shape to (1,n) if the shape is (n,1)
if ischar(elements)
    elements = {elements};
end

% fix shape for nFolders
folderShape = size(elements);
if folderShape(1) > folderShape(2)
    elements = transpose(elements);
end