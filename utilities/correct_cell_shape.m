function elements = correct_cell_shape(elements)
%[elements] = correct_cell_shape(elements)
% Takes a cell and converts its shape to (1,n) if the shape is (n,1)
if ischar(elements)
    elements = {elements};
end

% fix shape for nFolders
folderShape = size(elements);
if folderShape(1) > folderShape(2)
    elements = transpose(elements);
end

% remove training slash or backslash
for i = 1:size(elements,2)
    e = elements{i};
    if e(end) == filesep
        elements{i} = e(1:end-1);
    end
end