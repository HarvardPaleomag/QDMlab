function data = QDMreshape(dataStack, imgNumCols, imgNumRows) 
% data = QDMreshape(dataStack, imgNumCols, imgNumRows) 
% reshapes the 51xnumber_of_pixel data into Y x X X freq

    data = reshape(dataStack, [], imgNumCols, imgNumRows); % reshape into freq x col x rows
    data = permute(data,[3 2 1]); % permute the axis to rows x cols x freq
end