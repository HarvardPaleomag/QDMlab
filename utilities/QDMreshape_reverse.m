function data = QDMreshape_reverse(data, nFreq) 
% data = QDMreshape(dataStack, imgNumCols, imgNumRows) 
% reshapes the 51xnumber_of_pixel data into Y x X X freq
    data = permute(data,[3 2 1]); % permute the axis to rows x cols x freq
    data = reshape(data, nFreq, []); % reshape into freq x col x rows
end