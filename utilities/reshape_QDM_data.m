function data = reshape_QDM_data(data)
% reshapes data to freq:row:col
    data.imgStack1 = reshape(data.imgStack1, [], expData.imgNumCols, expData.imgNumRows); % reshape into freq x col x rows
    data = permute(data.imgStack1,[3 2 1]); % permute the axis to rows x cols x freq

	data.imgStack2 = reshape(data.imgStack2, [], expData.imgNumCols, expData.imgNumRows); % reshape into freq x col x rows
    data.imgStack2 = permute(data.imgStack2,[3 2 1]); % permute the axis to rows x cols x freq
end

