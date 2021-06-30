function data = reshape_QDM_data(data)
% reshapes data to freq:row:col
    data.imgStack1 = QDMreshape(data.imgStack1, expData.imgNumCols, expData.imgNumRows); % reshape into freq x col x rows
	data.imgStack2 = QDMreshape(data.imgStack2, expData.imgNumCols, expData.imgNumRows); % reshape into freq x col x rows
end

