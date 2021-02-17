function data = reshape_QDM_data(data)
% reshapes data to freq:row:col
    data.imgStack1 = QDMreshape(data.imgStack1, data.imgNumCols, data.imgNumRows);
	data.imgStack2 = QDMreshape(data.imgStack2, data.imgNumCols, data.imgNumRows);
end

